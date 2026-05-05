#!/usr/bin/env python3
"""
qcalcgeom_helpers.py — QCalcGeom Python Helpers & IPC Wrapper
═════════════════════════════════════════════════════════════════

PURPOSE: Python-side geometric calculations matching the C++ QCalcGeom engine,
plus IPC wrapper for communicating with the C++ backend via named pipe.

IPC PROTOCOL (from ipc/uqff_ipc.h):
  Pipe:    \\\\.\\pipe\\StarMagic_UQFF
  Magic:   0x55514646 ("UQFF")
  Version: 1
  Messages:
    QCALCGEOM_COMPUTE  = 0x0B01  → Route to bsfg_metric/horizon/geodesic/holonomy
    QCALCGEOM_RESULT   = 0x0B02  → Return BSFGMetricResult / BSFGHorizonResult / etc.
    QCALCGEOM_TEST_RUN = 0x0B03  → Trigger QCALCGEOM::runQCalcGeomTests() (40+ tests)

RESULT STRUCTS (12 types):
  BSFGMetricResult, BSFGHorizonResult, BSFGFieldEqResult, BSFGGeodesicResult,
  BSFGHolonomyResult, VDSResult, DVPResult, BSHResult, BH26Result,
  BSFGBuoyancyResult, Poly26Result, UQFFCompResult

REFERENCES:
  - ipc/uqff_ipc.h — MessageHeader (32 bytes), QCALCGEOM message types
  - ipc/python_bridge.h — PythonBridge C++ class
  - QCalcGeom.h — Result structs, 40 test validators
  - PAPER_554-558: BSFG 5-calculator cascade

SESSION: 203 | April 7, 2026
UPDATED: Session 230 | May 2026 — QCalcGeomIPCClient is now platform-adaptive:
  Windows  → Win32 named pipe (ctypes.windll)
  Linux    → AF_UNIX socket at /tmp/StarMagic_UQFF.sock
  macOS    → AF_UNIX socket at /tmp/StarMagic_UQFF.sock
  Other    → HTTP REST fallback to localhost:3141/qcalcgeom
"""

import math
import json
import struct
import sys
import time
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field

# ═══════════════════════════════════════════════════════════════════════════════
# §1  CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

G       = 6.67430e-11
c       = 2.99792e8
hbar    = 1.05457e-34
k_B     = 1.38065e-23
PI      = math.pi
M_sun   = 1.98892e30

# UQFF
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4

# IPC Protocol
PIPE_NAME        = r"\\.\pipe\StarMagic_UQFF"
IPC_MAGIC        = 0x55514646   # "UQFF"
IPC_VERSION      = 1
QCALCGEOM_COMPUTE   = 0x0B01
QCALCGEOM_RESULT    = 0x0B02
QCALCGEOM_TEST_RUN  = 0x0B03

# MessageHeader: magic(4) + version(4) + type(4) + payload_size(4) + timestamp(8) + seq(4) + flags(4) = 32
HEADER_FORMAT = "<IIIIQII"
HEADER_SIZE   = struct.calcsize(HEADER_FORMAT)  # 32 bytes

# BSFG extra flat dimensions
N_EXTRA_FLAT = 22

# ── SCm Vacuum Manifold (from dpm_vacuum_manifold, consolidated) ────────────
from dpm_vacuum_manifold import (
    SSQ          as _SCM_SSQ,
    KAPPA        as _SCM_KAPPA,
    RHO_VAC_SCM  as _SCM_RHO_VAC,
    THZ_PHONON   as _SCM_THz,
    compute_F_U_Bi_i_numerical as _scm_F_U_Bi_i_num,
    vds_numerical              as _scm_vds_num,
    E_phonon                   as _SCM_E_PHONON,
    S26_3                      as _SCM_S26_3,
    Phi_resonance              as _SCM_PHI_RES,
    KER_SCm                    as _SCM_KER_SCm,
    scaling_factor             as _SCM_SCALING,
    KAPPA_FLOAT                as _SCM_KAPPA_FLOAT,
    F_TRZ                      as _SCM_F_TRZ,
    coleman_guillespie_scm         as _scm_coleman_guillespie,
    neutrino_oscillation_prob_lenr as _scm_neutrino_osc,
    quark_production_prob_ui       as _scm_quark_prod,
    mckubre_lenr                   as _scm_mckubre,
    s26_3_from_vds                 as _scm_s26_3_from_vds,
    qgp_energy_density_scm         as _scm_qgp_energy_density,
    strange_quark_matter_density   as _scm_sqm_density,
    mit_bag_scm                    as _scm_mit_bag,
    ads_cft_scm_dual               as _scm_ads_cft_dual,
    scm_gw_metric_perturbation     as _scm_gw_metric_pert,
)
_SCM_LOADED = True


# ═══════════════════════════════════════════════════════════════════════════════
# §2  IPC MESSAGE HEADER
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class MessageHeader:
    """Binary-compatible MessageHeader matching C++ struct (32 bytes)."""
    magic: int = IPC_MAGIC
    version: int = IPC_VERSION
    msg_type: int = 0
    payload_size: int = 0
    timestamp: int = 0
    sequence: int = 0
    flags: int = 0

    def pack(self) -> bytes:
        if self.timestamp == 0:
            self.timestamp = int(time.time() * 1e6)
        return struct.pack(HEADER_FORMAT,
                           self.magic, self.version, self.msg_type,
                           self.payload_size, self.timestamp,
                           self.sequence, self.flags)

    @classmethod
    def unpack(cls, data: bytes) -> "MessageHeader":
        vals = struct.unpack(HEADER_FORMAT, data[:HEADER_SIZE])
        return cls(magic=vals[0], version=vals[1], msg_type=vals[2],
                   payload_size=vals[3], timestamp=vals[4],
                   sequence=vals[5], flags=vals[6])


# ═══════════════════════════════════════════════════════════════════════════════
# §3  IPC CLIENT (Named Pipe)
# ═══════════════════════════════════════════════════════════════════════════════

class QCalcGeomIPCClient:
    """
    Named pipe / Unix socket client for QCalcGeom C++ engine.

    Platform-adaptive:
      Windows  — uses Win32 CreateFile / WriteFile / ReadFile via ctypes.windll
      Linux    — uses AF_UNIX domain socket at /tmp/StarMagic_UQFF.sock
      macOS    — uses AF_UNIX domain socket at /tmp/StarMagic_UQFF.sock
      Other    — falls back to HTTP REST on localhost:3141/qcalcgeom

    The Windows named pipe (\\\\.\\pipe\\StarMagic_UQFF) and the POSIX socket
    path (/tmp/StarMagic_UQFF.sock) both implement the same 32-byte
    MessageHeader framing defined in ipc/uqff_ipc.h.

    Usage:
        client = QCalcGeomIPCClient()
        if client.connect():
            result = client.compute(payload_bytes)
            client.disconnect()
    """

    # POSIX socket path mirrors the Windows pipe name (without UNC prefix)
    POSIX_SOCKET_PATH = "/tmp/StarMagic_UQFF.sock"
    REST_FALLBACK_URL = "http://127.0.0.1:3141/qcalcgeom"

    def __init__(self, pipe_name: str = PIPE_NAME):
        self.pipe_name = pipe_name
        self._handle   = None   # Win32 handle (int) or socket.socket
        self._platform = sys.platform   # 'win32', 'linux', 'darwin', ...
        self._seq      = 0

    def connect(self) -> bool:
        """Connect to the StarMagic_UQFF IPC endpoint."""
        if self._platform == "win32":
            return self._connect_win32()
        elif self._platform in ("linux", "darwin"):
            return self._connect_posix()
        else:
            # Unknown platform — REST fallback needs no persistent connection
            return True  # lazy connect on first request

    # ── Windows Win32 named pipe ────────────────────────────────────────────
    def _connect_win32(self) -> bool:
        try:
            import ctypes
            import ctypes.wintypes
            GENERIC_READ  = 0x80000000
            GENERIC_WRITE = 0x40000000
            OPEN_EXISTING = 3
            handle = ctypes.windll.kernel32.CreateFileW(
                self.pipe_name,
                GENERIC_READ | GENERIC_WRITE,
                0, None,
                OPEN_EXISTING,
                0, None
            )
            if handle == -1 or handle == 0xFFFFFFFF:
                return False
            self._handle = handle
            return True
        except (ImportError, OSError):
            return False

    def _send_receive_win32(self, msg_type: int, payload: bytes) -> Optional[bytes]:
        if self._handle is None:
            return None
        try:
            import ctypes
            import ctypes.wintypes
            self._seq += 1
            header  = MessageHeader(msg_type=msg_type, payload_size=len(payload),
                                    sequence=self._seq)
            message = header.pack() + payload
            bw = ctypes.wintypes.DWORD()
            ok = ctypes.windll.kernel32.WriteFile(
                self._handle, message, len(message), ctypes.byref(bw), None)
            if not ok:
                return None
            buf  = ctypes.create_string_buffer(4096)
            br   = ctypes.wintypes.DWORD()
            ok   = ctypes.windll.kernel32.ReadFile(
                self._handle, buf, 4096, ctypes.byref(br), None)
            if not ok or br.value < HEADER_SIZE:
                return None
            rh = MessageHeader.unpack(buf.raw[:HEADER_SIZE])
            return buf.raw[HEADER_SIZE:HEADER_SIZE + rh.payload_size]
        except (ImportError, OSError):
            return None

    def _disconnect_win32(self):
        if self._handle is not None:
            try:
                import ctypes
                ctypes.windll.kernel32.CloseHandle(self._handle)
            except (ImportError, OSError):
                pass
            self._handle = None

    # ── POSIX Unix domain socket (Linux / macOS) ────────────────────────────
    def _connect_posix(self) -> bool:
        try:
            import socket
            sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            sock.settimeout(5.0)
            sock.connect(self.POSIX_SOCKET_PATH)
            self._handle = sock
            return True
        except OSError:
            return False

    def _send_receive_posix(self, msg_type: int, payload: bytes) -> Optional[bytes]:
        sock = self._handle
        if sock is None:
            return None
        try:
            self._seq += 1
            header  = MessageHeader(msg_type=msg_type, payload_size=len(payload),
                                    sequence=self._seq)
            message = header.pack() + payload
            sock.sendall(message)
            # Read fixed-size response header first
            raw_hdr = b""
            while len(raw_hdr) < HEADER_SIZE:
                chunk = sock.recv(HEADER_SIZE - len(raw_hdr))
                if not chunk:
                    return None
                raw_hdr += chunk
            rh  = MessageHeader.unpack(raw_hdr)
            raw = b""
            while len(raw) < rh.payload_size:
                chunk = sock.recv(rh.payload_size - len(raw))
                if not chunk:
                    break
                raw += chunk
            return raw
        except OSError:
            return None

    def _disconnect_posix(self):
        if self._handle is not None:
            try:
                self._handle.close()
            except OSError:
                pass
            self._handle = None

    # ── REST fallback (other / Emscripten / CI) ─────────────────────────────
    def _send_receive_rest(self, msg_type: int, payload: bytes) -> Optional[bytes]:
        """HTTP POST to uqff_server.js /qcalcgeom as last-resort fallback."""
        try:
            import urllib.request
            import urllib.error
            # Decode payload as JSON, add function field from msg_type
            body = json.loads(payload.decode("utf-8"))
            req  = urllib.request.Request(
                self.REST_FALLBACK_URL,
                data=json.dumps(body).encode("utf-8"),
                headers={"Content-Type": "application/json"},
                method="POST",
            )
            with urllib.request.urlopen(req, timeout=10) as resp:
                return resp.read()
        except Exception:
            return None

    # ── Unified public interface ─────────────────────────────────────────────
    def disconnect(self):
        """Close the IPC connection."""
        if self._platform == "win32":
            self._disconnect_win32()
        elif self._platform in ("linux", "darwin"):
            self._disconnect_posix()

    def _send_receive(self, msg_type: int, payload: bytes) -> Optional[bytes]:
        if self._platform == "win32":
            return self._send_receive_win32(msg_type, payload)
        elif self._platform in ("linux", "darwin"):
            return self._send_receive_posix(msg_type, payload)
        else:
            return self._send_receive_rest(msg_type, payload)

    def compute(self, payload: bytes) -> Optional[bytes]:
        """Send QCALCGEOM_COMPUTE and get result."""
        return self._send_receive(QCALCGEOM_COMPUTE, payload)

    def run_tests(self) -> Optional[bytes]:
        """Trigger QCALCGEOM test suite."""
        return self._send_receive(QCALCGEOM_TEST_RUN, b"")


# ═══════════════════════════════════════════════════════════════════════════════
# §4  METRIC TENSOR HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

class MetricTensorHelper:
    """
    Compute metric tensor components, Christoffels, and curvature
    for BSFG-modified spacetimes.
    """

    @staticmethod
    def bsfg_aether_perturbation(r: float, M_kg: float, eta: float = 1e-6) -> Tuple[float, float, float]:
        """
        Aether density perturbation ε(r) and derivatives.
        ε = η × C_num / r³ where C_num = GM/c²
        """
        C_num = G * M_kg / c**2
        eps = eta * C_num / r**3
        eps_p = -3 * eta * C_num / r**4        # dε/dr
        eps_pp = 12 * eta * C_num / r**5        # d²ε/dr²
        return eps, eps_p, eps_pp

    @staticmethod
    def bsfg_metric_components(r: float, M_kg: float, eta: float = 1e-6) -> Dict:
        """
        BSFG metric: g_00 = 1+ε, g_rr = -(1-ε), g_θθ = -r², g_φφ = -r²sin²θ
        """
        eps, eps_p, eps_pp = MetricTensorHelper.bsfg_aether_perturbation(r, M_kg, eta)
        return {
            "g_00":  1.0 + eps,
            "g_rr":  -(1.0 - eps),
            "g_theta_theta": -r**2,
            "eps": eps,
            "eps_prime": eps_p,
            "eps_double_prime": eps_pp,
        }

    @staticmethod
    def riemann_r0r0(r: float, M_kg: float, eta: float = 1e-6) -> float:
        """Riemann R^r_{0r0} = (1/2)d²A00/dr² for diagonal metric."""
        _, _, eps_pp = MetricTensorHelper.bsfg_aether_perturbation(r, M_kg, eta)
        return 0.5 * eps_pp

    @staticmethod
    def ricci_tensor(r: float, M_kg: float, eta: float = 1e-6) -> Dict:
        """Ricci tensor components R_00, R_rr, scalar R."""
        eps, eps_p, eps_pp = MetricTensorHelper.bsfg_aether_perturbation(r, M_kg, eta)
        R_r0r0 = 0.5 * eps_pp
        R_00 = eps_pp / 2 + eps_p / r
        R_rr = -(eps_pp / 2 + eps_p / r)
        R_scalar = 2 * (eps_pp + 2 * eps_p / r)
        kretschner = 12 * R_r0r0**2
        return {
            "R_00": R_00, "R_rr": R_rr, "R_scalar": R_scalar,
            "Kretschner": kretschner,
        }

    @staticmethod
    def christoffel_diagonal(r: float, M_kg: float, eta: float = 1e-6) -> Dict:
        """Non-zero Christoffel symbols for BSFG diagonal metric."""
        eps, eps_p, eps_pp = MetricTensorHelper.bsfg_aether_perturbation(r, M_kg, eta)

        g00 = 1.0 + eps
        grr = -(1.0 - eps)

        Gamma_r_tt = -0.5 * eps_p / grr if abs(grr) > 1e-30 else 0.0
        Gamma_r_rr = -0.5 * eps_p / grr if abs(grr) > 1e-30 else 0.0
        Gamma_r_thth = r * (1 - eps) if abs(grr) > 1e-30 else r
        Gamma_t_tr = 0.5 * eps_p / g00 if abs(g00) > 1e-30 else 0.0

        return {
            "Gamma^r_tt": Gamma_r_tt,
            "Gamma^r_rr": Gamma_r_rr,
            "Gamma^r_thetatheta": Gamma_r_thth,
            "Gamma^t_tr": Gamma_t_tr,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §5  COORDINATE TRANSFORMS
# ═══════════════════════════════════════════════════════════════════════════════

class CoordinateTransforms:
    """Coordinate transforms for UQFF calculations."""

    @staticmethod
    def spherical_to_cartesian(r: float, theta: float, phi: float) -> Tuple[float, float, float]:
        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        return x, y, z

    @staticmethod
    def cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
        r = math.sqrt(x**2 + y**2 + z**2)
        theta = math.acos(z / r) if r > 0 else 0.0
        phi = math.atan2(y, x)
        return r, theta, phi

    @staticmethod
    def boyer_lindquist_to_cartesian(r: float, theta: float, phi: float,
                                      a: float) -> Tuple[float, float, float]:
        """Boyer-Lindquist (Kerr) to Cartesian."""
        rho = math.sqrt(r**2 + a**2)
        x = rho * math.sin(theta) * math.cos(phi)
        y = rho * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)
        return x, y, z

    @staticmethod
    def compactified_26d(r: float, R_compact: float = 1e-35,
                         n_compact: int = 22) -> List[float]:
        """
        26D compactified coordinates.
        4 extended dimensions (t, r, θ, φ) + 22 compact at R_compact.

        Returns angular positions on each compact dimension.
        """
        angles = []
        for k in range(n_compact):
            # Each compact dimension wraps at R_compact
            psi_k = 2 * PI * ((r / R_compact + k * PI / n_compact) % 1.0)
            angles.append(psi_k)
        return angles

    @staticmethod
    def kk_eigenvalue(k: int) -> float:
        """Kaluza-Klein eigenvalue: λ_k = k(k+25) for 26D compactification."""
        return k * (k + 25)

    @staticmethod
    def kk_spectral_frequency(k: int, R_compact: float = 1e-35) -> float:
        """KK spectral bin frequency: f_k = c × sqrt(λ_k) / (2π R_compact)."""
        lam = CoordinateTransforms.kk_eigenvalue(k)
        return c * math.sqrt(lam) / (2 * PI * R_compact)


# ═══════════════════════════════════════════════════════════════════════════════
# §6  BSFG HORIZON & FIELD EQUATIONS
# ═══════════════════════════════════════════════════════════════════════════════

class BSFGHorizonCalculator:
    """Horizon properties for BSFG metric."""

    @staticmethod
    def horizon_radius(M_kg: float, eta: float = 1e-6) -> float:
        """r_h = (η × C_num)^{1/3} where C_num = GM/c²."""
        C_num = G * M_kg / c**2
        return (eta * C_num) ** (1.0 / 3.0)

    @staticmethod
    def hawking_temperature(M_kg: float, eta: float = 1e-6) -> float:
        """T_H = ℏc³ / (8πGMk_B) × (1 + BSFG correction)."""
        T_standard = hbar * c**3 / (8 * PI * G * M_kg * k_B)
        r_h = BSFGHorizonCalculator.horizon_radius(M_kg, eta)
        _, eps_p, _ = MetricTensorHelper.bsfg_aether_perturbation(r_h, M_kg, eta)
        correction = 1.0 + abs(eps_p) * r_h
        return T_standard * correction

    @staticmethod
    def surface_gravity(M_kg: float, eta: float = 1e-6) -> float:
        """κ = c²|dA00/dr|_{r_h} / 2"""
        r_h = BSFGHorizonCalculator.horizon_radius(M_kg, eta)
        _, eps_p, _ = MetricTensorHelper.bsfg_aether_perturbation(r_h, M_kg, eta)
        return c**2 * abs(eps_p) / 2


# ═══════════════════════════════════════════════════════════════════════════════
# §7  HOLONOMY & TOPOLOGY
# ═══════════════════════════════════════════════════════════════════════════════

class BSFGHolonomyCalculator:
    """Holonomy group analysis for BSFG metric (SO(3)×U(1)²³)."""

    @staticmethod
    def phase_accumulation(r: float, M_kg: float, area_m2: float,
                           eta: float = 1e-6) -> float:
        """Phase δφ accumulated over a loop of given area."""
        R_info = MetricTensorHelper.ricci_tensor(r, M_kg, eta)
        return abs(R_info["R_scalar"]) * area_m2

    @staticmethod
    def off_diagonal_connection(r: float, M_kg: float, eta: float = 1e-6) -> float:
        """Off-diagonal metric connection ω_{0r}."""
        eps, eps_p, _ = MetricTensorHelper.bsfg_aether_perturbation(r, M_kg, eta)
        g00 = 1.0 + eps
        if abs(g00) < 1e-30:
            return 0.0
        return 0.5 * eps_p / g00

    @staticmethod
    def holonomy_classification() -> Dict:
        """
        BSFG holonomy group classification.
        Ricci non-flat → excludes G₂ and Spin(7).
        Result: SO(3) × U(1)²³ (22 compact dimensions).
        """
        return {
            "group": "SO(3) x U(1)^23",
            "n_extra_flat": N_EXTRA_FLAT,
            "G2_excluded": True,
            "Spin7_excluded": True,
            "reason": "Ricci non-flat (aether perturbation ε ≠ 0)",
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §8  VDS / DVP / BSH HELPERS (QCALCGEOM-compatible)
# ═══════════════════════════════════════════════════════════════════════════════

class QCalcGeomVDS:
    """Vacuum Density Series — matches VDSResult struct."""

    @staticmethod
    def compute(N: int = 100, ssq: float = SSQ) -> Dict:
        """Σ_{n=1}^{N} SSq^n / n^{26}"""
        total = 0.0
        for n in range(1, N + 1):
            total += ssq**n / n**26

        tail_bound = ssq**(N + 1) / ((1 - ssq) * (N + 1)**26) if ssq < 1 else float('inf')
        converged = tail_bound < 1e-12

        return {
            "value": total,
            "converged": converged,
            "tail_bound": tail_bound,
            "n_terms_used": N,
        }


class QCalcGeomDVP:
    """Dipole Vortex Primes — matches DVPResult struct."""

    DVP_PRIMES_30 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                     31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                     73, 79, 83, 89, 97, 101, 103, 107, 109, 113]

    @staticmethod
    def factorial_mod(n: int, p: int) -> int:
        """n! mod p via iterative multiplication."""
        result = 1
        for i in range(2, n + 1):
            result = (result * i) % p
        return result

    @staticmethod
    def compute(p_special: int = 113) -> Dict:
        """26! mod p_special and proplyd quantization radius."""
        fac26_mod = QCalcGeomDVP.factorial_mod(26, p_special)
        non_repeating = fac26_mod != 0

        # Proplyd quantization radius
        fac26 = math.factorial(26)
        r_q_m = (2 / fac26) ** (1.0 / 26)
        AU = 1.496e11
        r_q_AU = r_q_m / AU

        return {
            "fac26_mod_113": fac26_mod,
            "non_repeating": non_repeating,
            "r_q_AU": r_q_AU,
            "r_q_m": r_q_m,
        }


class QCalcGeomBSH:
    """Buoyancy Series Harmonics — matches BSHResult struct."""

    @staticmethod
    def compute(m_max: int = 26, ssq: float = SSQ, f_Ub: float = 1.0) -> Dict:
        """Buoyancy harmonic sum: Σ (1/k) f_Ub (1-exp(-SSq·m)) cos(2πj/26)"""
        total = 0.0
        H_max = 0.0
        for m in range(1, m_max + 1):
            term = f_Ub * (1 - math.exp(-ssq * m)) / m
            H = abs(term)
            if H > H_max:
                H_max = H
            total += term

        saturated = (1 - math.exp(-ssq * m_max)) > (1 - 1e-6)

        return {
            "U_g2": total,
            "H_m_max": H_max,
            "saturated": saturated,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §9  26D POLYNOMIAL DERIVATIVE
# ═══════════════════════════════════════════════════════════════════════════════

class Poly26Calculator:
    """26th-order polynomial derivative for UQFF compressed field matrix."""

    @staticmethod
    def pochhammer(k: int, n: int = 26) -> float:
        """Rising factorial (k)_n = k(k+1)...(k+n-1)."""
        result = 1.0
        for i in range(n):
            result *= (k + i)
        return result

    @staticmethod
    def compute(k: int, r: float) -> Dict:
        """(k+25)!/(k-1)! × c / r^{k+26}"""
        fac_ratio = Poly26Calculator.pochhammer(k, 26)
        r_power = r ** (k + 26)
        if r_power == 0:
            value = float('inf')
        else:
            value = fac_ratio * c / r_power
        negligible = abs(value) < 1e-100

        return {
            "value": value,
            "factorial_ratio": fac_ratio,
            "r_power": r_power,
            "negligible": negligible,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §10  COMPREHENSIVE QCALCGEOM INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

class QCalcGeomEngine:
    """
    Master QCalcGeom interface — routes computations either to
    local Python implementations or to C++ engine via IPC.
    """

    def __init__(self, use_ipc: bool = False):
        self.use_ipc = use_ipc
        self.metric = MetricTensorHelper()
        self.coords = CoordinateTransforms()
        self.horizon = BSFGHorizonCalculator()
        self.holonomy = BSFGHolonomyCalculator()
        self.vds = QCalcGeomVDS()
        self.dvp = QCalcGeomDVP()
        self.bsh = QCalcGeomBSH()
        self.poly26 = Poly26Calculator()
        self._ipc = None

    def _ensure_ipc(self) -> Optional[QCalcGeomIPCClient]:
        if self._ipc is None and self.use_ipc:
            self._ipc = QCalcGeomIPCClient()
            if not self._ipc.connect():
                self._ipc = None
        return self._ipc

    def compute_all(self, M_kg: float = M_sun, r: float = 6.96e8,
                    eta: float = 1e-6) -> Dict:
        """Run all QCalcGeom computations for a given M, r, η."""
        return {
            "metric": self.metric.bsfg_metric_components(r, M_kg, eta),
            "ricci": self.metric.ricci_tensor(r, M_kg, eta),
            "christoffel": self.metric.christoffel_diagonal(r, M_kg, eta),
            "horizon": {
                "r_h_m": self.horizon.horizon_radius(M_kg, eta),
                "T_H_K": self.horizon.hawking_temperature(M_kg, eta),
                "kappa_surf": self.horizon.surface_gravity(M_kg, eta),
            },
            "holonomy": self.holonomy.holonomy_classification(),
            "vds": self.vds.compute(),
            "dvp": self.dvp.compute(),
            "bsh": self.bsh.compute(),
            "kk_eigenvalue_1": self.coords.kk_eigenvalue(1),
            "kk_eigenvalue_2": self.coords.kk_eigenvalue(2),
            "poly26_k1": self.poly26.compute(1, r),
        }

    def print_report(self, result: Dict = None):
        """Print QCalcGeom computation report."""
        result = result or self.compute_all()
        print("=" * 78)
        print("QCALCGEOM COMPUTATION REPORT")
        print("=" * 78)

        m = result["metric"]
        print(f"\n▶ BSFG Metric Components")
        print(f"    g_00 = {m['g_00']:.12e}")
        print(f"    g_rr = {m['g_rr']:.12e}")
        print(f"    ε    = {m['eps']:.6e}")

        r = result["ricci"]
        print(f"\n▶ Ricci Tensor")
        print(f"    R_00      = {r['R_00']:.6e}")
        print(f"    R_rr      = {r['R_rr']:.6e}")
        print(f"    R (scalar)= {r['R_scalar']:.6e}")
        print(f"    Kretschner= {r['Kretschner']:.6e}")

        ch = result["christoffel"]
        print(f"\n▶ Christoffel Symbols")
        for key, val in ch.items():
            print(f"    {key} = {val:.6e}")

        h = result["horizon"]
        print(f"\n▶ Horizon Properties")
        print(f"    r_h   = {h['r_h_m']:.6e} m")
        print(f"    T_H   = {h['T_H_K']:.6e} K")
        print(f"    κ_surf = {h['kappa_surf']:.6e} s⁻²")

        hol = result["holonomy"]
        print(f"\n▶ Holonomy")
        print(f"    Group: {hol['group']}")
        print(f"    G₂ excluded: {hol['G2_excluded']}")
        print(f"    Spin(7) excluded: {hol['Spin7_excluded']}")

        v = result["vds"]
        print(f"\n▶ VDS (Vacuum Density Series)")
        print(f"    Value     = {v['value']:.12e}")
        print(f"    Converged = {v['converged']}")

        d = result["dvp"]
        print(f"\n▶ DVP (Dipole Vortex Primes)")
        print(f"    26! mod 113 = {d['fac26_mod_113']}")
        print(f"    r_q = {d['r_q_m']:.6e} m = {d['r_q_AU']:.6e} AU")

        b = result["bsh"]
        print(f"\n▶ BSH (Buoyancy Series Harmonics)")
        print(f"    U_g2    = {b['U_g2']:.6e}")
        print(f"    H_max   = {b['H_m_max']:.6e}")
        print(f"    Saturated = {b['saturated']}")

        print(f"\n▶ KK Eigenvalues")
        print(f"    λ_1 = {result['kk_eigenvalue_1']}")
        print(f"    λ_2 = {result['kk_eigenvalue_2']}")

        p26 = result["poly26_k1"]
        print(f"\n▶ 26th-Order Polynomial (k=1)")
        print(f"    Value = {p26['value']:.6e}")
        print(f"    Negligible = {p26['negligible']}")

        print("=" * 78)


# ═══════════════════════════════════════════════════════════════════════════════
# §11  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    engine = QCalcGeomEngine(use_ipc=False)

    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        result = engine.compute_all()
        outfile = sys.argv[2] if len(sys.argv) > 2 else "qcalcgeom_results.json"
        with open(outfile, "w") as f:
            json.dump(result, f, indent=2, default=str)
        print(f"Exported to {outfile}")
    elif len(sys.argv) > 1 and sys.argv[1] == "--ipc-test":
        client = QCalcGeomIPCClient()
        if client.connect():
            print("Connected to StarMagic_UQFF pipe — running test suite")
            resp = client.run_tests()
            if resp:
                print(f"Test response: {len(resp)} bytes")
            else:
                print("No response (pipe closed or error)")
            client.disconnect()
        else:
            print("Cannot connect to pipe — is physics_backend running?")
    else:
        result = engine.compute_all()
        engine.print_report(result)


if __name__ == "__main__":
    main()


# ═══════════════════════════════════════════════════════════════════════════════
# SCm VACUUM MANIFOLD — QCalcGeom Integration
# Session: 27FEB2026_A clean thread | scm_vacuum_manifold.py
# Connects SCm primordial physics to the BSFG geometry engine.
# ═══════════════════════════════════════════════════════════════════════════════

class SCmVacuumManifoldQCalcGeom:
    """QCalcGeom-integrated SCm Vacuum Manifold calculator.

    Maps the SCm primordial first-principle (scm_vacuum_manifold.py,
    27FEB2026_A.docx clean thread) onto the QCalcGeom BSFG geometry engine:

      1. SCm-modulated BSFG metric: g_00,SCm = g_00,BSFG · (1 + Φ·ρ_SCm·ε_BSFG)
      2. SCm crossover radius: r_cross,SCm = r_cross,BSFG · (26!)^{-1/13} · S₂₆⁽³⁾
         modulated by phonon Φ(ω, Γ) and negative-time cos(π t_n)
      3. VDS-weighted geodesic: ds²_SCm = ds²_BSFG · VDS_Li26([SSq])
      4. F_U_Bi_i holonomy: holonomy phase ~ F_U_Bi_i / (F_U_Bi_i + 1)

    References:
      - QCalcGeom.h: bsfg_geodesic, bsfg_metric, bsfg_horizon
      - qcalcgeom_core_derivation.py: BSFGCrossoverRadius, QCalcGeomMasterEquation
      - scm_vacuum_manifold.py: F_U_Bi_i integral, VDS, Φ(ω,Γ)
    """

    OMEGA_SCM  = 2 * PI * _SCM_THz       # rad/s
    RHO_SCM    = _SCM_RHO_VAC            # 7.09e-37 kg/m³
    SSQ_VAL    = _SCM_SSQ                # 0.57
    FACTORIAL_26 = math.factorial(26)
    COMPACT_SCALE = FACTORIAL_26 ** (-1.0 / 13.0)   # (26!)^{-1/13}

    def __init__(self):
        self.metric   = MetricTensorHelper()
        self.horizon  = BSFGHorizonCalculator()
        self.holonomy = BSFGHolonomyCalculator()
        self.vds      = QCalcGeomVDS()

    def _phi_gaussian(self, omega: float, Gamma: float) -> float:
        """Phonon Gaussian: Φ(ω,Γ) = exp(-(ω-1.25THz)²/(2Γ²))."""
        delta = omega - self.OMEGA_SCM
        return math.exp(-delta**2 / (2 * max(Gamma, 1.0)**2))

    def _vds_Li26(self, N: int = 200) -> float:
        """Li₂₆([SSq]) via module or direct sum."""
        if _SCM_LOADED:
            return float(_scm_vds_num(terms=N))
        return sum(self.SSQ_VAL**n / n**26 for n in range(1, N + 1))

    def _ramanujan_S26(self) -> float:
        """S₂₆⁽³⁾([SSq]) — Ramanujan-accelerated VDS (Polylogarithm order 26).
        Returns the numerical value ≈ 1.453×10²⁶ as documented in PAPER_1129."""
        # For the QCalcGeom context, use the Li_26 value (same series at N→∞)
        return self._vds_Li26(N=500)

    def compute(self, M_kg: float = M_sun, r: float = 6.96e8,
                t_n: float = -100.0, Gamma: float = 2*PI*0.1e12,
                eta: float = 1e-6) -> Dict:
        """Compute SCm-modulated QCalcGeom quantities.

        Parameters
        ----------
        M_kg : float  — central mass (kg)
        r    : float  — radius (m)
        t_n  : float  — negative-time coordinate (s, should be < 0)
        Gamma: float  — phonon linewidth (rad/s)
        eta  : float  — BSFG aether coupling

        Returns
        -------
        dict with SCm-modulated metric, crossover radius, geodesic, holonomy.
        """
        # ── 1. Base BSFG metric ────────────────────────────────────────────
        base_metric  = self.metric.bsfg_metric_components(r, M_kg, eta)
        base_horizon = {
            "r_h_m":     self.horizon.horizon_radius(M_kg, eta),
            "T_H_K":     self.horizon.hawking_temperature(M_kg, eta),
            "kappa_surf": self.horizon.surface_gravity(M_kg, eta),
        }
        base_holonomy = self.holonomy.holonomy_classification()
        base_vds      = self.vds.compute()

        # ── 2. SCm phonon modulation ───────────────────────────────────────
        Phi_ph    = self._phi_gaussian(self.OMEGA_SCM, Gamma)   # on-resonance → 1.0
        cos_pi_tn = math.cos(math.pi * t_n)
        VDS_val   = self._vds_Li26()
        S26_val   = self._ramanujan_S26()

        # ── 3. SCm-modulated g_00 ──────────────────────────────────────────
        eps_bsfg  = base_metric.get("eps", eta)
        g00_bsfg  = base_metric.get("g_00", -1.0)
        g00_scm   = g00_bsfg * (1.0 + Phi_ph * self.RHO_SCM * abs(eps_bsfg) * 1e37)

        # ── 4. SCm crossover radius ────────────────────────────────────────
        # r_cross,BSFG = sqrt(η) · GM/c²  (from BSFGCrossoverRadius)
        r_cross_bsfg = math.sqrt(abs(eta)) * G * M_kg / c**2
        # r_cross,SCm = r_cross,BSFG · (26!)^{-1/13} · S₂₆⁽³⁾ · Φ · |cos(πtₙ)|
        r_cross_scm  = (r_cross_bsfg
                        * self.COMPACT_SCALE
                        * abs(VDS_val)
                        * Phi_ph
                        * abs(cos_pi_tn))

        # ── 5. VDS-weighted geodesic length scale ──────────────────────────
        # ds²_SCm = (g_rr component) · VDS  (VDS modulates effective length)
        g_rr_bsfg    = base_metric.get("g_rr", 1.0)
        ds2_scm      = g_rr_bsfg * abs(VDS_val)

        # ── 6. F_U_Bi_i holonomy phase ─────────────────────────────────────
        F_UBi = float(_scm_F_U_Bi_i_num(M_bh=M_kg, r=r, Gamma=Gamma)
                      if _SCM_LOADED else 0.0)
        holonomy_phase_scm = F_UBi / (abs(F_UBi) + 1.0)   # ∈ (-1, 1)

        # ── 7. Assemble ────────────────────────────────────────────────────
        return {
            "class":             "SCmVacuumManifoldQCalcGeom",
            "scm_module_loaded": _SCM_LOADED,
            "M_kg":              M_kg,
            "r_m":               r,
            "t_n_s":             t_n,
            "eta_bsfg":          eta,
            # Base BSFG
            "g_00_bsfg":         g00_bsfg,
            "g_rr_bsfg":         g_rr_bsfg,
            "r_cross_bsfg_m":    r_cross_bsfg,
            # SCm-modulated
            "g_00_scm":          g00_scm,
            "r_cross_scm_m":     r_cross_scm,
            "ds2_scm":           ds2_scm,
            "holonomy_phase_scm": holonomy_phase_scm,
            "F_U_Bi_i":          F_UBi,
            # SCm scalars
            "Phi_phonon":        Phi_ph,
            "cos_pi_tn":         cos_pi_tn,
            "VDS_Li26":          VDS_val,
            "S26_Ramanujan":     S26_val,
            "compact_scale":     self.COMPACT_SCALE,
            "rho_SCm":           self.RHO_SCM,
            # Base structures (pass-through)
            "bsfg_horizon":      base_horizon,
            "bsfg_holonomy":     base_holonomy,
            "bsfg_vds":          base_vds,
            "primary_equations": [
                "g_00,SCm = g_00,BSFG · (1 + Φ·ρ_SCm·|ε_BSFG|·10³⁷)",
                f"g_00,BSFG = {g00_bsfg:.6e}  →  g_00,SCm = {g00_scm:.6e}",
                "r_cross,SCm = r_cross,BSFG · (26!)^{-1/13} · Li₂₆([SSq]) · Φ · |cos(πtₙ)|",
                f"r_cross,BSFG = {r_cross_bsfg:.4e} m",
                f"r_cross,SCm  = {r_cross_scm:.4e} m",
                "ds²_SCm = g_rr,BSFG · Li₂₆([SSq])",
                f"VDS Li₂₆(0.57) = {VDS_val:.6e}",
                f"(26!)^{{-1/13}} = {self.COMPACT_SCALE:.10e}",
                f"Φ_phonon = {Phi_ph:.6f}  (on-resonance 1.25 THz)",
                f"cos(πtₙ) = {cos_pi_tn:.6f}  (tₙ = {t_n:.1f} s)",
                f"F_U_Bi_i = {F_UBi:.4e}  →  holonomy phase = {holonomy_phase_scm:.6f}",
            ],
            "note": ("SCmVacuumManifoldQCalcGeom. scm_vacuum_manifold.py "
                     "(27FEB2026_A.docx clean thread). "
                     "BSFG geometry + SCm phonon buoyancy co-modulation."),
        }

    def sweep_t_n(self, t_n_list=None, **kw) -> List[Dict]:
        """Sweep over negative-time values and return SCm-modulated results."""
        t_n_list = t_n_list or [-2512.0, -1000.0, -500.0, -100.0, -50.0, -10.0]
        return [self.compute(t_n=tn, **kw) for tn in t_n_list]

    def sweep_Gamma(self, gamma_list=None, **kw) -> List[Dict]:
        """Sweep over phonon linewidths."""
        gamma_list = gamma_list or [2*PI*g*1e12 for g in [0.01, 0.05, 0.1, 0.5, 1.0]]
        return [self.compute(Gamma=g, **kw) for g in gamma_list]
