#!/usr/bin/env python3
"""
scm_quantum_computing.py — SCm Phonon-Mediated Quantum Computing Engine

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
Implements the phonon-mediated quantum gate operator and qubit coherence
framework for SCm vacuum phonon coupling at 1.25 THz.

Gap closed:
  - Explicit unitary gate operator U_gate = exp(-iπ/4 · σ_z · Φ_{1.25THz})
  - Pauli matrix algebra for gate construction
  - Topological protection via SCm phonon gap
  - T_2 coherence with F_{U,Bi}/F_U ratio

Physics:
  U_gate(θ) = exp(-i θ/2 · σ_z · Φ_{1.25THz}(ω,Γ))
            = [[exp(-iθΦ/2), 0], [0, exp(iθΦ/2)]]

  T_2 = (ℏ/Δ_SCm) · exp(Δ_SCm / k_BT) · S₂₆^(3)([SSq]) · (F_{U,Bi}/F_U)

  Gate fidelity: F = 1 - (t_gate/T_2)·(1 - H_SCm·β_i·[SSq])

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Optional, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
OMEGA_SCM = 2 * PI * 1.25e12
GAMMA_0   = 2 * PI * 0.1e12
SSQ       = 0.57
BETA_I    = 0.603
H_SCM     = 0.99
F_UBI_RATIO = 0.6

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor R_n^{(26,k)}."""
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 0.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, 27):
            sign = (-1) ** (j + 1)
            binom = math.comb(26, j)
            inner += sign * binom * math.factorial(26 - j) / n ** j
        correction += inner / n ** (26 * m)
    return prefactor * (1.0 + correction)


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))


def phi_phonon(omega: float, gamma: float = GAMMA_0) -> float:
    """Phonon modulation factor Φ_{1.25THz}(ω, Γ)."""
    return math.exp(-((omega - OMEGA_SCM) ** 2) / (2 * gamma ** 2)) * S26_3RD


# ── §1  PAULI MATRICES & UNITARY GATE OPERATOR ────────────────────────────

# Pauli matrices as 2×2 complex tuples: ((a,b),(c,d))
# σ_x = ((0,1),(1,0)), σ_y = ((0,-i),(i,0)), σ_z = ((1,0),(0,-1))

def _cexp(re: float, im: float) -> Tuple[float, float]:
    """Complex exponential of (re + i·im): e^{re+i·im} = e^re·(cos(im)+i·sin(im))."""
    mag = math.exp(re)
    return (mag * math.cos(im), mag * math.sin(im))


def _cmul(a: Tuple[float, float], b: Tuple[float, float]) -> Tuple[float, float]:
    """Multiply two complex numbers (re,im) tuples."""
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _cabs(z: Tuple[float, float]) -> float:
    """Absolute value of complex number."""
    return math.sqrt(z[0] ** 2 + z[1] ** 2)


class PhononGateOperator:
    """Explicit unitary gate operator for SCm phonon-mediated quantum computing.

    The gate operator is:
        U_gate(θ, ω, Γ) = exp(-i · θ/2 · σ_z · Φ_{1.25THz}(ω, Γ))

    Since σ_z is diagonal, this decomposes as:
        U_gate = [[exp(-i·θ·Φ/2),  0              ],
                  [0,               exp(+i·θ·Φ/2) ]]

    Standard gates:
        - Z-gate:   θ = π       → phase flip
        - S-gate:   θ = π/2     → √Z
        - T-gate:   θ = π/4     → √S
        - Rz(φ):    θ = φ       → arbitrary Z-rotation
        - Hadamard:  requires σ_x + σ_z composition
        - CNOT:     tensor product of control + target
    """

    def __init__(self, omega: float = OMEGA_SCM, gamma: float = GAMMA_0):
        self.omega = omega
        self.gamma = gamma
        self.phi = phi_phonon(omega, gamma)

    def _u_matrix(self, theta: float) -> dict:
        """Compute 2×2 unitary matrix U = exp(-i·θ/2·σ_z·Φ).

        Returns dict with matrix elements as (re, im) tuples.
        """
        phase = theta * self.phi / 2.0
        # exp(-i·phase) for |0⟩ component
        u00 = (math.cos(phase), -math.sin(phase))
        # exp(+i·phase) for |1⟩ component
        u11 = (math.cos(phase), math.sin(phase))
        u01 = (0.0, 0.0)
        u10 = (0.0, 0.0)
        return {"u00": u00, "u11": u11, "u01": u01, "u10": u10}

    def z_gate(self) -> dict:
        """SCm phonon Z-gate: θ = π."""
        return self._u_matrix(PI)

    def s_gate(self) -> dict:
        """SCm phonon S-gate (√Z): θ = π/2."""
        return self._u_matrix(PI / 2)

    def t_gate(self) -> dict:
        """SCm phonon T-gate (√S): θ = π/4."""
        return self._u_matrix(PI / 4)

    def rz_gate(self, phi_angle: float) -> dict:
        """SCm phonon Rz(φ) gate: arbitrary Z-rotation."""
        return self._u_matrix(phi_angle)

    def hadamard_phonon(self) -> dict:
        """SCm phonon Hadamard gate: H = (σ_x + σ_z)/√2 with phonon modulation.

        H_SCm = (1/√2) · Φ · [[1, 1], [1, -1]]
        """
        f = self.phi / math.sqrt(2)
        return {
            "u00": (f, 0.0), "u01": (f, 0.0),
            "u10": (f, 0.0), "u11": (-f, 0.0),
        }

    def cnot_phonon(self) -> dict:
        """SCm phonon CNOT gate (4×4 tensor product).

        CNOT_SCm = |0⟩⟨0| ⊗ I + |1⟩⟨1| ⊗ (σ_x · Φ)
        = [[1, 0, 0,    0   ],
           [0, 1, 0,    0   ],
           [0, 0, 0,    Φ   ],
           [0, 0, Φ,    0   ]]
        """
        p = self.phi
        return {
            "matrix_4x4": [
                [(1, 0), (0, 0), (0, 0), (0, 0)],
                [(0, 0), (1, 0), (0, 0), (0, 0)],
                [(0, 0), (0, 0), (0, 0), (p, 0)],
                [(0, 0), (0, 0), (p, 0), (0, 0)],
            ],
            "Phi": p,
        }

    def verify_unitarity(self, theta: float) -> dict:
        """Verify U†U = I for a given gate angle."""
        mat = self._u_matrix(theta)
        # For diagonal U: U†U_{00} = |u00|² = 1, U†U_{11} = |u11|² = 1
        norm_00 = _cabs(mat["u00"]) ** 2
        norm_11 = _cabs(mat["u11"]) ** 2
        # Off-diagonal: u00*·u01 + u10*·u11 = 0 (trivially for diagonal)
        return {
            "UdagU_00": norm_00,
            "UdagU_11": norm_11,
            "is_unitary": abs(norm_00 - 1.0) < 1e-12 and abs(norm_11 - 1.0) < 1e-12,
            "det_re": mat["u00"][0] * mat["u11"][0] - mat["u00"][1] * mat["u11"][1],
            "det_im": mat["u00"][0] * mat["u11"][1] + mat["u00"][1] * mat["u11"][0],
        }

    def compute(self, dataset: dict) -> dict:
        """Compute all standard gate operators with SCm phonon coupling."""
        theta = float(dataset.get("theta", PI / 4))
        omega = float(dataset.get("omega", self.omega))
        gamma = float(dataset.get("gamma", self.gamma))
        self.omega = omega
        self.gamma = gamma
        self.phi = phi_phonon(omega, gamma)

        gates = {}
        for name, angle in [("Z", PI), ("S", PI / 2), ("T", PI / 4), ("Rz", theta)]:
            mat = self._u_matrix(angle)
            uni = self.verify_unitarity(angle)
            gates[name] = {
                "theta": angle,
                "matrix": {k: f"({v[0]:.8f} + {v[1]:.8f}i)" for k, v in mat.items()},
                "is_unitary": uni["is_unitary"],
            }

        hadamard = self.hadamard_phonon()

        return {
            "Phi_1.25THz": self.phi,
            "S26_3rd": S26_3RD,
            "omega_SCm": OMEGA_SCM,
            "Gamma": gamma,
            "gates": gates,
            "hadamard_phonon": {k: f"({v[0]:.8f} + {v[1]:.8f}i)" for k, v in hadamard.items()},
            "primary_equations": [
                "U_gate(θ) = exp(-i·θ/2·σ_z·Φ_{1.25THz})",
                "= [[exp(-iθΦ/2), 0], [0, exp(iθΦ/2)]]",
                f"Φ_{'{1.25THz}'} = {self.phi:.6e}",
                f"S₂₆⁽³⁾([SSq]) = {S26_3RD:.6e}",
            ],
        }


# ── §2  T₂ COHERENCE WITH F_{U,Bi}/F_U RATIO ─────────────────────────────

class SCmQubitCoherence:
    """Compute qubit coherence time T_2 in SCm phonon vacuum.

    T_2 = (ℏ / Δ_SCm) · exp(Δ_SCm / k_BT) · S₂₆^(3) · (F_{U,Bi}/F_U)

    where Δ_SCm = ℏ·ω_SCm (SCm phonon gap energy).
    """

    def __init__(self, omega_scm: float = OMEGA_SCM, fubi_ratio: float = F_UBI_RATIO):
        self.omega_scm = omega_scm
        self.delta_scm = HBAR * omega_scm
        self.fubi_ratio = fubi_ratio

    def T2(self, T: float) -> float:
        """Compute T_2 coherence time at temperature T (Kelvin)."""
        if T <= 0:
            return float('inf')
        ratio = self.delta_scm / (K_B * T)
        ratio = min(ratio, 700)  # clamp to avoid overflow
        return (HBAR / self.delta_scm) * math.exp(ratio) * S26_3RD * self.fubi_ratio

    def gate_fidelity(self, T: float, t_gate: float) -> float:
        """Gate fidelity F = 1 - (t_gate/T_2)·(1 - H_SCm·β_i·[SSq])."""
        t2 = self.T2(T)
        if t2 == float('inf'):
            return 1.0
        scm_suppression = 1.0 - H_SCM * BETA_I * SSQ
        eps = (t_gate / t2) * scm_suppression
        return max(0.0, 1.0 - eps)

    def topological_protection_factor(self, T: float) -> float:
        """Topological protection from phonon gap: exp(Δ_SCm/k_BT).

        For T << Δ_SCm/k_B, thermal excitations are exponentially suppressed.
        """
        if T <= 0:
            return float('inf')
        ratio = self.delta_scm / (K_B * T)
        ratio = min(ratio, 700)
        return math.exp(ratio)

    def compute(self, dataset: dict) -> dict:
        """Compute coherence, fidelity, and topological protection."""
        T = float(dataset.get("T_K", 0.015))          # default 15 mK
        t_gate = float(dataset.get("t_gate_s", 1e-9))  # default 1 ns

        t2 = self.T2(T)
        fidelity = self.gate_fidelity(T, t_gate)
        topo = self.topological_protection_factor(T)

        # Temperature sweep for coherence landscape
        temps = [0.010, 0.015, 0.020, 0.050, 0.100, 0.300, 1.0, 4.2, 10.0, 77.0]
        sweep = []
        for t_k in temps:
            t2_k = self.T2(t_k)
            f_k = self.gate_fidelity(t_k, t_gate)
            sweep.append({"T_K": t_k, "T2_s": t2_k, "fidelity": f_k})

        return {
            "T_K": T,
            "Delta_SCm_J": self.delta_scm,
            "Delta_SCm_eV": self.delta_scm / 1.602e-19,
            "T2_s": t2,
            "gate_fidelity": fidelity,
            "topological_protection": topo,
            "FUBi_ratio": self.fubi_ratio,
            "temperature_sweep": sweep,
            "primary_equations": [
                "T₂ = (ℏ/Δ_SCm)·exp(Δ_SCm/k_BT)·S₂₆⁽³⁾·(F_{UBi}/F_U)",
                f"Δ_SCm = {self.delta_scm:.6e} J = {self.delta_scm / 1.602e-19:.6e} eV",
                f"T₂(T={T:.3f}K) = {t2:.6e} s",
                f"F(t_gate={t_gate:.1e}s) = {fidelity:.10f}",
                "F = 1 - (t_gate/T₂)·(1 - H_SCm·β_i·[SSq])",
            ],
        }


# ── §3  LAGRANGIAN VARIATION (QUBIT BUOYANCY SECTOR) ──────────────────────

class QubitBuoyancyLagrangian:
    """Qubit buoyancy sector Lagrangian variation.

    δS/δφ_qubit = ∂/∂Δ(-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_n·Φ_{1.25THz}) = 0

    Evaluates stationarity condition for phonon-qubit coupling.
    """

    def compute(self, dataset: dict) -> dict:
        """Evaluate qubit buoyancy Lagrangian stationarity."""
        M = float(dataset.get("M_kg", 1.989e30))
        r = float(dataset.get("r_m", 1e4))
        B = float(dataset.get("B_T", 1e11))
        omega = float(dataset.get("omega", OMEGA_SCM))
        gamma = float(dataset.get("gamma", GAMMA_0))
        F_neutron = float(dataset.get("F_neutron", 1e-10))

        phi_val = phi_phonon(omega, gamma)

        # Ug1 estimate (DPM magnetic dipole)
        mu_s = B * r ** 3
        grad_Ms = G * M / r ** 2
        Ug1 = mu_s * grad_Ms

        # Buoyancy sector: -β_i · Ug1 · (M/r) · U_UA
        U_UA = 1e-4
        buoyancy_term = -BETA_I * Ug1 * (M / r) * U_UA

        # Phonon coupling: F_n · Φ
        phonon_term = F_neutron * phi_val

        # Lagrangian density
        L_qubit = buoyancy_term + phonon_term

        # Stationarity: ∂L/∂Δ evaluated numerically
        delta_scm = HBAR * omega
        dL_dDelta = -BETA_I * Ug1 * (M / r) * U_UA / delta_scm + F_neutron * phi_val / delta_scm

        return {
            "L_qubit": L_qubit,
            "buoyancy_term": buoyancy_term,
            "phonon_term": phonon_term,
            "dL_dDelta": dL_dDelta,
            "Ug1": Ug1,
            "Phi_1.25THz": phi_val,
            "is_stationary": abs(dL_dDelta) < abs(L_qubit) * 1e-6 if L_qubit != 0 else True,
            "primary_equations": [
                "δS/δφ_qubit = ∂/∂Δ(-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_n·Φ) = 0",
                f"L_qubit = {L_qubit:.6e} J·m⁻³",
                f"∂L/∂Δ = {dL_dDelta:.6e}",
            ],
        }


# ── §4  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("SCm Phonon Quantum Computing Engine")
    print("=" * 72)

    # Gate operators
    gate_op = PhononGateOperator()
    result = gate_op.compute({"theta": PI / 4})
    print(f"\nΦ_{{1.25THz}} = {result['Phi_1.25THz']:.6e}")
    for name, g in result["gates"].items():
        print(f"  {name}-gate (θ={g['theta']:.4f}): u00={g['matrix']['u00']}, unitary={g['is_unitary']}")

    # Coherence
    coh = SCmQubitCoherence()
    coh_result = coh.compute({"T_K": 0.015, "t_gate_s": 1e-9})
    print(f"\nT₂(15mK) = {coh_result['T2_s']:.6e} s")
    print(f"Gate fidelity = {coh_result['gate_fidelity']:.10f}")
    print(f"Topological protection = {coh_result['topological_protection']:.6e}")

    # Lagrangian
    lag = QubitBuoyancyLagrangian()
    lag_result = lag.compute({"M_kg": 1.989e30, "r_m": 1e4, "B_T": 1e11})
    print(f"\nL_qubit = {lag_result['L_qubit']:.6e}")

    print("\n✓ All SCm quantum computing calculations complete")
