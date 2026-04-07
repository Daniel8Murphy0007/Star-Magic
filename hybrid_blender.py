#!/usr/bin/env python3
"""
hybrid_blender.py — VDS/DVP/BSH Hybrid Blending Engine for 7-System Registry
══════════════════════════════════════════════════════════════════════════════

PURPOSE: Execute full hybrid blending runs on the canonical 7-system registry
(SGR1745, SagA*, Tapestry, Westerlund2, Pillars, Rings, StudentGuide) using
the three UQFF number systems:

  VDS  = Vacuum Density Series    (ρ_SCm/ρ_UA progression → [SSq] ≈ 0.57)
  DVP  = Dipole Vortex Primes     (p_special = 113 for hydrogen, Z-indexed)
  BSH  = Buoyancy Saturation Harmonics  (Σ(1/k) f_Ub (1-exp(-SSq·m)) cos(...))

Blending formula (triadic):
  g_hybrid = w_C · g_compressed + w_R · g_resonance + w_B · g_buoyancy
  where weights are VDS-modulated: w_X = f(ρ_SCm/ρ_UA, DVP(Z), BSH(k))

REFERENCES:
  - source4.cpp L1562-1710: Canonical 7-system MUGESystem definitions
  - CondensedPhysics4.py L8651: DipoleVortexPrimeEncodingCalculator
  - CondensedPhysics4.py L3445: YangMillsMassGapVacuumDensityEvolutionCalculator
  - CondensedPhysics3.py L3653: DPMHarmonicBuoyancySeriesCalculator
  - uqff_lagrangian_derivation.py: 9-sector Lagrangian constants

SESSION: 203 | April 7, 2026 | PAPER_859-877 integration
"""

import math
import json
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ═══════════════════════════════════════════════════════════════════════════════
# §1  CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

G       = 6.67430e-11
c       = 2.99792e8
hbar    = 1.05457e-34
PI      = math.pi
M_sun   = 1.98892e30

# UQFF calibrated
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_UA  = 7.09e-36
RHO_SCM = 7.09e-37

# VDS constants
RHO_VAC_SCM = 9.47e-27   # kg/m³ — SCm vacuum density (from CP4 L3445)
RHO_VAC_UA  = 5e-27      # kg/m³ — UA vacuum density

# DVP constants
P_SPECIAL_H = 113         # DVP prime for hydrogen


# ═══════════════════════════════════════════════════════════════════════════════
# §2  7-SYSTEM REGISTRY (mirrors source4.cpp L1562-1710)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class MUGESystem:
    """Python mirror of source4.cpp MUGESystem struct."""
    name: str
    I: float          # moment of inertia proxy
    A: float          # area / cross-section proxy
    omega1: float     # angular frequency 1
    omega2: float     # angular frequency 2
    Vsys: float       # system volume
    vexp: float       # expansion velocity
    t: float          # characteristic time
    z: float          # redshift
    ffluid: float     # fluid fraction
    M: float          # mass (kg)
    r: float          # radius (m)
    B: float          # magnetic field (T)
    Bcrit: float      # critical field (T)
    rho_fluid: float  # fluid density
    g_local: float    # local gravity
    M_DM: float       # dark matter mass
    delta_rho_rho: float  # density perturbation


SEVEN_SYSTEMS = [
    MUGESystem("Magnetar SGR 1745-2900",
               1e21, 3.142e8, 1e-3, -1e-3, 4.189e12, 1e3,
               3.799e10, 0.0009, 1.269e-14, 2.984e30, 1e4,
               1e10, 1e11, 1e-15, 10.0, 0.0, 1e-5),
    MUGESystem("Sagittarius A*",
               1e23, 2.813e30, 1e-5, -1e-5, 3.552e45, 5e6,
               3.786e14, 0.0009, 3.465e-8, 8.155e36, 1e12,
               1e-5, 1e-4, 1e-20, 1e-5, 1e37, 1e-3),
    MUGESystem("Tapestry of Blazing Starbirth",
               1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4,
               3.156e13, 0.0, 1e-12, 1.989e35, 3.086e17,
               1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4),
    MUGESystem("Westerlund 2",
               1e22, 1e35, 1e-4, -1e-4, 1e53, 1e4,
               3.156e13, 0.0, 1e-12, 1.989e35, 3.086e17,
               1e-4, 1e-3, 1e-21, 1e-8, 1e35, 1e-4),
    MUGESystem("Pillars of Creation",
               1e21, 2.813e32, 1e-3, -1e-3, 3.552e48, 2e3,
               3.156e13, 0.0, 8.457e-14, 1.989e32, 9.46e15,
               1e-4, 1e-3, 1e-21, 1e-8, 0.0, 1e-5),
    MUGESystem("Rings of Relativity",
               1e22, 1e35, 1e-4, -1e-4, 1e54, 1e5,
               3.156e14, 0.01, 1e-9, 1.989e36, 3.086e17,
               1e-5, 1e-4, 1e-20, 1e-5, 1e36, 1e-3),
    MUGESystem("Student's Guide to the Universe",
               1e24, 1e52, 1e-6, -1e-6, 1e80, 3e8,
               4.35e17, 0.0, 1e-18, 1e53, 1e26,
               1e-10, 1e-9, 1e-30, 1e-10, 1e53, 1e-6),
]


# ═══════════════════════════════════════════════════════════════════════════════
# §3  VDS — VACUUM DENSITY SERIES
# ═══════════════════════════════════════════════════════════════════════════════

class VacuumDensitySeries:
    """
    VDS engine: computes vacuum density ratio progression.
    Core: VDS_ratio = ρ_SCm / ρ_UA → converges to [SSq] ≈ 0.57
    """

    def __init__(self, rho_scm: float = RHO_VAC_SCM, rho_ua: float = RHO_VAC_UA):
        self.rho_scm = rho_scm
        self.rho_ua = rho_ua

    def ratio(self) -> float:
        """Base VDS ratio ρ_SCm/ρ_UA."""
        return self.rho_scm / self.rho_ua

    def ssq_from_vds(self) -> float:
        """Derive [SSq] from VDS ratio: [SSq] = 1 - exp(-VDS_ratio)."""
        return 1.0 - math.exp(-self.ratio())

    def series_term(self, n: int, t: float = 0.0) -> float:
        """
        n-th VDS series term: VDS_n = (ρ_SCm/ρ_UA)^n / n! × exp(-κt)
        Converges rapidly for n > 5.
        """
        r = self.ratio()
        return (r ** n) / math.factorial(n) * math.exp(-KAPPA * t)

    def truncated_sum(self, N: int = 10, t: float = 0.0) -> float:
        """Sum first N terms of VDS expansion."""
        return sum(self.series_term(n, t) for n in range(N))

    def weight_modulation(self, sys: MUGESystem) -> float:
        """
        VDS weight modulation for a given system.
        Stronger SCm → higher VDS weight → more buoyancy contribution.
        """
        local_ratio = self.ratio() * (1 + sys.B / max(sys.Bcrit, 1e-30))
        return min(local_ratio / (1 + local_ratio), 0.99)


# ═══════════════════════════════════════════════════════════════════════════════
# §4  DVP — DIPOLE VORTEX PRIMES
# ═══════════════════════════════════════════════════════════════════════════════

class DipoleVortexPrimes:
    """
    DVP engine: maps atomic Z to DVP prime encoding.
    Core: p_DVP(Z) = p_special(Z) with magnetic quantum number m_l weighting.
    Proto-H = Proto-Fe at Z_id = 26 (magnetic identity, PAPER_870).
    """

    # First 30 primes for DVP encoding
    PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]

    def __init__(self):
        self.p_special_H = P_SPECIAL_H

    def dvp_prime(self, Z: int) -> int:
        """Return the DVP prime for atomic number Z (1-indexed)."""
        idx = min(Z - 1, len(self.PRIMES) - 1)
        return self.PRIMES[max(idx, 0)]

    def dvp_encoding(self, Z: int) -> float:
        """
        DVP encoding value: p_DVP(Z) × (Z / 26) × SSq
        Normalized to proto-Fe Z_id=26 magnetic identity.
        """
        p = self.dvp_prime(Z)
        return p * (Z / 26.0) * SSQ

    def is_dvp_resonant(self, Z: int) -> bool:
        """Check if Z maps to a DVP prime > 26 (resonant threshold)."""
        return self.dvp_prime(Z) > 26

    def dvp_weight(self, Z: int = 26) -> float:
        """Normalized DVP weight for blending."""
        enc = self.dvp_encoding(Z)
        return enc / (enc + 1.0)

    def proto_fe_identity(self) -> Dict:
        """Proto-H = Proto-Fe at Z_id=26 (PAPER_870)."""
        return {
            "Z_id": 26,
            "element": "Fe (Proto-Hydrogen magnetic identity)",
            "dvp_prime": self.dvp_prime(26),
            "encoding": self.dvp_encoding(26),
            "magnetic_basis": "DPM coherent consciousness (PAPER_876)",
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §5  BSH — BUOYANCY SATURATION HARMONICS
# ═══════════════════════════════════════════════════════════════════════════════

class BuoyancySaturationHarmonics:
    """
    BSH engine: harmonic series with buoyancy saturation.
    Core: BSH(k) = Σ_{j=1}^{k} (1/j) × f_Ub × (1 - exp(-SSq·m)) × cos(2πj/26)

    The 1/j harmonic weights enforce diminishing returns at higher modes.
    Saturation (1 - exp(-SSq·m)) prevents divergence at large mass.
    cos(2πj/26) encodes the 26D layer projection.
    """

    def __init__(self, ssq: float = SSQ, beta_i: float = BETA_I):
        self.ssq = ssq
        self.beta_i = beta_i

    def harmonic_term(self, j: int, m_kg: float, f_Ub: float) -> float:
        """Single BSH harmonic term at mode j."""
        saturation = 1.0 - math.exp(-self.ssq * m_kg / M_sun)
        layer_proj = math.cos(2 * PI * j / 26.0)
        return (1.0 / j) * f_Ub * saturation * layer_proj

    def bsh_sum(self, k_max: int, m_kg: float, f_Ub: float = 1.0) -> float:
        """Sum BSH harmonics from j=1 to k_max."""
        return sum(self.harmonic_term(j, m_kg, f_Ub) for j in range(1, k_max + 1))

    def bsh_weight(self, sys: MUGESystem, k_max: int = 26) -> float:
        """BSH blending weight for a system."""
        f_Ub = self.beta_i * G * sys.M / max(sys.r, 1.0)**2
        bsh = self.bsh_sum(k_max, sys.M, f_Ub)
        return abs(bsh) / (abs(bsh) + 1.0)


# ═══════════════════════════════════════════════════════════════════════════════
# §6  GRAVITY COMPUTATIONS (mirrors source4.cpp)
# ═══════════════════════════════════════════════════════════════════════════════

def compute_compressed_g(sys: MUGESystem) -> float:
    """Compressed MUGE gravity (source4.cpp compute_compressed_MUGE)."""
    if sys.r == 0:
        return 0.0
    base = G * sys.M / (sys.r ** 2)
    H0 = 2.269e-18
    expansion = 1 + H0 * sys.t
    super_adj = 1 - sys.B / max(sys.Bcrit, 1e-30)
    adjusted = base * expansion * super_adj

    Lambda = 1.1e-52
    cosm = Lambda * c * c / 3.0
    quantum = (hbar / 1e-68) * 2.176e-18 * (2 * PI / 4.35e17)
    fluid = sys.rho_fluid * sys.Vsys * sys.g_local
    perturbation = (sys.M + sys.M_DM) * (sys.delta_rho_rho + 3 * G * sys.M / max(sys.r, 1.0)**3)

    return adjusted + cosm + quantum + fluid + perturbation


def compute_resonance_g(sys: MUGESystem) -> float:
    """Resonance MUGE gravity (source4.cpp compute_resonance_MUGE)."""
    FDPM = sys.I * sys.A * (sys.omega1 - sys.omega2)
    fDPM = 1e-3
    Evac_neb = 7.09e-36
    c_res = 1e6
    aDPM = FDPM * fDPM * Evac_neb * c_res * sys.Vsys

    fTHz = 1.25e12
    Evac_ISM = 1e-35
    aTHz = fTHz * Evac_neb * sys.vexp * aDPM / Evac_ISM / c_res

    Delta_Evac = abs(Evac_neb - Evac_ISM)
    avac_diff = Delta_Evac * sys.vexp**2 * aDPM / Evac_neb / c_res**2

    Fsuper = 1e6
    asuper = Fsuper * fTHz * aDPM / Evac_neb / c_res

    UA_SCM = U_UA
    omega_i = 1e-3
    fTRZ = 0.1
    aaether = UA_SCM * omega_i * fTHz * aDPM * (1 + fTRZ)

    k4_res = 1e-20
    Ereact = 1e46 * math.exp(-KAPPA * sys.t)
    freact = 1e-10
    Ug4i = k4_res * Ereact * freact * aDPM / Evac_neb * c_res

    fquantum = 1e-15
    aq = fquantum * Evac_neb * aDPM / Evac_ISM / c_res
    fAether = 1e-12
    aA = fAether * Evac_neb * aDPM / Evac_ISM / c_res
    afluid = sys.ffluid * Evac_neb * sys.Vsys / Evac_ISM / c_res

    H_z = 2.270e-18
    fexp = 2 * PI * H_z * sys.t
    aexp = fexp * Evac_neb * aDPM / Evac_ISM / c_res

    b_worm = 1.0
    a_worm = Evac_neb / (b_worm**2 + sys.r**2)

    return aDPM + aTHz + avac_diff + asuper + aaether + Ug4i + aq + aA + afluid + aexp + fTRZ + a_worm


def compute_buoyancy_g(sys: MUGESystem) -> float:
    """Buoyancy gravity contribution (Ubi from Lagrangian sector 6)."""
    if sys.r == 0:
        return 0.0
    Ug_base = G * sys.M / sys.r**2
    Omega_g = 1.0
    Ubi = -BETA_I * Ug_base * Omega_g * sys.M / max(sys.r, 1.0) * U_UA
    return Ubi


# ═══════════════════════════════════════════════════════════════════════════════
# §7  HYBRID BLENDING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

class HybridBlender:
    """
    Blends compressed, resonance, and buoyancy gravity using VDS/DVP/BSH weights.

    g_hybrid = w_C · g_compressed + w_R · g_resonance + w_B · g_buoyancy

    Weights are VDS-modulated with DVP prime correction and BSH saturation.
    """

    def __init__(self, Z: int = 26):
        self.vds = VacuumDensitySeries()
        self.dvp = DipoleVortexPrimes()
        self.bsh = BuoyancySaturationHarmonics()
        self.Z = Z

    def compute_weights(self, sys: MUGESystem) -> Tuple[float, float, float]:
        """Compute VDS/DVP/BSH-modulated triadic weights."""
        w_vds = self.vds.weight_modulation(sys)
        w_dvp = self.dvp.dvp_weight(self.Z)
        w_bsh = self.bsh.bsh_weight(sys)

        # Triadic normalization
        raw_C = (1 - w_vds) * (1 - w_bsh)
        raw_R = w_dvp * w_vds
        raw_B = w_bsh * w_vds
        total = raw_C + raw_R + raw_B
        if total == 0:
            return (1.0 / 3, 1.0 / 3, 1.0 / 3)
        return (raw_C / total, raw_R / total, raw_B / total)

    def blend(self, sys: MUGESystem) -> Dict:
        """Execute hybrid blending for one system."""
        g_comp = compute_compressed_g(sys)
        g_res = compute_resonance_g(sys)
        g_buoy = compute_buoyancy_g(sys)

        w_C, w_R, w_B = self.compute_weights(sys)
        g_hybrid = w_C * g_comp + w_R * g_res + w_B * g_buoy

        return {
            "system": sys.name,
            "g_compressed": g_comp,
            "g_resonance": g_res,
            "g_buoyancy": g_buoy,
            "w_compressed": w_C,
            "w_resonance": w_R,
            "w_buoyancy": w_B,
            "g_hybrid": g_hybrid,
            "vds_ratio": self.vds.ratio(),
            "dvp_encoding": self.dvp.dvp_encoding(self.Z),
            "bsh_weight": self.bsh.bsh_weight(sys),
        }

    def blend_all(self, systems: List[MUGESystem] = None) -> List[Dict]:
        """Run hybrid blending on all 7 systems."""
        systems = systems or SEVEN_SYSTEMS
        return [self.blend(sys) for sys in systems]

    def print_report(self, results: List[Dict] = None):
        """Print formatted blending report."""
        results = results or self.blend_all()
        print("=" * 90)
        print("VDS/DVP/BSH HYBRID BLENDING — 7-SYSTEM REGISTRY")
        print("=" * 90)
        print(f"  VDS ratio  = {self.vds.ratio():.4f}")
        print(f"  DVP prime  = {self.dvp.dvp_prime(self.Z)} (Z={self.Z})")
        print(f"  [SSq]      = {SSQ}")
        print(f"  β_i        = {BETA_I}")
        print()

        for r in results:
            print(f"▶ {r['system']}")
            print(f"    g_comp    = {r['g_compressed']:.6e}")
            print(f"    g_res     = {r['g_resonance']:.6e}")
            print(f"    g_buoy    = {r['g_buoyancy']:.6e}")
            print(f"    w(C/R/B)  = ({r['w_compressed']:.4f}, {r['w_resonance']:.4f}, {r['w_buoyancy']:.4f})")
            print(f"    g_hybrid  = {r['g_hybrid']:.6e}")
            print()
        print("=" * 90)


# ═══════════════════════════════════════════════════════════════════════════════
# §8  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    blender = HybridBlender(Z=26)

    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        results = blender.blend_all()
        outfile = sys.argv[2] if len(sys.argv) > 2 else "hybrid_blending_results.json"
        with open(outfile, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Exported to {outfile}")
    else:
        blender.print_report()


if __name__ == "__main__":
    main()
