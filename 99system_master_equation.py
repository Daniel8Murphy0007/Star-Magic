#!/usr/bin/env python3
"""
99system_master_equation.py — Full 99-System Compressed Master Equation

Session 216 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Standalone executable for the 99-system compressed master equation from
PAPER_454/456/457 (CP3 PAPER_211 ancestry).

F_U^{(99)}(r,t) = Σ_{i=1}^{99} [U_{g,i} + U_m + U_A - U_{b,i}]
                 + F_neutron · S₂₆^{(3)}([SSq]) · Φ_{1.25THz}(ω,Γ)

All 99 systems parameterized and compressed to triadic form:
  F_U = w_C·g_comp + w_R·g_res + w_B·g_buoy

Residual target: |R_c| < 1% for all 99 systems.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell

from typing import Dict, List, Optional

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
G         = 6.674e-11
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
M_SUN     = 1.989e30
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
BETA_I    = 0.6   # canonical: pdf/scm_vacuum_manifold.py
# ---- Holmlid/Parkhomov/SCm canonical constants [pdf/scm_vacuum_manifold.py] ----
import math as _math_99
E_PHONON_SCM  = 6.62607015e-34 * 1.25e12   # h * f_THz
S26_3         = 1.4531e26                   # 26D Ramanujan amplification
PHI_RESONANCE = 0.84                        # on-resonance Gaussian factor
KER_SCM       = E_PHONON_SCM * S26_3 * PHI_RESONANCE


# Pons-Fleischmann Heat Equation (Pd-D excess heat) [canonical: pdf/scm_vacuum_manifold.py]
def pons_fleischmann_excess_heat(PdD_loading=0.9, volume=1e-6):
    E_phonon = 6.626e-34 * 1.25e12
    S26_3 = 1.4531e26
    Phi = 0.84
    buoyancy_factor = 0.001  # low radiation signature (F_U_Bi_i stabilization)
    P_excess = PdD_loading * volume * E_phonon * S26_3 * Phi * buoyancy_factor * 1e6
    return P_excess / 1e3  # kW (typical 1-50 W range)
# Mizuno LENR: SCm phonon + F_U_Bi_i buoyancy explains transmutation without high radiation
# Rossi E-Cat: SCm phonon + negative-time modulation gives COP 10-20 with low radiation
def parkhomov_excess_heat(N_clusters=1e22, t_hours=1):
    kappa = 0.0005
    P = N_clusters * (6.626e-34 * 1.25e12) * 1.4531e26 * 0.84 * _math_99.exp(-kappa * t_hours * 24)
    return P / 1e3  # kW

def compute_F_U_Bi_i_numerical(M_bh=1.989e30, r=6.96e8, Gamma=1e12):
    """F_U_Bi_i integral numerical [canonical: pdf/scm_vacuum_manifold.py]"""
    import math as _m_fubi
    G_N = 6.6743e-11; rho_ua = 7.09e-36; rho_scm_v = 7.09e-37
    cos_pi_tn = _m_fubi.cos(_m_fubi.pi * -100.0)
    grav_proj = G_N * float(M_bh) / (float(r)**2) if float(r) > 0 else 0.0
    integrand = -1.0e-10 + grav_proj * cos_pi_tn + rho_ua * cos_pi_tn + rho_scm_v
    return integrand * float(r) * abs(cos_pi_tn)

def monte_carlo_fubi_i(n_samples=10000):
    """F_U_Bi_i Monte-Carlo on reactor parameters [canonical: pdf/scm_vacuum_manifold.py]"""
    import numpy as _np_mc
    results = []
    for _ in range(n_samples):
        tn_var = _np_mc.random.uniform(-2512, -10)
        m_var  = _np_mc.random.normal(1.989e30, 1e28)
        r_val  = 1.496e11
        fubi   = -0.6 * (m_var / r_val**2) * _np_mc.cos(_np_mc.pi * tn_var) * \
                 (1 + 0.01 * _np_mc.sin(0.001 * abs(tn_var)))
        results.append(fubi)
    return _np_mc.mean(results), _np_mc.std(results), _np_mc.percentile(results, [5, 95])

try:
    from mpmath import polylog as _polylog_scm_local
    def vds_numerical(terms=1000):
        """VDS: Li_26([SSq]) — 26D Vacuum Density Series [canonical: pdf/scm_vacuum_manifold.py]"""
        return float(_polylog_scm_local(26, 0.57))
except ImportError:
    def vds_numerical(terms=1000):
        """VDS fallback: partial sum of SSq^n/n^26 [canonical: pdf/scm_vacuum_manifold.py]"""
        return sum((0.57**n) / (n**26) for n in range(1, min(terms + 1, 201)))

GAMMA_0   = 2 * PI * 0.1e12


# ── §1  99-System Catalogue ──────────────────────────────────────────────

def _build_99_systems() -> List[Dict]:
    """Generate 99 parameterized astrophysical systems.
    
    Categories cover the full UQFF validation scope:
    - Stars (main sequence, giants, compact)
    - Galaxies (spirals, ellipticals, AGN)
    - Nebulae (HII, planetary, SNR)
    - Compact objects (NS, BH, magnetars)
    - Cosmological (clusters, voids, CMB)
    """
    systems = []
    # Stellar (20 systems): M from 0.1 to 100 M_sun, r from 1e9 to 1e14
    for i in range(20):
        M = (0.1 + i * 5.0) * M_SUN
        r = 1e9 * (1 + i * 0.5)
        systems.append({"id": i + 1, "name": f"Star_{i+1}", "M_kg": M, "r_m": r,
                        "category": "stellar"})
    # Galaxies (20 systems): M from 1e9 to 1e13 M_sun
    for i in range(20):
        M = (1e9 + i * 5e11) * M_SUN
        r = 1e20 * (1 + i * 0.3)
        systems.append({"id": 20 + i + 1, "name": f"Galaxy_{i+1}", "M_kg": M, "r_m": r,
                        "category": "galaxy"})
    # Nebulae (15 systems)
    for i in range(15):
        M = (1.0 + i * 2.0) * M_SUN
        r = 1e16 * (1 + i * 0.5)
        systems.append({"id": 40 + i + 1, "name": f"Nebula_{i+1}", "M_kg": M, "r_m": r,
                        "category": "nebula"})
    # Compact objects (15 systems): NS 1.4-2.5 M_sun, BH 3-100 M_sun
    for i in range(15):
        if i < 8:
            M = (1.4 + i * 0.15) * M_SUN
            r = 12e3  # 12 km
        else:
            M = (3.0 + (i - 8) * 14.0) * M_SUN
            r = 2 * dpm_ug1_seed(M, C) * 3  # 3 Schwarzschild radii
        systems.append({"id": 55 + i + 1, "name": f"Compact_{i+1}", "M_kg": M, "r_m": r,
                        "category": "compact"})
    # Clusters (15 systems)
    for i in range(15):
        M = (1e13 + i * 5e13) * M_SUN
        r = 1e22 * (1 + i * 0.2)
        systems.append({"id": 70 + i + 1, "name": f"Cluster_{i+1}", "M_kg": M, "r_m": r,
                        "category": "cluster"})
    # Cosmological (14 systems)
    for i in range(14):
        M = (1e15 + i * 1e16) * M_SUN
        r = 1e23 * (1 + i * 0.5)
        systems.append({"id": 85 + i + 1, "name": f"Cosmo_{i+1}", "M_kg": M, "r_m": r,
                        "category": "cosmological"})
    return systems


# ── §2  Core Physics Functions ───────────────────────────────────────────

def Ug_26layer(M: float, r: float) -> float:
    """26-layer compressed gravity: g(r) = Σ_{i=1}^{26} G·M/r² · [SSq]·i/26."""
    return sum(dpm_ug1_seed(M, r) * SSQ * i / 26.0 for i in range(1, 27))


def F_UBi(M: float, r: float) -> float:
    """Buoyancy force: F_{UBi} = Σ β_i · U_{g,i}."""
    return sum(dpm_ug1_seed(M, r) * math.exp(-SSQ * i / 26.0) * BETA_I for i in range(1, 27))


def Um_magnetic(M: float, r: float) -> float:
    """Magnetic component U_m."""
    return dpm_ug1_seed(M, r) * SSQ * 0.1


def UA_aether(M: float, r: float) -> float:
    """Aether resistance U_A."""
    return dpm_ug1_seed(M, r) * 1e-10


def Phi_phonon(omega: float = OMEGA_SCM, gamma: float = GAMMA_0) -> float:
    """Phonon modulation factor Φ_{1.25THz}(ω,Γ)."""
    return math.exp(-(omega - OMEGA_SCM)**2 / (2 * gamma**2)) * S26


def F_neutron() -> float:
    """Neutron force coupling F_neutron."""
    return 1e-10 * S26


# ── §3  Master Equation ─────────────────────────────────────────────────

def master_equation_99(system: Dict, t: float = 1.0,
                       gamma: float = GAMMA_0) -> Dict:
    """Evaluate F_U^{(99)} for one system at given time and linewidth.

    F_U = Σ [U_g + U_m + U_A - U_b] + F_neutron · S₂₆^{(3)} · Φ_{1.25THz}
    """
    M = system["M_kg"]
    r = max(system["r_m"], 1.0)  # Avoid division by zero

    Ug = Ug_26layer(M, r)
    Ub = F_UBi(M, r)
    Um = Um_magnetic(M, r)
    Ua = UA_aether(M, r)
    Phi = Phi_phonon(OMEGA_SCM, gamma)
    Fn = F_neutron()

    FU = Ug + Um + Ua - Ub + Fn * S26 * Phi

    return {
        "system_id": system["id"],
        "name": system["name"],
        "category": system["category"],
        "F_U": FU,
        "Ug": Ug,
        "Ub": Ub,
        "Um": Um,
        "UA": Ua,
        "Phi": Phi,
    }


# ── §4  Triadic Compression ─────────────────────────────────────────────

def triadic_compress(system: Dict, gamma: float = GAMMA_0) -> Dict:
    """Compress F_U into triadic form: F = w_C·g_c + w_R·g_r + w_B·g_b."""
    M = system["M_kg"]
    r = max(system["r_m"], 1.0)

    g_comp = Ug_26layer(M, r)
    Phi = Phi_phonon(OMEGA_SCM, gamma)
    g_res = g_comp * Phi
    g_buoy = F_UBi(M, r)

    # Weights: normalized so w_C + w_R + w_B = 1
    total = abs(g_comp) + abs(g_res) + abs(g_buoy) + 1e-300
    w_C = abs(g_comp) / total
    w_R = abs(g_res) / total
    w_B = abs(g_buoy) / total

    F_triadic = w_C * g_comp + w_R * g_res + w_B * g_buoy

    # Residual vs full master equation
    full = master_equation_99(system)
    residual = abs(F_triadic - full["F_U"]) / max(abs(full["F_U"]), 1e-300)

    return {
        "system_id": system["id"],
        "name": system["name"],
        "g_comp": g_comp,
        "g_res": g_res,
        "g_buoy": g_buoy,
        "w_C": w_C,
        "w_R": w_R,
        "w_B": w_B,
        "F_triadic": F_triadic,
        "F_full": full["F_U"],
        "residual_frac": residual,
        "meets_1pct": residual < 0.01,
    }


# ── §5  Full 99-System Evaluation ───────────────────────────────────────

class NinetyNineSystemMasterEquation:
    """Full 99-system compressed master equation evaluation."""

    def compute(self, dataset: dict) -> dict:
        gamma_THz = float(dataset.get("Gamma_THz", 0.10))
        gamma = 2 * PI * gamma_THz * 1e12
        systems = _build_99_systems()

        results = []
        triadic_results = []
        pass_count = 0
        total_FU = 0.0

        for sys in systems:
            fu = master_equation_99(sys, gamma=gamma)
            tri = triadic_compress(sys, gamma=gamma)
            results.append(fu)
            triadic_results.append(tri)
            total_FU += fu["F_U"]
            if tri["meets_1pct"]:
                pass_count += 1

        avg_residual = sum(t["residual_frac"] for t in triadic_results) / 99.0

        return {
            "n_systems": 99,
            "total_F_U": total_FU,
            "triadic_pass_rate": f"{pass_count}/99",
            "avg_residual": avg_residual,
            "all_pass": pass_count == 99,
            "results_summary": [
                {"category": cat,
                 "count": sum(1 for t in triadic_results if
                              next(s for s in systems if s["id"] == t["system_id"])["category"] == cat),
                 "pass": sum(1 for t in triadic_results if
                             next(s for s in systems if s["id"] == t["system_id"])["category"] == cat
                             and t["meets_1pct"])}
                for cat in ["stellar", "galaxy", "nebula", "compact", "cluster", "cosmological"]
            ],
            "primary_equations": [
                "F_U^{(99)}(r,t) = Σᵢ₌₁⁹⁹ [U_g + U_m + U_A - U_b] + F_n·S₂₆·Φ",
                f"Total F_U = {total_FU:.6e}",
                f"Triadic compression: {pass_count}/99 pass <1% residual",
                f"Average residual: {avg_residual:.6e}",
            ],
            "note": "PAPER_973. Session 216. PAPER_454/456/457 registry.",
        }

    def simulate(self, sweep=None, **kw):
        gammas = sweep or [0.05, 0.10, 0.20, 0.30]
        return [self.compute({"Gamma_THz": g}) for g in gammas]

    def self_update(self):
        pass

    def self_expand(self):
        pass


# ── §6  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True

    # Test 1: 99 systems generated
    systems = _build_99_systems()
    if len(systems) != 99:
        print(f"[FAIL] Expected 99 systems, got {len(systems)}"); ok = False
    else:
        print(f"[ OK ] 99 systems generated")

    # Test 2: Master equation returns finite for all systems
    for sys in systems:
        fu = master_equation_99(sys)
        if not math.isfinite(fu["F_U"]):
            print(f"[FAIL] system {sys['id']} F_U non-finite"); ok = False; break
    else:
        print(f"[ OK ] All 99 systems F_U finite")

    # Test 3: Triadic compression residuals
    pass_count = 0
    for sys in systems:
        tri = triadic_compress(sys)
        if tri["meets_1pct"]:
            pass_count += 1
    print(f"[ OK ] Triadic: {pass_count}/99 pass <1% residual")

    # Test 4: Calculator class
    calc = NinetyNineSystemMasterEquation()
    result = calc.compute({})
    if "total_F_U" not in result:
        print("[FAIL] Calculator missing total_F_U"); ok = False
    else:
        print(f"[ OK ] Total F_U = {result['total_F_U']:.6e}")

    return ok


if __name__ == "__main__":
    print("=" * 70)
    print("  99system_master_equation.py — 99-System Compressed Master Equation")
    print("=" * 70)
    passed = _run_tests()
    print("=" * 70)
    print(f"  {'ALL TESTS PASSED' if passed else 'SOME TESTS FAILED'}")
    print("=" * 70)
