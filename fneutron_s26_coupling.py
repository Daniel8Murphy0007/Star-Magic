"""
F_neutron × S₂₆ Buoyancy-Polylog Coupling

Session 204 | Daniel Murphy
PURPOSE: Two forces exist independently in the codebase:

         1. F_neutron^SCm (inspiral_kozima_coupling.py, kozima_scm_cross_section.py):
              F_neutron = N_n · σ_n^SCm(ω,n) · Φ_phonon · (F_{U,Bi}/F_U − 1)
            Kozima neutron-drop force modulated by SCm lattice phonon flux.

         2. S₂₆([SSq]) (ramanujan_polylog_s26.py — new):
              S₂₆(z) = Li₂₆(z) via Euler-Ramanujan acceleration
            26-dimensional vacuum density polylogarithm.

         MISSING: Their product as a buoyancy-polylog force channel.

         Physical motivation: The neutron-drop force occurs WITHIN the
         26-level VDS vacuum structure. Each VDS level n contributes a
         different phonon resonance width and cross-section. The polylog
         S₂₆ encodes the cumulative vacuum density across all 26 levels.
         The coupled force:

           F_coupled(ω) = F_neutron(ω, n) × S₂₆([SSq] × (1 + n/26))

         gives the buoyancy force weighted by the full 26-level vacuum
         density spectrum, integrated over all VDS layers.

         This module:
           1. Computes F_neutron at each VDS level n=0..26
           2. Multiplies by the Ramanujan-accelerated S₂₆ at that level
           3. Integrates the coupled force over all 26 levels
           4. Compares coupled vs uncoupled total forces

ARCHITECTURE: Pure calculator. No hardcoded systems.
"""

import json
import math
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
M_SUN = 1.989e30

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36
RHO_VAC_SCM = 9.47e-27

# VDS
N_LEVELS = 26

# LENR / Kozima parameters (from kozima_scm_cross_section.py)
F_LENR_THZ = 1.25e12
OMEGA_SCM = 2 * PI * F_LENR_THZ
GAMMA_DEFAULT = 0.1e12 * 2 * PI
SIGMA_0 = 1e-4


# ── §1  NEUTRON FORCE AT VDS LEVEL ──────────────────────────────────────

def _f_neutron_at_level(omega: float, n: int,
                        N_n: float, phi_phonon: float,
                        F_U_Bi: float, F_U: float) -> Dict[str, Any]:
    """
    F_neutron^SCm at VDS level n:
      σ_n^SCm(ω,n) = σ₀ · exp[−(ω−ω_SCm)²/(2Γ²)] · (1 + [SSq]·n/26)
      F_neutron = N_n · σ · Φ · (F_{U,Bi}/F_U − 1)
    """
    exponent = -((omega - OMEGA_SCM)**2) / (2 * GAMMA_DEFAULT**2)
    gaussian = math.exp(exponent)
    vds_factor = 1.0 + (SSQ * n) / N_LEVELS
    sigma_scm = SIGMA_0 * gaussian * vds_factor

    buoyancy_rev = (F_U_Bi / F_U - 1.0) if F_U != 0 else 0.0
    F_n = N_n * sigma_scm * phi_phonon * buoyancy_rev

    return {
        "F_neutron_N": F_n,
        "sigma_n_scm": sigma_scm,
        "vds_level": n,
    }


# ── §2  ETA-ACCELERATED S₂₆ (inline, no circular import) ────────────────

def _S26(z: float, N_terms: int = 55) -> float:
    """
    Eta-function accelerated Li_{26}(z) via Euler transform.

    Li_s(z) = η_s(z) + 2^{1-s} Li_s(z²)

    where η_s(z) = -Li_s(-z) is alternating and Euler-accelerable:
      η_s(z) = Σ_{n=0}^N (1/2^{n+1}) × d_n(z,s)
      d_n(z,s) = Σ_{j=0}^n (-1)^{n-j} C(n,j) z^{j+1}/(j+1)^s

    Returns scalar value for coupling computation.
    """
    if abs(z) >= 1.0 or z == 0.0:
        return 0.0

    s = 26

    # Component 1: Euler-accelerated eta function
    eta = 0.0
    inv2 = 0.5  # 1/2^{n+1}
    for n in range(N_terms):
        d_n = 0.0
        for j in range(n + 1):
            sign = (-1)**j
            binom = math.comb(n, j)
            d_n += sign * binom * z**(j + 1) / (j + 1)**s
        term = d_n * inv2
        eta += term
        inv2 *= 0.5
        if n > 2 and abs(term) < 1e-16:
            break

    # Component 2: 2^{1-s} × Li_s(z²)  (negligible for s=26)
    z_sq = z * z
    correction = 0.0
    if abs(z_sq) < 1.0:
        factor = 2.0**(1 - s)
        for k in range(1, 6):
            correction += factor * z_sq**k / k**s

    return eta + correction


# ── §3  COUPLED FORCE CALCULATOR ────────────────────────────────────────

class FneutronS26Coupling:
    """
    Buoyancy-polylog coupled force: F_neutron(n) × S₂₆(z_n).

    At each VDS level n:
      z_n = [SSq] × (1 + n/26)          (clamped to |z| < 1)
      F_coupled(n) = F_neutron(ω, n) × S₂₆(z_n) × (2π)^{n/6}

    The (2π)^{n/6} volume scaling comes from PAPER_855
    (PseudoMonopole26StateVacuumDensityCalc).

    Total coupled force:
      F_total = Σ_{n=0}^{26} F_coupled(n)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:      angular frequency (rad/s, default ω_SCm)
          N_n:        neutron number density (1/m³)
          phi_phonon: phonon flux (1/m²·s)
          F_U_Bi:     buoyancy force (N)
          F_U:        unified field force (N)
          SSq:        string squeezing (default 0.57)
          N_levels:   VDS dimensionality (default 26)
        """
        omega = dataset.get('omega', OMEGA_SCM)
        N_n = dataset.get('N_n', 1e28)
        phi_phonon = dataset.get('phi_phonon', 1e16)
        F_U_Bi = dataset.get('F_U_Bi', 1.1)
        F_U = dataset.get('F_U', 1.0)
        ssq = dataset.get('SSq', SSQ)
        n_levels = dataset.get('N_levels', N_LEVELS)

        levels = []
        F_total_coupled = 0.0
        F_total_uncoupled = 0.0

        for n in range(n_levels + 1):
            # F_neutron at this level
            fn = _f_neutron_at_level(omega, n, N_n, phi_phonon, F_U_Bi, F_U)
            F_n = fn["F_neutron_N"]

            # VDS-modulated argument
            z_n = ssq * (1.0 + n / n_levels)
            if z_n >= 1.0:
                z_n = 0.999

            # Ramanujan-accelerated S₂₆
            s26_val = _S26(z_n)

            # Volume scaling (PAPER_855)
            delta_n = (2 * PI)**(n / 6.0)

            # Coupled force at this level
            F_coupled_n = F_n * s26_val * delta_n

            F_total_coupled += F_coupled_n
            F_total_uncoupled += F_n

            levels.append({
                "n": n,
                "F_neutron_n": F_n,
                "z_n": z_n,
                "S_26_n": s26_val,
                "delta_n": delta_n,
                "F_coupled_n": F_coupled_n,
                "sigma_n": fn["sigma_n_scm"],
            })

        # Coupling amplification ratio
        amp_ratio = (F_total_coupled / F_total_uncoupled
                     if F_total_uncoupled != 0 else 0.0)

        return {
            "F_total_coupled": F_total_coupled,
            "F_total_uncoupled": F_total_uncoupled,
            "amplification_ratio": amp_ratio,
            "N_levels": n_levels,
            "omega": omega,
            "SSq": ssq,
            "levels": levels,
            "equation": (
                "F_coupled(n) = F_neutron(ω,n) × S₂₆([SSq]·(1+n/26)) "
                "× (2π)^{n/6}"
            ),
            "total_equation": "F_total = Σ_{n=0}^{26} F_coupled(n)",
        }

    def compute_damping_extension(self, dataset: dict) -> Dict[str, Any]:
        """
        Extend 5-channel GW damping (inspiral_kozima_coupling.py) using
        the S₂₆-coupled force instead of the uncoupled F_neutron.

          D_Kozima_S26 = 1 / (1 + |F_total_coupled| / F_GW)

        This replaces the simple D_Kozima from inspiral_kozima_coupling.py
        with a force that includes the polylogarithmic VDS weighting.
        """
        coupled = self.compute(dataset)
        F_GW = dataset.get('F_GW_N', 1e40)

        D_AETHER = 1.000
        D_SCM_CH = 1.000
        D_TRZ = 0.900
        D_STRING = 0.370
        D_4ch = D_AETHER * D_SCM_CH * D_TRZ * D_STRING

        ratio = abs(coupled["F_total_coupled"]) / F_GW if F_GW != 0 else 0
        D_Kozima_S26 = 1.0 / (1.0 + ratio)
        D_5ch_S26 = D_4ch * D_Kozima_S26

        return {
            "D_Kozima_S26": D_Kozima_S26,
            "D_total_5ch_S26": D_5ch_S26,
            "D_total_4ch": D_4ch,
            "F_total_coupled": coupled["F_total_coupled"],
            "F_GW_N": F_GW,
            "ratio": ratio,
            "strain_reduction_5ch_S26_pct": (1.0 - D_5ch_S26) * 100,
            "equation": "D_Kozima_S26 = 1/(1 + |F_coupled_total|/F_GW)",
        }


# ── §4  SELF-TEST ────────────────────────────────────────────────────────

def self_test() -> Dict[str, Any]:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    calc = FneutronS26Coupling()

    # Test 1: default dataset (at resonance)
    result = calc.compute({
        'omega': OMEGA_SCM,
        'N_n': 1e28,
        'phi_phonon': 1e16,
        'F_U_Bi': 1.1,
        'F_U': 1.0,
    })

    # Test 2: damping extension
    damping = calc.compute_damping_extension({
        'omega': OMEGA_SCM,
        'N_n': 1e28,
        'phi_phonon': 1e16,
        'F_U_Bi': 1.1,
        'F_U': 1.0,
        'F_GW_N': 1e40,
    })

    # Test 3: S₂₆ consistency check
    s26_ssq = _S26(SSQ)
    # Compare with naive Li₂₆
    naive_li26 = sum(SSQ**k / k**26 for k in range(1, 100))

    checks = {
        "levels_computed": len(result["levels"]),
        "F_total_coupled_nonzero": result["F_total_coupled"] != 0,
        "F_total_uncoupled_nonzero": result["F_total_uncoupled"] != 0,
        "amplification_ratio": result["amplification_ratio"],
        "D_Kozima_S26_valid": 0 < damping["D_Kozima_S26"] <= 1,
        "S26_SSq": s26_ssq,
        "naive_Li26_SSq": naive_li26,
        "S26_naive_agree": abs(s26_ssq - naive_li26) / abs(naive_li26) < 1e-10
                           if naive_li26 != 0 else True,
    }

    all_pass = (
        checks["levels_computed"] == 27
        and checks["F_total_coupled_nonzero"]
        and checks["D_Kozima_S26_valid"]
        and checks["S26_naive_agree"]
    )

    return {
        "module": "fneutron_s26_coupling",
        "timestamp": ts,
        "status": "PASS" if all_pass else "FAIL",
        "checks": checks,
        "coupled_force": {
            "F_total_coupled": result["F_total_coupled"],
            "F_total_uncoupled": result["F_total_uncoupled"],
            "amplification_ratio": result["amplification_ratio"],
        },
        "damping": {
            "D_Kozima_S26": damping["D_Kozima_S26"],
            "D_total_5ch_S26": damping["D_total_5ch_S26"],
        },
    }


# ── §5  CLI ENTRY ────────────────────────────────────────────────────────

if __name__ == "__main__":
    result = self_test()
    print(json.dumps(result, indent=2, default=str))
