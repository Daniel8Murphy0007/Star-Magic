#!/usr/bin/env python3
"""
bsm_bounds_derivation.py — BSM Mass Bound Derivation from SCm Phonon Gap

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
First-principles derivation of BSM particle mass bounds from the SCm phonon
gap energy, with full error propagation.

Gap closed:
  - Standalone derivation of m_BSM = (ℏω_SCm/c²)·S₂₆^(3)·(2F_{U,Bi}/F_U - 1)
  - Error propagation from ω_SCm, [SSq], F_{U,Bi}/F_U uncertainties
  - Comparison against experimental bounds (tau/CKM/LFV/VLQ from arXiv 2506.*)

Physics:
  The SCm phonon gap Δ_SCm = ℏω_SCm sets a natural mass scale for BSM
  particles via the UQFF vacuum structure. The buoyancy ratio (2F_{U,Bi}/F_U - 1)
  encodes the asymmetry between emergent and fundamental force sectors.

  m_BSM = (ℏω_SCm / c²) · S₂₆^(3)([SSq]) · (2·F_{U,Bi}/F_U - 1)

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
eV        = 1.602e-19       # J per eV
GeV       = 1.602e-10       # J per GeV

OMEGA_SCM = 2 * PI * 1.25e12
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

# Experimental BSM bounds (from bsm_physics_validation.py / arXiv June 2025)
EXPERIMENTAL_BOUNDS = {
    "tau_g-2": {
        "observable": "Re(a_τ)",
        "range": (-4.5e-3, 6.9e-3),
        "SM_prediction": 1.17721e-3,
        "source": "arXiv:2506.15245",
    },
    "CKM_Vcb": {
        "observable": "|V_cb|",
        "value": 39.2e-3,
        "error": 0.9e-3,
        "source": "arXiv:2506.15256",
    },
    "LFV_B_Kstar": {
        "observable": "BR(B→K*τe)",
        "upper_limit": 5.9e-6,
        "CL": 0.90,
        "source": "arXiv:2506.15347",
    },
    "VLQ_mass": {
        "observable": "m_VLQ",
        "range_GeV": (1150, 2600),
        "kappa_range": (0.22, 0.52),
        "source": "arXiv:2506.15515",
    },
    "neutrino_polarizability": {
        "observable": "P_ν",
        "upper_limit_cm3": 1e-32,
        "source": "arXiv:2506.14881",
    },
}


# ── §1  BSM MASS BOUND DERIVATION ─────────────────────────────────────────

class BSMMassBound:
    """Derive BSM particle mass bounds from SCm phonon gap.

    Derivation:
    ───────────
    1. The SCm vacuum has a phonon gap at ω_SCm = 2π × 1.25 THz.
    2. This sets the natural energy scale: Δ_SCm = ℏω_SCm ≈ 8.29e-22 J ≈ 5.18e-3 eV.
    3. The UQFF mass scale is modulated by:
       - S₂₆^(3)([SSq]): Ramanujan-accelerated vacuum spectral weight
       - (2·F_{U,Bi}/F_U - 1): buoyancy asymmetry factor

    4. BSM mass bound:
       m_BSM = (ℏω_SCm / c²) · S₂₆^(3) · (2·F_{U,Bi}/F_U - 1)

    5. Error propagation via quadrature:
       (δm/m)² = (δω/ω)² + (δS₂₆/S₂₆)² + (2·δr/(2r-1))²
       where r = F_{U,Bi}/F_U
    """

    def __init__(self, omega_scm: float = OMEGA_SCM, fubi_ratio: float = F_UBI_RATIO):
        self.omega_scm = omega_scm
        self.fubi_ratio = fubi_ratio
        self.delta_scm = HBAR * omega_scm

    def mass_bound(self, fubi_ratio: float = None) -> float:
        """Compute m_BSM in kg."""
        r = fubi_ratio if fubi_ratio is not None else self.fubi_ratio
        return (self.delta_scm / C ** 2) * S26_3RD * (2 * r - 1)

    def mass_bound_eV(self, fubi_ratio: float = None) -> float:
        """Compute m_BSM in eV/c²."""
        return self.mass_bound(fubi_ratio) * C ** 2 / eV

    def error_propagation(self, delta_omega_frac: float = 0.01,
                          delta_ssq: float = 0.03,
                          delta_fubi: float = 0.05) -> dict:
        """Propagate uncertainties to m_BSM.

        Parameters:
            delta_omega_frac: fractional uncertainty in ω_SCm
            delta_ssq: absolute uncertainty in [SSq]
            delta_fubi: absolute uncertainty in F_{U,Bi}/F_U ratio
        """
        m0 = self.mass_bound()
        r = self.fubi_ratio

        # Partial derivatives (numerical)
        # ∂m/∂ω via chain rule
        dm_domega = m0 / self.omega_scm  # linear dependence
        sigma_omega = delta_omega_frac * self.omega_scm

        # ∂m/∂[SSq] via numerical diff
        ssq_hi = SSQ + delta_ssq
        ssq_lo = SSQ - delta_ssq
        s26_hi = sum((ssq_hi ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))
        s26_lo = sum((ssq_lo ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))
        dm_dssq = (self.delta_scm / C ** 2) * (s26_hi - s26_lo) / (2 * delta_ssq) * (2 * r - 1)

        # ∂m/∂r
        dm_dr = (self.delta_scm / C ** 2) * S26_3RD * 2

        # Total uncertainty (quadrature)
        sigma_m = math.sqrt(
            (dm_domega * sigma_omega) ** 2 +
            (dm_dssq * delta_ssq) ** 2 +
            (dm_dr * delta_fubi) ** 2
        )

        return {
            "m_BSM_kg": m0,
            "m_BSM_eV": m0 * C ** 2 / eV,
            "sigma_m_kg": sigma_m,
            "sigma_m_eV": sigma_m * C ** 2 / eV,
            "relative_error": sigma_m / abs(m0) if m0 != 0 else float('inf'),
            "contributions": {
                "omega_frac": (dm_domega * sigma_omega) ** 2 / sigma_m ** 2 if sigma_m > 0 else 0,
                "SSq_frac": (dm_dssq * delta_ssq) ** 2 / sigma_m ** 2 if sigma_m > 0 else 0,
                "FUBi_frac": (dm_dr * delta_fubi) ** 2 / sigma_m ** 2 if sigma_m > 0 else 0,
            },
        }

    def compare_experimental(self) -> dict:
        """Compare m_BSM against experimental BSM bounds."""
        m_bsm_eV = self.mass_bound_eV()
        m_bsm_GeV = m_bsm_eV / 1e9

        comparisons = {}
        for name, data in EXPERIMENTAL_BOUNDS.items():
            comparison = {"source": data["source"]}
            if "range_GeV" in data:
                lo, hi = data["range_GeV"]
                comparison["experimental_range_GeV"] = (lo, hi)
                comparison["m_BSM_GeV"] = m_bsm_GeV
                comparison["within_range"] = lo <= m_bsm_GeV <= hi
                comparison["ratio_to_lower"] = m_bsm_GeV / lo if lo > 0 else float('inf')
            elif "upper_limit" in data:
                comparison["upper_limit"] = data["upper_limit"]
                comparison["observable"] = data["observable"]
            elif "value" in data:
                comparison["measured"] = data["value"]
                comparison["error"] = data.get("error", 0)
                comparison["observable"] = data["observable"]
            comparisons[name] = comparison

        return comparisons

    def fubi_ratio_sweep(self, r_min: float = 0.1, r_max: float = 1.0, n_points: int = 20) -> list:
        """Sweep F_{U,Bi}/F_U ratio and compute mass at each point."""
        results = []
        for i in range(n_points + 1):
            r = r_min + (r_max - r_min) * i / n_points
            m = self.mass_bound(r)
            results.append({
                "FUBi_ratio": r,
                "m_BSM_kg": m,
                "m_BSM_eV": m * C ** 2 / eV,
                "asymmetry_factor": 2 * r - 1,
            })
        return results

    def compute(self, dataset: dict) -> dict:
        """Full BSM mass bound derivation with error propagation."""
        r = float(dataset.get("FUBi_ratio", self.fubi_ratio))
        self.fubi_ratio = r

        m = self.mass_bound()
        m_eV = self.mass_bound_eV()
        errors = self.error_propagation()
        comparisons = self.compare_experimental()
        sweep = self.fubi_ratio_sweep()

        return {
            "m_BSM_kg": m,
            "m_BSM_eV": m_eV,
            "m_BSM_GeV": m_eV / 1e9,
            "Delta_SCm_J": self.delta_scm,
            "Delta_SCm_eV": self.delta_scm / eV,
            "S26_3rd": S26_3RD,
            "FUBi_ratio": r,
            "asymmetry_factor": 2 * r - 1,
            "error_propagation": errors,
            "experimental_comparison": comparisons,
            "fubi_sweep_len": len(sweep),
            "primary_equations": [
                "m_BSM = (ℏω_SCm/c²)·S₂₆⁽³⁾·(2·F_{UBi}/F_U - 1)",
                f"ℏω_SCm = {self.delta_scm:.6e} J = {self.delta_scm / eV:.6e} eV",
                f"S₂₆⁽³⁾ = {S26_3RD:.6e}",
                f"2·F_{{UBi}}/F_U - 1 = {2 * r - 1:.4f}",
                f"m_BSM = {m:.6e} kg = {m_eV:.6e} eV/c²",
                f"σ_m/m = {errors['relative_error']:.4f} ({errors['relative_error'] * 100:.1f}%)",
            ],
        }


# ── §2  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("BSM Mass Bound Derivation from SCm Phonon Gap")
    print("=" * 72)

    bsm = BSMMassBound()
    result = bsm.compute({"FUBi_ratio": 0.6})

    for eq in result["primary_equations"]:
        print(f"  {eq}")

    print(f"\nError budget:")
    for k, v in result["error_propagation"]["contributions"].items():
        print(f"  {k}: {v * 100:.1f}%")

    print(f"\nExperimental comparisons:")
    for name, comp in result["experimental_comparison"].items():
        print(f"  {name}: {comp['source']}")

    print("\n✓ BSM mass bound derivation complete")
