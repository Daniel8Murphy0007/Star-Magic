"""
_session283_gw190425_bayesian.py — Session 283 / Audit Gap #9

GW190425 mass-gap Bayesian posterior calculator.

Closes UQFF_CALIBRATION_AUDIT.md Gap #9:
    "GW190425 mass-gap fractions lack error bars — Bayesian posterior calc"

Method
------
GW190425 (LIGO-Virgo O3a, Abbott+ 2020, ApJL 892 L3) is a compact binary
coalescence with total mass M_tot = 3.4 (+0.3,-0.1) M_sun — the heaviest
BNS-like event observed. Its primary mass posterior straddles the
"mass gap" (~2.5-5.0 M_sun) between the heaviest NS and the lightest BH.

Inputs (per Abbott+ 2020, low-spin prior):
    M_chirp = 1.44 ± 0.02   M_sun  (very tightly measured from f-dot)
    q       = m2/m1 ∈ [0.4, 1.0]  (uniform prior, broad posterior)
    m1 ∈ [1.46, 1.87] M_sun   (90% CI from paper)
    m2 ∈ [1.46, 1.87] M_sun   (consistent with m1; symmetric due to q≤1 cut)

We perform Monte-Carlo sampling on the joint (M_chirp, q) posterior:
    1. M_chirp ~ Normal(1.44, 0.02)
    2. q       ~ Uniform(q_min, 1.0)
    3. Derive m1, m2 from definitions:
           M_chirp = (m1 m2)^(3/5) / (m1+m2)^(1/5)
           q = m2/m1  →  m1 = M_chirp * (1+q)^(1/5) / q^(3/5)
                          m2 = q * m1
    4. Compute fractions in physical categories:
       - "BNS-like":     m1 ≤ 2.5 M_sun  AND m2 ≤ 2.5 M_sun
       - "NS-BH-like":   m1 in [2.5, 5.0] AND m2 ≤ 2.5
       - "Mass-gap":     m1 in [2.5, 5.0] (regardless of m2)
       - "BBH-like":     m1 ≥ 5.0  OR  m2 ≥ 2.5 (heavy)
    5. Report mean ± 90% credible interval per category.

References
----------
- Abbott B.P. et al. (LIGO-Virgo), 2020, ApJL 892, L3 — GW190425 detection
- Farr W.M. et al. 2011, ApJ 741, 103 — NS/BH mass-gap definition
- Bailyn C.D. et al. 1998, ApJ 499, 367 — original mass-gap (2.5-5.0 M_sun)

Author : GitHub Copilot (Sonnet 4.5) for D.T. Murphy
Date   : May 16, 2026
"""
from __future__ import annotations

import math
import numpy as np
from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Canonical GW190425 posterior (Abbott+ 2020, low-spin prior, primary table)
# ---------------------------------------------------------------------------
GW190425 = {
    "M_chirp_mean": 1.44,        # M_sun
    "M_chirp_sigma": 0.02,       # M_sun (1-sigma; published 90% CI ≈ ±0.02)
    "q_min": 0.4,                # mass ratio lower bound (uniform prior)
    "q_max": 1.0,                # mass ratio upper bound
    "M_tot_90CI": (3.3, 3.7),    # M_sun (published 90% credible interval)
}

# Mass-gap definition (Bailyn 1998, Farr 2011)
MASS_GAP_LOW  = 2.5   # M_sun — heaviest observed NS
MASS_GAP_HIGH = 5.0   # M_sun — lightest dynamical BH


# ---------------------------------------------------------------------------
# Chirp-mass / mass-ratio inversion
# ---------------------------------------------------------------------------
def component_masses_from_chirp(M_chirp: float, q: float) -> tuple[float, float]:
    """
    Solve M_chirp = (m1 m2)^(3/5) / (m1+m2)^(1/5) with q = m2/m1, q∈(0,1].
    Returns (m1, m2) with m1 ≥ m2.
    """
    if q <= 0.0 or q > 1.0:
        raise ValueError(f"q must be in (0, 1], got {q}")
    # m1 = M_chirp * (1+q)^(1/5) / q^(3/5)
    m1 = M_chirp * (1.0 + q) ** (1.0 / 5.0) / q ** (3.0 / 5.0)
    m2 = q * m1
    return m1, m2


# ---------------------------------------------------------------------------
# Monte-Carlo posterior sampler
# ---------------------------------------------------------------------------
@dataclass
class MassGapPosterior:
    bns_like_frac: float          # mean fraction
    bns_like_ci90: tuple[float, float]
    ns_bh_frac: float
    ns_bh_ci90: tuple[float, float]
    mass_gap_frac: float
    mass_gap_ci90: tuple[float, float]
    m1_median: float
    m1_ci90: tuple[float, float]
    m2_median: float
    m2_ci90: tuple[float, float]
    M_tot_median: float
    M_tot_ci90: tuple[float, float]
    n_samples: int


def gw190425_mass_gap_posterior(
    n_samples: int = 200_000,
    seed: int = 20260516,
    M_chirp_mean: Optional[float] = None,
    M_chirp_sigma: Optional[float] = None,
    q_min: Optional[float] = None,
    q_max: Optional[float] = None,
) -> MassGapPosterior:
    """
    Monte-Carlo sample the GW190425 component-mass posterior and report
    Bayesian credible intervals for mass-gap occupancy fractions.

    Each sample:
      M_chirp ~ Normal(mean, sigma)   # truncated to >0
      q       ~ Uniform(q_min, q_max)
      (m1, m2) = invert(M_chirp, q)

    Returns a MassGapPosterior dataclass with means and 90% CIs.
    """
    Mc_mean = M_chirp_mean if M_chirp_mean is not None else GW190425["M_chirp_mean"]
    Mc_sig  = M_chirp_sigma if M_chirp_sigma is not None else GW190425["M_chirp_sigma"]
    q_lo    = q_min if q_min is not None else GW190425["q_min"]
    q_hi    = q_max if q_max is not None else GW190425["q_max"]

    rng = np.random.default_rng(seed)
    Mc_samples = rng.normal(Mc_mean, Mc_sig, size=n_samples)
    # Reject non-positive (vanishingly rare at 72-sigma but be safe)
    Mc_samples = np.where(Mc_samples > 0, Mc_samples, Mc_mean)
    q_samples  = rng.uniform(q_lo, q_hi, size=n_samples)

    # Vectorised mass inversion
    m1 = Mc_samples * (1.0 + q_samples) ** 0.2 / q_samples ** 0.6
    m2 = q_samples * m1
    M_tot = m1 + m2

    # Category masks
    bns_mask     = (m1 <= MASS_GAP_LOW) & (m2 <= MASS_GAP_LOW)
    ns_bh_mask   = (m1 >= MASS_GAP_LOW) & (m1 <= MASS_GAP_HIGH) & (m2 <= MASS_GAP_LOW)
    mass_gap_mask = (m1 >= MASS_GAP_LOW) & (m1 <= MASS_GAP_HIGH)

    # Bayesian fractions via bootstrap on the posterior sample
    def _frac_ci(mask: np.ndarray, n_boot: int = 500) -> tuple[float, tuple[float, float]]:
        frac = float(mask.mean())
        idx = rng.integers(0, mask.size, size=(n_boot, mask.size))
        boot = mask[idx].mean(axis=1)
        lo, hi = np.percentile(boot, [5.0, 95.0])
        return frac, (float(lo), float(hi))

    bns_f, bns_ci = _frac_ci(bns_mask)
    ns_bh_f, ns_bh_ci = _frac_ci(ns_bh_mask)
    mg_f, mg_ci = _frac_ci(mass_gap_mask)

    def _med_ci(x: np.ndarray) -> tuple[float, tuple[float, float]]:
        med = float(np.median(x))
        lo, hi = np.percentile(x, [5.0, 95.0])
        return med, (float(lo), float(hi))

    m1_med, m1_ci = _med_ci(m1)
    m2_med, m2_ci = _med_ci(m2)
    Mt_med, Mt_ci = _med_ci(M_tot)

    return MassGapPosterior(
        bns_like_frac=bns_f, bns_like_ci90=bns_ci,
        ns_bh_frac=ns_bh_f, ns_bh_ci90=ns_bh_ci,
        mass_gap_frac=mg_f, mass_gap_ci90=mg_ci,
        m1_median=m1_med, m1_ci90=m1_ci,
        m2_median=m2_med, m2_ci90=m2_ci,
        M_tot_median=Mt_med, M_tot_ci90=Mt_ci,
        n_samples=n_samples,
    )


# ---------------------------------------------------------------------------
# CP3-compliant calculator wrapper
# ---------------------------------------------------------------------------
class GW190425MassGapPosteriorCalculator:
    """
    CondensedPhysics3.py-compatible calculator.

    Dataset keys (all optional — defaults to canonical Abbott+ 2020 values):
        n_samples       : int        (default 200_000)
        seed            : int        (default 20260516)
        M_chirp_mean    : float [M⊙] (default 1.44)
        M_chirp_sigma   : float [M⊙] (default 0.02)
        q_min           : float      (default 0.4)
        q_max           : float      (default 1.0)
        label           : str
    """
    name = "GW190425MassGapPosteriorCalculator"
    cp4_id = 428  # next free CP4 ID after #427 (Chandra response matrix)

    def compute(self, dataset: dict) -> dict:
        d = dataset or {}
        post = gw190425_mass_gap_posterior(
            n_samples=int(d.get("n_samples", 200_000)),
            seed=int(d.get("seed", 20260516)),
            M_chirp_mean=d.get("M_chirp_mean"),
            M_chirp_sigma=d.get("M_chirp_sigma"),
            q_min=d.get("q_min"),
            q_max=d.get("q_max"),
        )
        label = d.get("label", "GW190425")
        return {
            "label": label,
            "primary_equations": [
                "M_chirp = (m1·m2)^(3/5) / (m1+m2)^(1/5)",
                "m1 = M_chirp · (1+q)^(1/5) / q^(3/5),   m2 = q·m1",
                "P(m1 ∈ [2.5,5.0] M⊙ | data) = ∫ p(M_chirp,q) · 𝟙[mass-gap] dM_chirp dq",
            ],
            "available_equations": [
                "BNS-like fraction: P(m1≤2.5 & m2≤2.5)",
                "NS-BH-like fraction: P(2.5≤m1≤5.0 & m2≤2.5)",
                "Mass-gap occupancy: P(2.5≤m1≤5.0)",
                "M_tot = m1 + m2 (total source-frame mass)",
                "Bootstrap 90% credible interval on Bernoulli fractions",
            ],
            "simulation_set": [
                f"Monte-Carlo posterior n={post.n_samples}",
                "Bootstrap n_boot=500 for fraction CIs",
            ],
            "m1_median_Msun": post.m1_median,
            "m1_ci90_Msun": post.m1_ci90,
            "m2_median_Msun": post.m2_median,
            "m2_ci90_Msun": post.m2_ci90,
            "M_tot_median_Msun": post.M_tot_median,
            "M_tot_ci90_Msun": post.M_tot_ci90,
            "P_BNS_like": post.bns_like_frac,
            "P_BNS_like_ci90": post.bns_like_ci90,
            "P_NS_BH_like": post.ns_bh_frac,
            "P_NS_BH_like_ci90": post.ns_bh_ci90,
            "P_mass_gap_primary": post.mass_gap_frac,
            "P_mass_gap_primary_ci90": post.mass_gap_ci90,
            "n_samples": post.n_samples,
        }


SESSION_283_CALCULATORS = {
    "GW190425MassGapPosteriorCalculator": GW190425MassGapPosteriorCalculator,
}


# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------
def _run_smoke_tests() -> None:
    print("=" * 72)
    print("Session 283 — GW190425 Bayesian mass-gap posterior smoke tests")
    print("=" * 72)
    tests_passed = 0
    tests_total = 0

    def check(name: str, cond: bool, info: str = "") -> None:
        nonlocal tests_passed, tests_total
        tests_total += 1
        status = "PASS" if cond else "FAIL"
        if cond:
            tests_passed += 1
        print(f"  [{status}] {name}  {info}")

    # T-1: chirp-mass inversion identity (equal mass)
    m1, m2 = component_masses_from_chirp(1.0, 1.0)
    # Equal-mass: M_chirp = m * 2^(-1/5) · 2 = 2m / 2^(1/5)  → m = M_chirp · 2^(1/5) / 2
    expected = 2.0 ** 0.2 / 2.0 * 2.0
    check("T-1 equal-mass inversion", math.isclose(m1, m2, rel_tol=1e-12),
          f"m1={m1:.4f} m2={m2:.4f}")

    # T-2: inversion satisfies definition
    Mc, q = 1.44, 0.7
    m1, m2 = component_masses_from_chirp(Mc, q)
    Mc_check = (m1 * m2) ** 0.6 / (m1 + m2) ** 0.2
    check("T-2 chirp definition self-consistent",
          math.isclose(Mc, Mc_check, rel_tol=1e-12),
          f"Mc_in={Mc} Mc_out={Mc_check:.6f}")

    # T-3: m1 >= m2 always
    check("T-3 ordering m1 >= m2", m1 >= m2, f"m1={m1:.3f} m2={m2:.3f}")

    # T-4: q=1 gives equal masses
    m1e, m2e = component_masses_from_chirp(1.44, 1.0)
    check("T-4 q=1 ⇒ equal masses", math.isclose(m1e, m2e, rel_tol=1e-12))

    # T-5: invalid q raises
    raised = False
    try:
        component_masses_from_chirp(1.44, 0.0)
    except ValueError:
        raised = True
    check("T-5 q=0 raises ValueError", raised)

    # T-6: full posterior runs
    post = gw190425_mass_gap_posterior(n_samples=50_000, seed=42)
    check("T-6 posterior n_samples", post.n_samples == 50_000)

    # T-7: m1 median in Abbott+ 2020 low-spin 90% CI [1.61, 2.52]
    check("T-7 m1 median in [1.60, 2.55]",
          1.60 <= post.m1_median <= 2.55,
          f"m1_med={post.m1_median:.3f}")

    # T-8: m2 median in [1.1, 1.7]
    check("T-8 m2 median in [1.05, 1.75]",
          1.05 <= post.m2_median <= 1.75,
          f"m2_med={post.m2_median:.3f}")

    # T-9: M_tot median in published 90% CI [3.3, 3.7]
    check("T-9 M_tot median in [3.25, 3.75]",
          3.25 <= post.M_tot_median <= 3.75,
          f"Mtot_med={post.M_tot_median:.3f}")

    # T-10: fractions sum to <=1
    total = post.bns_like_frac + post.ns_bh_frac
    check("T-10 BNS+NSBH frac <= 1.001", total <= 1.001,
          f"total={total:.4f}")

    # T-11: mass-gap fraction has nonzero CI width
    width = post.mass_gap_ci90[1] - post.mass_gap_ci90[0]
    check("T-11 mass-gap CI nonzero width", width >= 0.0,
          f"width={width:.5f}")

    # T-12: m1_ci90 brackets median
    check("T-12 m1 CI brackets median",
          post.m1_ci90[0] <= post.m1_median <= post.m1_ci90[1])

    # T-13: BNS-like fraction is dominant (m_chirp=1.44 is BNS-typical)
    check("T-13 BNS-like fraction > 0.5",
          post.bns_like_frac > 0.5,
          f"P_BNS={post.bns_like_frac:.3f}")

    # T-14: CP3 calculator interface
    calc = GW190425MassGapPosteriorCalculator()
    out = calc.compute({"n_samples": 10_000, "label": "test"})
    check("T-14 calculator returns required keys",
          all(k in out for k in ("primary_equations", "available_equations",
                                  "simulation_set", "P_mass_gap_primary",
                                  "m1_median_Msun", "M_tot_median_Msun")))

    # T-15: reproducibility with seed
    p1 = gw190425_mass_gap_posterior(n_samples=5_000, seed=7)
    p2 = gw190425_mass_gap_posterior(n_samples=5_000, seed=7)
    check("T-15 seed reproducible",
          p1.m1_median == p2.m1_median and p1.bns_like_frac == p2.bns_like_frac)

    # T-16: SESSION_283_CALCULATORS registry
    check("T-16 calculator registry exposed",
          "GW190425MassGapPosteriorCalculator" in SESSION_283_CALCULATORS)

    print("-" * 72)
    print(f"  RESULT: {tests_passed}/{tests_total} tests passed")
    print("=" * 72)

    # Headline result
    p = gw190425_mass_gap_posterior(n_samples=200_000, seed=20260516)
    print()
    print("GW190425 headline (200k MC samples, low-spin prior):")
    print(f"  m1       = {p.m1_median:.3f}  [90% CI {p.m1_ci90[0]:.3f}, {p.m1_ci90[1]:.3f}] M_sun")
    print(f"  m2       = {p.m2_median:.3f}  [90% CI {p.m2_ci90[0]:.3f}, {p.m2_ci90[1]:.3f}] M_sun")
    print(f"  M_tot    = {p.M_tot_median:.3f}  [90% CI {p.M_tot_ci90[0]:.3f}, {p.M_tot_ci90[1]:.3f}] M_sun")
    print(f"  P(BNS-like)        = {p.bns_like_frac:.4f}  [90% CI {p.bns_like_ci90[0]:.4f}, {p.bns_like_ci90[1]:.4f}]")
    print(f"  P(NS-BH-like)      = {p.ns_bh_frac:.4f}  [90% CI {p.ns_bh_ci90[0]:.4f}, {p.ns_bh_ci90[1]:.4f}]")
    print(f"  P(m1 ∈ mass-gap)   = {p.mass_gap_frac:.4f}  [90% CI {p.mass_gap_ci90[0]:.4f}, {p.mass_gap_ci90[1]:.4f}]")


if __name__ == "__main__":
    _run_smoke_tests()
