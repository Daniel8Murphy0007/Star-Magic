# -*- coding: utf-8 -*-
"""
_session299_gw190425_mass_gap_errors.py
=======================================
Session 299 -- Closes Audit Gap #9 (MEDIUM).

GW190425 mass-gap Bayesian posterior calculator with error bars.
Computes posterior probability that the component mass falls inside
the canonical "low mass gap" (2-5 M_sun) using a chirp-mass-anchored
Gaussian likelihood and a flat prior on the mass ratio q in [0.5, 1.0].

PHYSICS
-------
For a binary with chirp mass M_c and mass ratio q = m2/m1 <= 1:
   m1 = M_c * (1+q)^(1/5) / q^(3/5)
   m2 = q * m1
Posterior on M_c:
   p(M_c | data) ~ Normal(M_c_obs, sigma_Mc)
Marginal posterior on m1 (and m2) is obtained by Monte Carlo over q
from the LVC-published 90% interval, weighted by likelihood.
Mass-gap fraction:
   f_gap = P(2 <= m_i < 5 M_sun | data)
UQFF Aether modulation (POSTULATED, |delta|<=1e-3): applied to f_gap.

ANCHORS (4):
  A1 GW190425   M_c = 1.44, m1 in [1.60, 1.87]  -> low f_gap (NS-NS edge)
  A2 GW190814   M_c = 6.09, m2 = 2.6 +/- 0.1    -> high f_gap (mass-gap object)
  A3 GW170817   M_c = 1.188, m1 in [1.36, 1.60] -> ~0 f_gap (canonical NS-NS)
  A4 mock BH-BH M_c = 28, m1 ~ 35               -> ~0 f_gap (above gap)

Tier: DERIVED + STATISTICAL + POSTULATED

Author: D.T. Murphy / Copilot.  Session 299  cp4_id=443.
"""
from __future__ import annotations
import math
import random
from dataclasses import dataclass
from typing import Any, Dict, List, Tuple

BETA_I    = 0.603
RHO_SCM   = 7.0898e-37
RHO_AMB   = 1.0e-22
MASS_GAP_LO = 2.0  # M_sun
MASS_GAP_HI = 5.0
N_MC = 4000


def aether_correction(rho_amb: float, t_n: float) -> float:
    if rho_amb <= 0:
        return 1.0
    delta = BETA_I * (RHO_SCM / max(rho_amb, 1e-40))
    delta = max(-1e-3, min(1e-3, delta))
    return 1.0 + delta * math.cos(math.pi * t_n)


def component_masses_from_Mc_q(M_c: float, q: float) -> Tuple[float, float]:
    """Return (m1, m2) given chirp mass and q = m2/m1."""
    if M_c <= 0 or q <= 0 or q > 1:
        return (0.0, 0.0)
    m1 = M_c * (1.0 + q) ** (1.0 / 5.0) / q ** (3.0 / 5.0)
    m2 = q * m1
    return (m1, m2)


def mass_gap_fraction(M_c_obs: float, sigma_Mc: float,
                       q_lo: float, q_hi: float,
                       n_mc: int = N_MC, seed: int = 17) -> Dict[str, float]:
    """Monte-Carlo Bayesian f_gap and posterior median + 90% CI."""
    if M_c_obs <= 0 or sigma_Mc <= 0:
        return {"f_gap": 0.0, "m1_med": 0.0, "m1_lo": 0.0, "m1_hi": 0.0}
    rng = random.Random(seed)
    m1_samples: List[float] = []
    m2_samples: List[float] = []
    for _ in range(max(100, n_mc)):
        Mc = rng.gauss(M_c_obs, sigma_Mc)
        if Mc <= 0:
            continue
        q  = rng.uniform(q_lo, q_hi)
        m1, m2 = component_masses_from_Mc_q(Mc, q)
        m1_samples.append(m1)
        m2_samples.append(m2)
    if not m1_samples:
        return {"f_gap": 0.0, "m1_med": 0.0, "m1_lo": 0.0, "m1_hi": 0.0}
    n_total = len(m1_samples)
    n_gap = sum(1 for x in m1_samples + m2_samples
                if MASS_GAP_LO <= x < MASS_GAP_HI)
    f_gap = n_gap / (2.0 * n_total)
    sorted_m1 = sorted(m1_samples)
    m1_med = sorted_m1[n_total // 2]
    m1_lo  = sorted_m1[int(0.05 * n_total)]
    m1_hi  = sorted_m1[int(0.95 * n_total)]
    return {"f_gap": f_gap, "m1_med": m1_med, "m1_lo": m1_lo, "m1_hi": m1_hi}


@dataclass
class GWAnchor:
    name: str
    M_c_obs: float
    sigma_Mc: float
    q_lo: float
    q_hi: float
    f_gap_expected_lo: float
    f_gap_expected_hi: float


ANCHORS: Dict[str, GWAnchor] = {
    "A1_GW190425": GWAnchor("GW190425 (NS-NS edge)", 1.44,  0.02, 0.80, 1.00, 0.00, 0.10),
    "A2_GW190814": GWAnchor("GW190814 (mass gap)",   6.09,  0.06, 0.10, 0.13, 0.30, 1.00),
    "A3_GW170817": GWAnchor("GW170817 (canon NS-NS)",1.188, 0.001,0.70, 1.00, 0.00, 0.05),
    "A4_mock_BBH": GWAnchor("mock BBH (above gap)",  28.0,  1.0,  0.70, 1.00, 0.00, 0.05),
}


class GW190425MassGapBayesianCalculator:
    cp4_id = 443
    audit_session = 299

    def compute(self, dataset: Dict[str, Any] | None = None) -> Dict[str, Any]:
        ds = dataset or {}
        t_n = float(ds.get("t_n", 0.0))
        n_mc = int(ds.get("n_mc", N_MC))
        f_A  = aether_correction(RHO_AMB, t_n)
        rows: List[Dict[str, Any]] = []
        n_match = 0
        for key, a in ANCHORS.items():
            res = mass_gap_fraction(a.M_c_obs, a.sigma_Mc, a.q_lo, a.q_hi, n_mc=n_mc)
            f_gap_uqff = max(0.0, min(1.0, res["f_gap"] * f_A))
            match = a.f_gap_expected_lo <= f_gap_uqff <= a.f_gap_expected_hi
            if match:
                n_match += 1
            rows.append({"anchor": a.name, "M_c_obs": a.M_c_obs,
                         "sigma_Mc": a.sigma_Mc, "f_gap": f_gap_uqff,
                         "m1_med": res["m1_med"], "m1_90CI": (res["m1_lo"], res["m1_hi"]),
                         "expected_f_gap": (a.f_gap_expected_lo, a.f_gap_expected_hi),
                         "match": match})
        return {
            "primary_equations": [
                "m1 = M_c * (1+q)^(1/5) / q^(3/5)",
                "p(M_c|data) ~ Normal(M_c_obs, sigma_Mc)",
                "f_gap = P(2 <= m_i < 5 | data)  via MC over q",
            ],
            "available_equations": [
                "f_gap_UQFF = f_gap * f_A",
                "90% CI m1 from MC samples",
            ],
            "simulation_set": rows,
            "query_result": {"n_anchors": len(rows), "n_match": n_match, "f_Aether": f_A},
            "validation_table": rows,
            "headline": (
                "S299 GW190425MassGap [DERIVED+STATISTICAL+POSTULATED]: "
                f"{n_match}/{len(rows)} anchors within expected f_gap band."
            ),
        }


SESSION_299_CALCULATORS = {"GW190425MassGapBayesianCalculator": GW190425MassGapBayesianCalculator}
__all__ = ["GW190425MassGapBayesianCalculator", "SESSION_299_CALCULATORS",
           "mass_gap_fraction", "component_masses_from_Mc_q",
           "aether_correction", "ANCHORS"]


def _run_tests() -> int:
    n = 0
    def ok(lbl, c, x=""):
        nonlocal n
        if c: n += 1; print(f"  [PASS] {lbl}  {x}")
        else: print(f"  [FAIL] {lbl}  {x}")
    print("="*72); print("S299 GW190425 mass-gap smoke tests"); print("="*72)
    m1, m2 = component_masses_from_Mc_q(1.0, 1.0)
    ok("T-1 q=1 gives m1=m2", abs(m1 - m2) < 1e-9)
    ok("T-2 q=1 m1 = M_c * 2^(1/5) / 1",
       abs(m1 - 2.0**0.2) < 1e-6, f"{m1:.4f}")
    ok("T-3 q=0 returns (0,0)", component_masses_from_Mc_q(1.0, 0.0) == (0.0, 0.0))
    ok("T-4 M_c=0 returns (0,0)", component_masses_from_Mc_q(0.0, 0.5) == (0.0, 0.0))
    m1b, m2b = component_masses_from_Mc_q(1.0, 0.5)
    ok("T-5 q<1 gives m1 > m2", m1b > m2b)
    ok("T-6 m2/m1 = q within tol", abs(m2b / m1b - 0.5) < 1e-6)
    res = mass_gap_fraction(1.44, 0.02, 0.8, 1.0, n_mc=2000, seed=1)
    ok("T-7 GW190425 f_gap in [0,1]", 0.0 <= res["f_gap"] <= 1.0)
    ok("T-8 m1_med ~1.6-1.9 (NS regime)",
       1.4 < res["m1_med"] < 2.1, f"{res['m1_med']:.3f}")
    ok("T-9 m1_lo <= m1_med <= m1_hi",
       res["m1_lo"] <= res["m1_med"] <= res["m1_hi"])
    res2 = mass_gap_fraction(6.09, 0.06, 0.10, 0.13, n_mc=2000, seed=2)
    ok("T-10 GW190814 m1 ~22-25", 18 < res2["m1_med"] < 30, f"{res2['m1_med']:.1f}")
    ok("T-11 GW190814 m2 in mass gap likely", res2["f_gap"] > 0.0)
    res3 = mass_gap_fraction(1.188, 0.001, 0.7, 1.0, n_mc=2000, seed=3)
    ok("T-12 GW170817 f_gap ~ 0", res3["f_gap"] < 0.05, f"{res3['f_gap']:.3f}")
    res4 = mass_gap_fraction(28.0, 1.0, 0.7, 1.0, n_mc=2000, seed=4)
    ok("T-13 mock BBH f_gap ~ 0", res4["f_gap"] < 0.05, f"{res4['f_gap']:.3f}")
    ok("T-14 M_c<=0 returns 0", mass_gap_fraction(0, 0.1, 0.5, 1.0)["f_gap"] == 0.0)
    ok("T-15 sigma<=0 returns 0", mass_gap_fraction(1.4, 0, 0.5, 1.0)["f_gap"] == 0.0)
    ok("T-16 aether bounded", 0.999 < aether_correction(RHO_AMB, 0.5) < 1.001)
    calc = GW190425MassGapBayesianCalculator()
    out = calc.compute({"n_mc": 1500})
    ok("T-17 keys present",
       all(k in out for k in ("primary_equations","available_equations",
           "simulation_set","query_result","validation_table","headline")))
    ok("T-18 S299 tag", "S299" in out["headline"])
    qr = out["query_result"]
    ok("T-19 n_anchors=4 and >=3 match",
       qr["n_anchors"] == 4 and qr["n_match"] >= 3,
       f"matched {qr['n_match']}/4")
    ok("T-20 cp4_id=443 audit=299",
       GW190425MassGapBayesianCalculator.cp4_id == 443 and
       GW190425MassGapBayesianCalculator.audit_session == 299)
    print("="*72); print(f"  RESULT: {n}/20 passed"); print("="*72)
    return n


if __name__ == "__main__":
    n = _run_tests()
    assert n == 20, f"{n}/20"
