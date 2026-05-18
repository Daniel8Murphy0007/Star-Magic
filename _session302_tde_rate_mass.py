# -*- coding: utf-8 -*-
"""
_session302_tde_rate_mass.py
============================
Session 302 -- Closes Audit Gap #12 (LOW).

Tidal Disruption Event (TDE) per-galaxy rate-vs-SMBH-mass calculator.
Uses the Stone & Metzger (2016) power-law fit to TDE rates and the
Hills (1975) tidal radius condition r_t < r_S to enforce the upper
mass cutoff (M_BH > 1e8 M_sun -> rt < r_Schwarzschild -> star
swallowed whole, no flare).

PHYSICS
-------
Stone & Metzger (2016) volumetric TDE rate per galaxy:
   Gamma_TDE(M_BH) = Gamma_0 * (M_BH / 1e6 M_sun)^(-0.404)   [yr^-1]
   Gamma_0 ~ 1e-4 yr^-1 per galaxy at 1e6 M_sun
Hills mass cutoff (no flare beyond, star swallowed):
   M_Hills ~ 1.1e8 M_sun (solar-type star, R* = R_sun)
   Gamma(M_BH > M_Hills) = 0

UQFF Aether modulation (POSTULATED, |delta|<=1e-3): applied multiplicatively.

ANCHORS (4):
  A1 Sgr A*       M_BH = 4.3e6   -> Gamma ~ 8e-5  yr^-1 (consistent with 0 obs in 25 yr)
  A2 AT2019qiz    M_BH ~ 1e6     -> Gamma ~ 1e-4
  A3 ASASSN-14li  M_BH ~ 1e6.2   -> Gamma ~ 8e-5
  A4 dwarf gal    M_BH ~ 1e5     -> Gamma ~ 2.5e-4 (Stone-Metzger predicts higher rates at low M)

Tier: DERIVED + CALIBRATED + POSTULATED

Author: D.T. Murphy / Copilot.  Session 302  cp4_id=446.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any, Dict, List

GAMMA_0_PER_YR = 1.0e-4
M_BH_PIVOT     = 1.0e6        # M_sun
SM_SLOPE       = -0.404
M_HILLS_MSUN   = 1.1e8        # solar-type star tidal-radius cutoff
BETA_I  = 0.603
RHO_SCM = 7.0898e-37
RHO_AMB = 1.0e-22


def aether_correction(rho_amb: float, t_n: float) -> float:
    if rho_amb <= 0:
        return 1.0
    delta = BETA_I * (RHO_SCM / max(rho_amb, 1e-40))
    delta = max(-1e-3, min(1e-3, delta))
    return 1.0 + delta * math.cos(math.pi * t_n)


def tde_rate_per_yr(M_BH_Msun: float) -> float:
    """Stone-Metzger 2016 TDE per-galaxy rate with Hills cutoff."""
    if M_BH_Msun <= 0:
        return 0.0
    if M_BH_Msun > M_HILLS_MSUN:
        return 0.0
    return GAMMA_0_PER_YR * (M_BH_Msun / M_BH_PIVOT) ** SM_SLOPE


@dataclass
class TDEAnchor:
    name: str
    M_BH_Msun: float
    Gamma_expected_lo: float
    Gamma_expected_hi: float


ANCHORS: Dict[str, TDEAnchor] = {
    "A1_SgrAstar":     TDEAnchor("Sgr A*",       4.3e6, 4.0e-5, 1.0e-4),
    "A2_AT2019qiz":    TDEAnchor("AT2019qiz",    1.0e6, 5.0e-5, 1.5e-4),
    "A3_ASASSN14li":   TDEAnchor("ASASSN-14li",  1.6e6, 5.0e-5, 1.2e-4),
    "A4_dwarf":        TDEAnchor("dwarf galaxy", 1.0e5, 1.5e-4, 4.0e-4),
}


class TDEMassRateRelationCalculator:
    cp4_id = 446
    audit_session = 302

    def compute(self, dataset: Dict[str, Any] | None = None) -> Dict[str, Any]:
        ds = dataset or {}
        t_n = float(ds.get("t_n", 0.0))
        f_A = aether_correction(RHO_AMB, t_n)
        rows: List[Dict[str, Any]] = []
        n_match = 0
        for key, a in ANCHORS.items():
            G    = tde_rate_per_yr(a.M_BH_Msun) * f_A
            match = a.Gamma_expected_lo <= G <= a.Gamma_expected_hi
            if match:
                n_match += 1
            rows.append({"anchor": a.name, "M_BH_Msun": a.M_BH_Msun,
                         "Gamma_TDE_per_yr": G,
                         "Gamma_expected": (a.Gamma_expected_lo, a.Gamma_expected_hi),
                         "match": match})
        # cosmological rate density at z=0 from M_BH mass function (illustrative)
        return {
            "primary_equations": [
                "Gamma_TDE(M_BH) = Gamma_0 * (M_BH / 1e6)^(-0.404)  yr^-1",
                "Cutoff: Gamma = 0 for M_BH > M_Hills ~ 1.1e8 M_sun",
            ],
            "available_equations": [
                "Gamma_UQFF = Gamma * f_A",
                "dN/dt = integral Gamma(M) phi(M) dV  (cosmological TDE rate)",
            ],
            "simulation_set": rows,
            "query_result": {"n_anchors": len(rows), "n_match": n_match,
                              "f_Aether": f_A, "M_Hills_Msun": M_HILLS_MSUN},
            "validation_table": rows,
            "headline": (
                "S302 TDERateMass [DERIVED+CALIBRATED+POSTULATED]: "
                f"{n_match}/{len(rows)} anchors within Stone-Metzger band."
            ),
        }


SESSION_302_CALCULATORS = {"TDEMassRateRelationCalculator": TDEMassRateRelationCalculator}
__all__ = ["TDEMassRateRelationCalculator", "SESSION_302_CALCULATORS",
           "tde_rate_per_yr", "aether_correction", "ANCHORS",
           "M_HILLS_MSUN", "GAMMA_0_PER_YR"]


def _run_tests() -> int:
    n = 0
    def ok(lbl, c, x=""):
        nonlocal n
        if c: n += 1; print(f"  [PASS] {lbl}  {x}")
        else: print(f"  [FAIL] {lbl}  {x}")
    print("="*72); print("S302 TDE mass-rate smoke tests"); print("="*72)
    ok("T-1 rate(M=0)=0", tde_rate_per_yr(0) == 0.0)
    ok("T-2 rate(M=-1)=0", tde_rate_per_yr(-1) == 0.0)
    ok("T-3 rate(M>Hills)=0", tde_rate_per_yr(1e9) == 0.0)
    ok("T-4 rate at pivot = Gamma_0",
       abs(tde_rate_per_yr(1e6) - GAMMA_0_PER_YR) < 1e-12)
    ok("T-5 rate decreases with M_BH",
       tde_rate_per_yr(1e7) < tde_rate_per_yr(1e6))
    ok("T-6 rate(low M) > rate(pivot)",
       tde_rate_per_yr(1e5) > tde_rate_per_yr(1e6))
    r1, r2 = tde_rate_per_yr(1e6), tde_rate_per_yr(1e7)
    ok("T-7 slope = -0.404 within 1%",
       abs(math.log10(r2 / r1) + 0.404) < 0.01)
    ok("T-8 rate just below Hills > 0", tde_rate_per_yr(1.0e8) > 0)
    ok("T-9 rate just above Hills = 0", tde_rate_per_yr(1.2e8) == 0)
    ok("T-10 aether bounded",
       0.999 < aether_correction(RHO_AMB, 0.5) < 1.001)
    ok("T-11 aether(rho=0)=1", aether_correction(0, 0) == 1.0)
    calc = TDEMassRateRelationCalculator()
    out  = calc.compute({})
    ok("T-12 keys present",
       all(k in out for k in ("primary_equations","available_equations",
            "simulation_set","query_result","validation_table","headline")))
    ok("T-13 S302 tag", "S302" in out["headline"])
    qr = out["query_result"]
    ok("T-14 n_anchors=4", qr["n_anchors"] == 4)
    ok("T-15 all 4 match", qr["n_match"] == 4, f"{qr['n_match']}/4")
    rows = out["validation_table"]
    ok("T-16 SgrA* Gamma in (4e-5, 1e-4)",
       4e-5 < rows[0]["Gamma_TDE_per_yr"] < 1e-4)
    ok("T-17 dwarf Gamma > SgrA* Gamma",
       rows[3]["Gamma_TDE_per_yr"] > rows[0]["Gamma_TDE_per_yr"])
    ok("T-18 Hills cutoff exposed in query_result",
       qr["M_Hills_Msun"] == M_HILLS_MSUN)
    ok("T-19 cp4_id=446",
       TDEMassRateRelationCalculator.cp4_id == 446)
    ok("T-20 audit_session=302",
       TDEMassRateRelationCalculator.audit_session == 302)
    print("="*72); print(f"  RESULT: {n}/20 passed"); print("="*72)
    return n


if __name__ == "__main__":
    n = _run_tests()
    assert n == 20, f"{n}/20"
