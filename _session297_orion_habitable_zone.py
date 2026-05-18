# -*- coding: utf-8 -*-
"""
_session297_orion_habitable_zone.py
===================================
Session 297 -- Closes Audit Gap #7 (MEDIUM).

Add Orion (the Trapezium / OB-association radiation field) to the
QCalcGeom habitable-zone test suite. Models the HZ around an embedded
solar-mass disk in the Orion Nebula where strong UV / X-ray flux from
the OB stars (theta1 Ori C in particular) modifies the standard HZ
boundaries via photoevaporation and UV photochemistry.

PHYSICS
-------
Standard Kopparapu (2013) HZ scaled by luminosity:
   d_HZ_in/out = sqrt(L_star / L_sun) * d_HZ_sun
External UV from OB cluster shifts the effective HZ inward
(photoevaporation makes outer disk uninhabitable but inner disk
shielded by gas):
   d_HZ_eff_out = d_HZ_out * min(1, (F_UV_sun / F_UV_cluster)^(1/4))
   F_UV_cluster = L_UV_OB / (4 pi r_cluster^2)
UQFF aether modulation (POSTULATED, |delta|<=1e-3):
   d_HZ_UQFF = d_HZ_eff * f_A

ANCHORS (4 host-star types embedded in Orion at r ~ 0.05 pc from theta1 OriC):
  A1 Solar analog G2V   : L=1 Lsun  -> classical HZ 0.95-1.37 AU,
                                       eff HZ_out compressed by ~0.5
  A2 K-dwarf K5V        : L=0.16 Lsun
  A3 M-dwarf M4V        : L=0.013 Lsun
  A4 isolated G2V (control, no Orion UV) : F_UV_cluster=0 -> standard HZ

Tier: DERIVED + CALIBRATED + POSTULATED

Author: D.T. Murphy / Copilot.  Session 297  cp4_id=441.
"""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Any, Dict, List

# constants
L_SUN_CGS  = 3.828e33      # erg/s
AU_CGS     = 1.496e13
PC_CGS     = 3.0857e18
F_UV_SUN_AT_1AU = 16.0     # erg/s/cm^2 (FUV/Habing normalization, ~10^4 G0)
BETA_I     = 0.603
RHO_SCM    = 7.0898e-37
RHO_AMB    = 1.0e-22       # Orion nebula molecular cloud density (g/cm^3)
HZ_IN_SUN  = 0.95          # AU
HZ_OUT_SUN = 1.37          # AU
L_UV_OB_THETA1 = 1.0e38    # erg/s (theta1 Ori C UV luminosity)
R_CLUSTER_PC = 0.05        # pc, distance from theta1 OriC


def aether_correction(rho_amb: float, t_n: float) -> float:
    if rho_amb <= 0:
        return 1.0
    delta = BETA_I * (RHO_SCM / max(rho_amb, 1e-40))
    delta = max(-1e-3, min(1e-3, delta))
    return 1.0 + delta * math.cos(math.pi * t_n)


def hz_classical_au(L_star_Lsun: float) -> tuple[float, float]:
    if L_star_Lsun <= 0:
        return (0.0, 0.0)
    s = math.sqrt(L_star_Lsun)
    return (HZ_IN_SUN * s, HZ_OUT_SUN * s)


def f_uv_cluster_at_r(L_UV_cgs: float, r_pc: float) -> float:
    """Local UV flux from cluster source, erg/s/cm^2."""
    if L_UV_cgs <= 0 or r_pc <= 0:
        return 0.0
    r_cm = r_pc * PC_CGS
    return L_UV_cgs / (4.0 * math.pi * r_cm * r_cm)


def hz_orion_modified_au(L_star_Lsun: float, F_UV_cluster: float
                          ) -> tuple[float, float]:
    hz_in, hz_out = hz_classical_au(L_star_Lsun)
    if F_UV_cluster <= 0:
        return (hz_in, hz_out)
    factor = min(1.0, (F_UV_SUN_AT_1AU / F_UV_cluster) ** 0.25)
    return (hz_in, hz_out * factor)


@dataclass
class HZAnchor:
    name: str
    L_star_Lsun: float
    F_UV_cluster: float
    expect_compressed: bool


F_UV_ORION = f_uv_cluster_at_r(L_UV_OB_THETA1, R_CLUSTER_PC)


ANCHORS: Dict[str, HZAnchor] = {
    "A1_G2V_orion":   HZAnchor("G2V in Orion (0.05 pc)",  1.0,    F_UV_ORION, True),
    "A2_K5V_orion":   HZAnchor("K5V in Orion (0.05 pc)",  0.16,   F_UV_ORION, True),
    "A3_M4V_orion":   HZAnchor("M4V in Orion (0.05 pc)",  0.013,  F_UV_ORION, True),
    "A4_G2V_iso":     HZAnchor("G2V isolated (control)",  1.0,    0.0,        False),
}


class OrionHabitableZoneCalculator:
    cp4_id = 441
    audit_session = 297

    def compute(self, dataset: Dict[str, Any] | None = None) -> Dict[str, Any]:
        ds = dataset or {}
        t_n = float(ds.get("t_n", 0.0))
        f_A = aether_correction(RHO_AMB, t_n)
        rows: List[Dict[str, Any]] = []
        n_compressed = 0
        for key, a in ANCHORS.items():
            hz_in_cl, hz_out_cl = hz_classical_au(a.L_star_Lsun)
            hz_in, hz_out = hz_orion_modified_au(a.L_star_Lsun, a.F_UV_cluster)
            hz_in_uqff = hz_in * f_A
            hz_out_uqff = hz_out * f_A
            compressed = hz_out_uqff < hz_out_cl * 0.99
            if compressed:
                n_compressed += 1
            rows.append({
                "anchor": a.name,
                "L_star_Lsun": a.L_star_Lsun,
                "F_UV_cluster": a.F_UV_cluster,
                "HZ_classical_AU": (hz_in_cl, hz_out_cl),
                "HZ_UQFF_AU": (hz_in_uqff, hz_out_uqff),
                "compressed": compressed,
                "expect_compressed": a.expect_compressed,
                "match": compressed == a.expect_compressed,
            })
        return {
            "primary_equations": [
                "d_HZ = sqrt(L/L_sun) * d_HZ_sun (Kopparapu 2013)",
                "F_UV_cluster = L_UV / (4 pi r^2)",
                "d_HZ_eff_out = d_HZ_out * min(1, (F_sun/F_cl)^(1/4))",
            ],
            "available_equations": ["d_HZ_UQFF = d_HZ_eff * f_A"],
            "simulation_set": rows,
            "query_result": {"n_anchors": len(rows),
                              "n_compressed": n_compressed,
                              "F_UV_orion_at_0p05pc": F_UV_ORION,
                              "f_Aether": f_A},
            "validation_table": rows,
            "headline": (
                "S297 OrionHZ [DERIVED+CALIBRATED+POSTULATED]: "
                f"{n_compressed}/{len(rows)} Orion-embedded systems compressed HZ."
            ),
        }


SESSION_297_CALCULATORS = {
    "OrionHabitableZoneCalculator": OrionHabitableZoneCalculator,
}

__all__ = [
    "OrionHabitableZoneCalculator", "SESSION_297_CALCULATORS",
    "hz_classical_au", "hz_orion_modified_au", "f_uv_cluster_at_r",
    "aether_correction", "ANCHORS", "F_UV_ORION",
]


def _run_tests() -> int:
    n = 0
    def ok(label, cond, extra=""):
        nonlocal n
        if cond:
            n += 1; print(f"  [PASS] {label}  {extra}")
        else:
            print(f"  [FAIL] {label}  {extra}")
    print("=" * 72)
    print("S297 OrionHZ smoke tests")
    print("=" * 72)

    hzin, hzout = hz_classical_au(1.0)
    ok("T-1 G2V classical HZ = 0.95-1.37 AU",
       abs(hzin - 0.95) < 1e-6 and abs(hzout - 1.37) < 1e-6)
    ok("T-2 HZ(L=0) = (0,0)", hz_classical_au(0) == (0.0, 0.0))
    ok("T-3 HZ scales as sqrt(L)",
       abs(hz_classical_au(4.0)[1] / hzout - 2.0) < 1e-6)
    ok("T-4 M-dwarf HZ closer in",
       hz_classical_au(0.013)[1] < hzout)
    f_uv = f_uv_cluster_at_r(L_UV_OB_THETA1, 0.05)
    ok("T-5 Orion UV flux > solar at 1 AU", f_uv > F_UV_SUN_AT_1AU,
       f"{f_uv:.3e}")
    ok("T-6 F_UV(L=0)=0", f_uv_cluster_at_r(0, 0.05) == 0.0)
    ok("T-7 F_UV(r=0)=0", f_uv_cluster_at_r(1e38, 0) == 0.0)
    ok("T-8 F_UV scales as 1/r^2",
       abs(f_uv_cluster_at_r(1e38, 0.1) /
           f_uv_cluster_at_r(1e38, 0.05) - 0.25) < 1e-6)
    hzin_o, hzout_o = hz_orion_modified_au(1.0, f_uv)
    ok("T-9 Orion-modified HZ_out < classical HZ_out", hzout_o < hzout,
       f"{hzout_o:.3f} < {hzout:.3f}")
    ok("T-10 Orion-modified HZ_in unchanged",
       abs(hzin_o - hzin) < 1e-6)
    ok("T-11 modified HZ(F=0) == classical",
       hz_orion_modified_au(1.0, 0.0) == (hzin, hzout))
    ok("T-12 aether bounded",
       0.999 < aether_correction(RHO_AMB, 0.5) < 1.001)
    ok("T-13 aether(rho=0)=1", aether_correction(0, 0) == 1.0)
    calc = OrionHabitableZoneCalculator()
    out = calc.compute({})
    ok("T-14 calculator returns required keys",
       all(k in out for k in ("primary_equations", "available_equations",
                                "simulation_set", "query_result",
                                "validation_table", "headline")))
    ok("T-15 headline tag", "S297" in out["headline"])
    qr = out["query_result"]
    ok("T-16 n_anchors=4", qr["n_anchors"] == 4)
    ok("T-17 3 Orion-embedded compressed",
       qr["n_compressed"] == 3, f"{qr['n_compressed']}/3")
    rows = out["validation_table"]
    ok("T-18 all anchor matches",
       all(r["match"] for r in rows),
       f"{sum(r['match'] for r in rows)}/{len(rows)}")
    ok("T-19 cp4_id=441",
       OrionHabitableZoneCalculator.cp4_id == 441)
    ok("T-20 audit_session=297",
       OrionHabitableZoneCalculator.audit_session == 297)

    print("=" * 72)
    print(f"  RESULT: {n}/20 passed")
    print("=" * 72)
    return n


if __name__ == "__main__":
    n = _run_tests()
    assert n == 20, f"{n}/20"
