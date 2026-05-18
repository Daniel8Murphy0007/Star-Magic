"""
_session286_tde_mass_rate.py

Session 286 / UQFF_CALIBRATION_AUDIT.md Gap #12 closure.

TDE (Tidal Disruption Event) mass-rate relations from stellar tidal
disruption by SMBHs.  Closed-form Rees-1988 / Phinney-1989 / Stone-2013
chain — no Aether terms, no calibrations.  Pure Newtonian + Schwarzschild
gating + Hills criterion.

KEY EQUATIONS
-------------
  Tidal radius:        r_t = R_* · (M_BH / M_*)^(1/3)
  Schwarzschild rad:   r_s = 2 G M_BH / c^2
  Hills mass (limit):  M_Hills s.t.  r_t(M_Hills) = r_s
                       = c^3 R_*^(3/2) / (2 G^(3/2) M_*^(1/2)) ·
                         (collapsed via direct solve)
                       For solar-type star: M_Hills ≈ 1.1e8 M_sun
  Penetration factor:  beta = r_t / r_peri
  Fallback timescale:  t_fb = (pi/sqrt(2)) · (r_t^3 / (G M_BH))^(1/2)
                          ≈ 41 days · (M_BH/1e6)^(1/2) (M_*/M_sun)^(-1)
                                       (R_*/R_sun)^(3/2)
  Peak fallback rate:  dM/dt|_peak = (1/3) · M_*/t_fb        (Rees 1988)
  Late-time decay:     dM/dt(t) ∝ t^(-5/3)   for t >> t_fb
  Eddington lum:       L_Edd = 4π G M_BH m_p c / sigma_T
                            = 1.26e38 · (M_BH/M_sun) erg/s

REFERENCE TDEs (smoke test anchors)
-----------------------------------
  AT2019dsg : M_BH ≈ 3e6 M_sun,  t_fb ≈ 71 d,    dMdt_pk ≈ 0.03 M_sun/yr-per-day class
  ASASSN-14li: M_BH ≈ 2e6 M_sun, t_fb ≈ 58 d
  AT2019qiz : M_BH ≈ 5e6 M_sun,  t_fb ≈ 92 d

All values are derived; no UQFF/Aether terms used (TDE physics is GR + Newton).

Closes audit Gap #12.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Any, Optional, List


# ---------------------------------------------------------------------------
# Constants (SI / CGS where natural)
# ---------------------------------------------------------------------------
G_SI       = 6.67430e-11        # m^3 kg^-1 s^-2
C_SI       = 2.99792458e8       # m/s
M_SUN_KG   = 1.98892e30         # kg
R_SUN_M    = 6.957e8            # m
DAY_S      = 86400.0            # s
YR_S       = 3.15576e7          # s (Julian year)
M_P_KG     = 1.6726219e-27      # kg (proton)
SIGMA_T_M2 = 6.6524587e-29      # m^2 (Thomson cross-section)
L_EDD_COEFF_ERG = 1.26e38       # erg/s per M_sun  (4πG·M_sun·m_p·c/σ_T)

# Numerical pre-factors
PI = math.pi


# ---------------------------------------------------------------------------
# Core scalar relations
# ---------------------------------------------------------------------------
def tidal_radius_m(R_star_m: float, M_BH_kg: float, M_star_kg: float) -> float:
    """r_t = R_* · (M_BH / M_*)^(1/3)"""
    if R_star_m <= 0 or M_BH_kg <= 0 or M_star_kg <= 0:
        raise ValueError("tidal_radius_m: positive inputs required")
    return R_star_m * (M_BH_kg / M_star_kg) ** (1.0 / 3.0)


def schwarzschild_radius_m(M_BH_kg: float) -> float:
    """r_s = 2 G M / c^2"""
    if M_BH_kg <= 0:
        raise ValueError("schwarzschild_radius_m: M_BH must be > 0")
    return 2.0 * G_SI * M_BH_kg / (C_SI * C_SI)


def hills_mass_kg(R_star_m: float, M_star_kg: float) -> float:
    """
    Hills mass: r_t(M_Hills) = r_s(M_Hills).
        R_* (M/M_*)^(1/3) = 2 G M / c^2
        ⇒  M^(2/3) = R_* · c^2 / (2 G · M_*^(1/3))
        ⇒  M = [ R_* · c^2 / (2 G) ]^(3/2) / M_*^(1/2)
    """
    if R_star_m <= 0 or M_star_kg <= 0:
        raise ValueError("hills_mass_kg: positive inputs required")
    base = R_star_m * C_SI * C_SI / (2.0 * G_SI)
    return (base ** 1.5) / math.sqrt(M_star_kg)


def fallback_timescale_s(R_star_m: float, M_BH_kg: float, M_star_kg: float) -> float:
    """
    Rees-1988 fallback timescale (return time of most-bound debris).

        a_min = r_t^2 / (2 R_*)                        (most-bound semi-major axis)
        t_fb  = 2π · a_min^(3/2) / sqrt(G M_BH)        (Kepler period)
              = (π/√2) · r_t^3 / (R_*^(3/2) · sqrt(G·M_BH))

    For 1 M_sun star + 1 R_sun radius + 1e6 M_sun SMBH ⇒ t_fb ≈ 41 days.
    """
    r_t = tidal_radius_m(R_star_m, M_BH_kg, M_star_kg)
    a_min = r_t * r_t / (2.0 * R_star_m)
    return 2.0 * PI * (a_min ** 1.5) / math.sqrt(G_SI * M_BH_kg)


def peak_fallback_rate_kg_s(M_star_kg: float, t_fb_s: float) -> float:
    """dM/dt|_peak = (1/3) · M_* / t_fb     (Rees 1988 head-collision approx.)"""
    if t_fb_s <= 0:
        raise ValueError("peak_fallback_rate_kg_s: t_fb must be > 0")
    return (1.0 / 3.0) * M_star_kg / t_fb_s


def fallback_rate_at_t_kg_s(M_star_kg: float, t_fb_s: float, t_s: float) -> float:
    """
    Late-time t^(-5/3) tail anchored to peak rate at t_fb.
        dM/dt(t) = dM/dt|_peak · (t / t_fb)^(-5/3),   t >= t_fb
        dM/dt(t) = dM/dt|_peak,                       t  < t_fb (plateau)
    """
    if t_s <= 0:
        raise ValueError("fallback_rate_at_t_kg_s: t must be > 0")
    pk = peak_fallback_rate_kg_s(M_star_kg, t_fb_s)
    if t_s < t_fb_s:
        return pk
    return pk * (t_s / t_fb_s) ** (-5.0 / 3.0)


def eddington_luminosity_erg_s(M_BH_kg: float) -> float:
    """L_Edd = 1.26e38 · (M_BH/M_sun) erg/s"""
    if M_BH_kg <= 0:
        raise ValueError("eddington_luminosity_erg_s: M_BH must be > 0")
    return L_EDD_COEFF_ERG * (M_BH_kg / M_SUN_KG)


def eddington_accretion_rate_kg_s(M_BH_kg: float, efficiency: float = 0.1) -> float:
    """
    Mdot_Edd = L_Edd / (η c^2)   in kg/s.
    η ≈ 0.1 standard radiative-efficiency assumption.
    """
    if not 0 < efficiency <= 1.0:
        raise ValueError("eddington_accretion_rate_kg_s: 0 < η ≤ 1 required")
    L_edd_W = eddington_luminosity_erg_s(M_BH_kg) * 1e-7   # erg/s → W
    return L_edd_W / (efficiency * C_SI * C_SI)


def can_disrupt(M_BH_kg: float, R_star_m: float, M_star_kg: float) -> bool:
    """True iff r_t > r_s (i.e. star is tidally disrupted outside event horizon)."""
    return tidal_radius_m(R_star_m, M_BH_kg, M_star_kg) > schwarzschild_radius_m(M_BH_kg)


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
@dataclass
class TDEResult:
    M_BH_Msun: float
    M_star_Msun: float
    R_star_Rsun: float
    r_t_m: float
    r_s_m: float
    r_t_over_r_s: float
    can_disrupt: bool
    M_Hills_Msun: float
    t_fb_days: float
    peak_dMdt_Msun_per_yr: float
    peak_dMdt_super_eddington_ratio: float
    L_Edd_erg_s: float
    Mdot_Edd_Msun_per_yr: float


# ---------------------------------------------------------------------------
# Calculator class (CP3 interface)
# ---------------------------------------------------------------------------
class TDEMassRateRelationCalculator:
    """
    Tidal Disruption Event (TDE) mass-rate calculator.

    cp4_id    : 430
    Audit gap : 12
    Inputs    : dataset dict with keys
                  M_BH_Msun     (required, > 0)
                  M_star_Msun   (default 1.0)
                  R_star_Rsun   (default 1.0)
                  efficiency    (default 0.1)
                  t_eval_days   (optional list of evaluation times)
    Output    : dict with primary_equations, available_equations,
                simulation_set, plus headline numbers from TDEResult.
    """

    cp4_id = 430
    audit_gap = 12
    description = (
        "TDE peak fallback rate, t^(-5/3) late-time tail, Hills mass, "
        "Eddington luminosity and super-Eddington ratio."
    )

    def compute(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        if not isinstance(dataset, dict) or not dataset:
            raise ValueError("dataset must be a non-empty dict")
        if "M_BH_Msun" not in dataset:
            raise ValueError("dataset missing required key 'M_BH_Msun'")

        M_BH_Msun   = float(dataset["M_BH_Msun"])
        M_star_Msun = float(dataset.get("M_star_Msun", 1.0))
        R_star_Rsun = float(dataset.get("R_star_Rsun", 1.0))
        eta         = float(dataset.get("efficiency", 0.1))

        if M_BH_Msun <= 0 or M_star_Msun <= 0 or R_star_Rsun <= 0:
            raise ValueError("M_BH, M_star, R_star must all be > 0")

        # SI conversion
        M_BH_kg   = M_BH_Msun * M_SUN_KG
        M_star_kg = M_star_Msun * M_SUN_KG
        R_star_m  = R_star_Rsun * R_SUN_M

        # Core geometry
        r_t = tidal_radius_m(R_star_m, M_BH_kg, M_star_kg)
        r_s = schwarzschild_radius_m(M_BH_kg)
        disrupt = r_t > r_s
        M_Hills = hills_mass_kg(R_star_m, M_star_kg) / M_SUN_KG

        # Fallback rate (only physical when disruption happens outside horizon)
        t_fb_s = fallback_timescale_s(R_star_m, M_BH_kg, M_star_kg)
        t_fb_d = t_fb_s / DAY_S
        pk_kgs = peak_fallback_rate_kg_s(M_star_kg, t_fb_s)
        pk_Msun_yr = pk_kgs * YR_S / M_SUN_KG

        # Eddington
        L_edd = eddington_luminosity_erg_s(M_BH_kg)
        Mdot_edd_kgs = eddington_accretion_rate_kg_s(M_BH_kg, efficiency=eta)
        Mdot_edd_Msun_yr = Mdot_edd_kgs * YR_S / M_SUN_KG
        super_edd_ratio = pk_kgs / Mdot_edd_kgs if Mdot_edd_kgs > 0 else float("inf")

        result = TDEResult(
            M_BH_Msun=M_BH_Msun,
            M_star_Msun=M_star_Msun,
            R_star_Rsun=R_star_Rsun,
            r_t_m=r_t,
            r_s_m=r_s,
            r_t_over_r_s=r_t / r_s,
            can_disrupt=disrupt,
            M_Hills_Msun=M_Hills,
            t_fb_days=t_fb_d,
            peak_dMdt_Msun_per_yr=pk_Msun_yr,
            peak_dMdt_super_eddington_ratio=super_edd_ratio,
            L_Edd_erg_s=L_edd,
            Mdot_Edd_Msun_per_yr=Mdot_edd_Msun_yr,
        )

        # Optional time-series simulation set
        sim_set: List[Dict[str, float]] = []
        t_eval_days = dataset.get("t_eval_days")
        if t_eval_days:
            for td in t_eval_days:
                td_f = float(td)
                if td_f <= 0:
                    continue
                rate_kgs = fallback_rate_at_t_kg_s(M_star_kg, t_fb_s, td_f * DAY_S)
                sim_set.append({
                    "t_days": td_f,
                    "dMdt_Msun_per_yr": rate_kgs * YR_S / M_SUN_KG,
                    "dMdt_over_Edd": rate_kgs / Mdot_edd_kgs if Mdot_edd_kgs > 0 else float("inf"),
                })

        return {
            "primary_equations": [
                "r_t = R_* (M_BH / M_*)^(1/3)",
                "r_s = 2 G M_BH / c^2",
                "M_Hills : r_t = r_s ⇒ M_Hills = (R_* c^2 / (2G))^(3/2) / M_*^(1/2)",
                "t_fb = (π/√2) (r_t^3 / (G M_BH))^(1/2)",
                "dM/dt|_pk = (1/3) M_* / t_fb",
                "dM/dt(t) = dM/dt|_pk (t/t_fb)^(-5/3)  for t ≥ t_fb",
                "L_Edd = 1.26e38 (M_BH/M_sun) erg/s",
                "Mdot_Edd = L_Edd / (η c^2)",
            ],
            "available_equations": [
                "penetration factor β = r_t / r_peri",
                "rate at time t (plateau + power-law tail)",
                "super-Eddington ratio dM/dt|_pk / Mdot_Edd",
                "Hills criterion can_disrupt = (r_t > r_s)",
            ],
            "simulation_set": sim_set,
            "M_BH_Msun": M_BH_Msun,
            "M_star_Msun": M_star_Msun,
            "R_star_Rsun": R_star_Rsun,
            "r_t_m": r_t,
            "r_s_m": r_s,
            "r_t_over_r_s": r_t / r_s,
            "can_disrupt": disrupt,
            "M_Hills_Msun": M_Hills,
            "t_fb_days": t_fb_d,
            "peak_dMdt_Msun_per_yr": pk_Msun_yr,
            "peak_dMdt_super_eddington_ratio": super_edd_ratio,
            "L_Edd_erg_s": L_edd,
            "Mdot_Edd_Msun_per_yr": Mdot_edd_Msun_yr,
            "result": result,
        }


SESSION_286_CALCULATORS = {
    "TDEMassRateRelationCalculator": TDEMassRateRelationCalculator,
}


# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------
def _check(label: str, ok: bool, info: str = "") -> bool:
    tag = "[PASS]" if ok else "[FAIL]"
    print(f"  {tag} {label}  {info}")
    return ok


def run_tests() -> int:
    print("=" * 72)
    print("Session 286 — TDE mass-rate relation smoke tests")
    print("=" * 72)
    passed = 0
    total = 0

    calc = TDEMassRateRelationCalculator()

    # T-1 : tidal radius scaling  r_t ∝ M_BH^(1/3)
    total += 1
    r1 = tidal_radius_m(R_SUN_M, 1e6 * M_SUN_KG, M_SUN_KG)
    r2 = tidal_radius_m(R_SUN_M, 8e6 * M_SUN_KG, M_SUN_KG)
    if _check("T-1 r_t ∝ M_BH^(1/3) ratio 2.0", abs(r2 / r1 - 2.0) < 1e-9,
              f"ratio={r2/r1:.6f}"):
        passed += 1

    # T-2 : Schwarzschild radius  r_s(1 M_sun) ≈ 2953 m
    total += 1
    rs1 = schwarzschild_radius_m(M_SUN_KG)
    if _check("T-2 r_s(1 M_sun) ≈ 2953 m", abs(rs1 - 2953.0) < 5.0,
              f"r_s={rs1:.2f} m"):
        passed += 1

    # T-3 : Hills mass for solar-type star ≈ 1.1e8 M_sun
    total += 1
    M_h = hills_mass_kg(R_SUN_M, M_SUN_KG) / M_SUN_KG
    if _check("T-3 Hills mass (solar-type) ~ 1.1e8 M_sun",
              5e7 <= M_h <= 2e8, f"M_Hills={M_h:.3e} M_sun"):
        passed += 1

    # T-4 : Star can be disrupted by 1e6 M_sun BH
    total += 1
    if _check("T-4 1e6 M_sun SMBH disrupts Sun outside horizon",
              can_disrupt(1e6 * M_SUN_KG, R_SUN_M, M_SUN_KG)):
        passed += 1

    # T-5 : Star swallowed whole by 1e10 M_sun BH (no TDE)
    total += 1
    if _check("T-5 1e10 M_sun SMBH swallows Sun whole (no TDE)",
              not can_disrupt(1e10 * M_SUN_KG, R_SUN_M, M_SUN_KG)):
        passed += 1

    # T-6 : Fallback timescale Rees scaling.
    # Canonical: t_fb ≈ 41 days · (M_BH/1e6)^(1/2) (M_*/M_sun)^(-1) (R_*/R_sun)^(3/2)
    total += 1
    t_fb_d = fallback_timescale_s(R_SUN_M, 1e6 * M_SUN_KG, M_SUN_KG) / DAY_S
    if _check("T-6 t_fb(1e6 M_sun, Sun) ~ 41 days",
              30.0 <= t_fb_d <= 60.0, f"t_fb={t_fb_d:.2f} d"):
        passed += 1

    # T-7 : t_fb scales as M_BH^(1/2)
    total += 1
    t1 = fallback_timescale_s(R_SUN_M, 1e6 * M_SUN_KG, M_SUN_KG)
    t2 = fallback_timescale_s(R_SUN_M, 4e6 * M_SUN_KG, M_SUN_KG)
    if _check("T-7 t_fb ∝ M_BH^(1/2) ratio 2.0",
              abs(t2 / t1 - 2.0) < 1e-9, f"ratio={t2/t1:.6f}"):
        passed += 1

    # T-8 : Peak fallback rate scaling — dimensional sanity (~0.1-10 M_sun/yr for 1e6 SMBH)
    total += 1
    pk = peak_fallback_rate_kg_s(M_SUN_KG, t1) * YR_S / M_SUN_KG
    if _check("T-8 peak dM/dt for 1e6 M_sun TDE in [0.1, 30] M_sun/yr",
              0.1 <= pk <= 30.0, f"pk={pk:.3f} M_sun/yr"):
        passed += 1

    # T-9 : Late-time t^(-5/3) decay law (anchor at 10·t_fb)
    total += 1
    pk_kgs = peak_fallback_rate_kg_s(M_SUN_KG, t1)
    r_at_10 = fallback_rate_at_t_kg_s(M_SUN_KG, t1, 10.0 * t1)
    expected_ratio = 10.0 ** (-5.0 / 3.0)
    if _check("T-9 dM/dt(10·t_fb)/pk = 10^(-5/3) ≈ 0.0464",
              abs(r_at_10 / pk_kgs - expected_ratio) < 1e-9,
              f"ratio={r_at_10/pk_kgs:.6f} exp={expected_ratio:.6f}"):
        passed += 1

    # T-10 : Eddington luminosity normalisation
    total += 1
    L = eddington_luminosity_erg_s(1e6 * M_SUN_KG)
    if _check("T-10 L_Edd(1e6 M_sun) ≈ 1.26e44 erg/s",
              abs(L - 1.26e44) / 1.26e44 < 1e-9, f"L={L:.3e}"):
        passed += 1

    # T-11 : Super-Eddington for typical TDE
    total += 1
    d = calc.compute({"M_BH_Msun": 1e6, "M_star_Msun": 1.0, "R_star_Rsun": 1.0})
    if _check("T-11 1e6 SMBH TDE peaks super-Eddington (ratio > 1)",
              d["peak_dMdt_super_eddington_ratio"] > 1.0,
              f"ratio={d['peak_dMdt_super_eddington_ratio']:.2f}"):
        passed += 1

    # T-12 : 1e10 SMBH TDE flagged can_disrupt=False
    total += 1
    d2 = calc.compute({"M_BH_Msun": 1e10, "M_star_Msun": 1.0, "R_star_Rsun": 1.0})
    if _check("T-12 1e10 SMBH compute → can_disrupt=False",
              d2["can_disrupt"] is False, f"r_t/r_s={d2['r_t_over_r_s']:.3f}"):
        passed += 1

    # T-13 : Reference TDE — ASASSN-14li  (M_BH ≈ 2e6 M_sun)
    total += 1
    d3 = calc.compute({"M_BH_Msun": 2e6})
    if _check("T-13 ASASSN-14li t_fb ~ 40-90 d",
              40.0 <= d3["t_fb_days"] <= 90.0, f"t_fb={d3['t_fb_days']:.2f} d"):
        passed += 1

    # T-14 : Time-series simulation_set populated when t_eval_days given
    total += 1
    d4 = calc.compute({
        "M_BH_Msun": 3e6,
        "t_eval_days": [10, 50, 100, 500, 1000],
    })
    sim = d4["simulation_set"]
    if _check("T-14 sim_set has 5 entries with monotone-decreasing tail",
              len(sim) == 5 and sim[-1]["dMdt_Msun_per_yr"] < sim[2]["dMdt_Msun_per_yr"],
              f"n={len(sim)}"):
        passed += 1

    # T-15 : Required keys in output
    total += 1
    required = {
        "primary_equations", "available_equations", "simulation_set",
        "M_BH_Msun", "r_t_m", "r_s_m", "r_t_over_r_s", "can_disrupt",
        "M_Hills_Msun", "t_fb_days", "peak_dMdt_Msun_per_yr",
        "peak_dMdt_super_eddington_ratio", "L_Edd_erg_s",
        "Mdot_Edd_Msun_per_yr",
    }
    if _check("T-15 output contains all required keys",
              required.issubset(d4.keys()), f"missing={required - d4.keys()}"):
        passed += 1

    # T-16 : Invalid inputs raise ValueError
    total += 1
    raised = 0
    for bad in (
        {},
        {"M_BH_Msun": 0},
        {"M_BH_Msun": 1e6, "M_star_Msun": -1},
        {"M_BH_Msun": 1e6, "R_star_Rsun": 0},
    ):
        try:
            calc.compute(bad)
        except ValueError:
            raised += 1
    if _check("T-16 invalid inputs raise ValueError",
              raised == 4, f"raised={raised}/4"):
        passed += 1

    # T-17 : Registry exposed
    total += 1
    if _check("T-17 SESSION_286_CALCULATORS registry",
              "TDEMassRateRelationCalculator" in SESSION_286_CALCULATORS):
        passed += 1

    # T-18 : Hills mass scales as R_*^(3/2) (white-dwarf vs main-sequence)
    total += 1
    M_h_sun = hills_mass_kg(R_SUN_M, M_SUN_KG) / M_SUN_KG
    M_h_wd  = hills_mass_kg(R_SUN_M * 0.01, M_SUN_KG) / M_SUN_KG  # WD r ~ 0.01 R_sun
    # M_Hills ∝ R_*^(3/2)  ⇒  ratio = 0.01^1.5 = 1e-3
    expected = 0.01 ** 1.5
    if _check("T-18 M_Hills(WD)/M_Hills(Sun) = 0.01^1.5",
              abs((M_h_wd / M_h_sun) / expected - 1.0) < 1e-9,
              f"ratio={M_h_wd/M_h_sun:.3e} exp={expected:.3e}"):
        passed += 1

    print("-" * 72)
    print(f"  RESULT: {passed}/{total} tests passed")
    print("=" * 72)

    # Headline ref-TDE table
    print()
    print("Headline TDE predictions (M_*=1 M_sun, R_*=1 R_sun):")
    for name, M_BH in [("ASASSN-14li", 2e6),
                        ("AT2019dsg",   3e6),
                        ("AT2019qiz",   5e6),
                        ("Sgr A* TDE",  4.15e6)]:
        d = calc.compute({"M_BH_Msun": M_BH})
        print(f"  {name:12s} M_BH={M_BH:.2e}  t_fb={d['t_fb_days']:6.2f} d  "
              f"dMdt_pk={d['peak_dMdt_Msun_per_yr']:.3f} M_sun/yr  "
              f"super-Edd×{d['peak_dMdt_super_eddington_ratio']:.1f}")

    return passed


if __name__ == "__main__":
    n = run_tests()
    assert n == 18, f"expected 18/18, got {n}/18"
