# -*- coding: utf-8 -*-
"""
_session288_universal_buoyancy_simultaneous_solver.py
=====================================================
QCalcGeom — Session 288 — Universal Buoyancy ⇌ Universal Gravity solver

Adds a new level on top of the existing QCalcGeom v2.0.0 FUBi/FUBii crossing
machinery: a fully simultaneous symbolic + numeric 5-DOF closure for the
Aether-UA Universal Buoyancy system as the Greek φύσις conceived it —
counter-balanced outward Aether pressure (F_U_Bi) against inward collapse
(F_U_Bi_i), with mass and habitable-zone radius both EMERGING from the same
algebraic crossing.

SIMULTANEOUS SYSTEM (5 equations, 5 unknowns)
---------------------------------------------
Unknowns:  r, T_eq, F_UBi, F_UBi_i, g_eff

 E1.  F_UBi   −     (ρ_UA / c²)·V_body·g_eff                = 0    (Archimedes; UA)
 E2.  F_UBi_i − β_i·(ρ_SCm/c²)·V_body·g_eff·(1 − F_TRZ)     = 0    (collapsing SCm)
 E3.  g_eff   − G·M_⋆ / r²                                  = 0    (Newton emergent)
 E4.  σ·T_eq⁴ − L_⋆·(1−A) / (16π·r²)                        = 0    (radiative)
 E5.  T_eq    − T_target                                    = 0    (HZ closure)

UNIVERSAL GRAVITY OUTPUT
------------------------
    g_U(r) = g_eff(r) · [1 − χ_UB]
    χ_UB   = K_UB · ρ_SCm / (c² · ρ_body)
    K_UB   = 10 − 9·β_i / 10           (G7 ‖ G9 ‖ F_TRZ algebraic closure)
    with  ρ_UA = 10·ρ_SCm  (G7),  F_TRZ = 1/10  (1/|SO(5)|),
          ρ_SCm = 4√π · 10⁻³⁷  (G9 axiom)

The (E1−E2) buoyancy-net force divided by F_grav is exactly χ_UB at any r;
this is proved symbolically in T-2 below and is *independent of r*.

USAGE
-----
    from _session288_universal_buoyancy_simultaneous_solver import (
        UniversalBuoyancySimultaneousSolver,
        solve_universal_buoyancy,
        chi_UB,
    )
    calc = UniversalBuoyancySimultaneousSolver()
    out  = calc.compute({
        "M_star_kg": 1.989e30, "L_star_W": 3.828e26,
        "M_body_kg": 5.972e24, "R_body_m": 6.371e6, "albedo": 0.30,
        "T_target_K": 288.0,
    })

Closes UQFF_CALIBRATION_AUDIT.md Gap-related deepening of QCalcGeom v2 → v2.1.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Dict, List

import sympy as sp


# ---------------------------------------------------------------------------
# SECTION 1 — CANONICAL CONSTANTS (frozen per /memories/repo/uqff_calibration_constants.md)
# ---------------------------------------------------------------------------
G_SI         = 6.6743e-11           # m^3 / (kg s^2)
C_LIGHT      = 2.998e8              # m / s
SIGMA_SB     = 5.670374419e-8       # W / (m^2 K^4)
PI           = math.pi
AU_M         = 1.49597870700e11     # m
L_SUN        = 3.828e26             # W
M_SUN        = 1.989e30             # kg

# UQFF canonical (G7 ‖ G9 ‖ F_TRZ — see Session 287)
RHO_SCm      = 4.0 * math.sqrt(PI) * 1.0e-37     # J / m^3   (G9 axiom)
RHO_UA       = 10.0 * RHO_SCm                     # J / m^3   (G7 companion)
F_TRZ        = 0.10                               # 1/|SO(5)|
BETA_I       = 0.6029                             # locked /uqff_program.py

# Universal-Buoyancy lump constant
K_UB         = 10.0 - 9.0 * BETA_I / 10.0         # = 9.45739


# ---------------------------------------------------------------------------
# SECTION 2 — SYMBOLIC SYSTEM (sympy) — derivation chain & χ_UB identity
# ---------------------------------------------------------------------------
def symbolic_system() -> Dict[str, Any]:
    """Return the symbolic 5-equation simultaneous system + derived χ_UB.

    Proves algebraically that
        (F_UBi - F_UBi_i) / F_grav  =  K_UB · ρ_SCm / (c² · ρ_body)
    independent of r and g_eff (cancellation is exact).
    """
    # symbols
    r, T_eq, F_UBi, F_UBi_i, g_eff = sp.symbols(
        "r T_eq F_UBi F_UBi_i g_eff", positive=True, real=True
    )
    M_star, L_star, M_body, R_body, A_alb, T_target = sp.symbols(
        "M_star L_star M_body R_body A T_target", positive=True, real=True
    )
    G_s, c_s, sigma_s = sp.symbols("G c sigma", positive=True, real=True)
    rho_SCm_s, rho_UA_s, beta_s, fTRZ_s = sp.symbols(
        "rho_SCm rho_UA beta_i F_TRZ", positive=True, real=True
    )

    V_body = sp.Rational(4, 3) * sp.pi * R_body**3
    rho_body = M_body / V_body

    # 5-equation system
    E1 = sp.Eq(F_UBi,  (rho_UA_s / c_s**2) * V_body * g_eff)
    E2 = sp.Eq(F_UBi_i, beta_s * (rho_SCm_s / c_s**2) * V_body * g_eff * (1 - fTRZ_s))
    E3 = sp.Eq(g_eff,  G_s * M_star / r**2)
    E4 = sp.Eq(sigma_s * T_eq**4, L_star * (1 - A_alb) / (16 * sp.pi * r**2))
    E5 = sp.Eq(T_eq, T_target)

    # χ_UB derivation: (F_UBi - F_UBi_i) / F_grav, F_grav = M_body·g_eff
    F_UB_net_expr = (rho_UA_s / c_s**2 - beta_s * (rho_SCm_s / c_s**2) * (1 - fTRZ_s)) * V_body * g_eff
    F_grav_expr   = M_body * g_eff
    chi_UB_sym    = sp.simplify(F_UB_net_expr / F_grav_expr)
    # canonical closed form with ρ_UA = 10·ρ_SCm, F_TRZ = 1/10:
    chi_UB_canon  = sp.simplify(
        chi_UB_sym.subs({rho_UA_s: 10 * rho_SCm_s, fTRZ_s: sp.Rational(1, 10)})
    )
    K_UB_sym = sp.simplify(chi_UB_canon * rho_body * c_s**2 / rho_SCm_s)

    # Solve E4∧E5 for r (habitable-zone radius from radiative balance)
    r_hz = sp.solve([E4, E5], (r, T_eq), dict=True)
    return {
        "equations": [E1, E2, E3, E4, E5],
        "chi_UB_general": chi_UB_sym,
        "chi_UB_canonical": chi_UB_canon,
        "K_UB_symbolic": K_UB_sym,
        "r_hz_symbolic": r_hz,
        "symbols": {
            "r": r, "T_eq": T_eq, "F_UBi": F_UBi, "F_UBi_i": F_UBi_i,
            "g_eff": g_eff, "M_star": M_star, "L_star": L_star,
            "M_body": M_body, "R_body": R_body, "A": A_alb,
            "T_target": T_target, "rho_SCm": rho_SCm_s, "rho_UA": rho_UA_s,
            "beta_i": beta_s, "F_TRZ": fTRZ_s,
        },
    }


# ---------------------------------------------------------------------------
# SECTION 3 — NUMERIC SOLVER
# ---------------------------------------------------------------------------
@dataclass
class UBSolution:
    r_hz_m:        float
    T_eq_K:        float
    F_UBi_N:       float
    F_UBi_i_N:     float
    g_eff_m_s2:    float
    g_universal_m_s2: float
    delta_g_aether_m_s2: float   # = g_eff * χ_UB  (explicit; survives float64)
    chi_UB:        float
    K_UB:          float
    F_grav_N:      float
    rho_body_kg_m3: float
    r_hz_AU:       float
    habitable:     bool
    primary_equations: List[str] = field(default_factory=list)


def chi_UB(rho_body_kg_m3: float,
           beta_i: float = BETA_I,
           rho_SCm: float = RHO_SCm,
           c: float = C_LIGHT) -> float:
    """χ_UB(ρ_body) = K_UB · ρ_SCm / (c² · ρ_body)  — r-independent.

    Algebraic identity:
        (F_UBi − F_UBi_i)/F_grav  ≡  χ_UB
    Proven in symbolic_system().
    """
    if rho_body_kg_m3 <= 0:
        raise ValueError("rho_body must be > 0")
    K = 10.0 - 9.0 * beta_i / 10.0
    return K * rho_SCm / (c * c * rho_body_kg_m3)


def habitable_radius_radiative(L_star_W: float,
                                albedo: float,
                                T_target_K: float) -> float:
    """Closed-form Eq4∧Eq5 solution:  r = sqrt( L·(1−A) / (16π·σ·T⁴) )."""
    if L_star_W <= 0 or T_target_K <= 0:
        raise ValueError("L_star and T_target must be > 0")
    if not (0.0 <= albedo < 1.0):
        raise ValueError("albedo must be in [0, 1)")
    return math.sqrt(L_star_W * (1.0 - albedo) / (16.0 * PI * SIGMA_SB * T_target_K**4))


def solve_universal_buoyancy(M_star_kg:   float,
                             L_star_W:    float,
                             M_body_kg:   float,
                             R_body_m:    float,
                             albedo:      float,
                             T_target_K:  float = 288.0,
                             beta_i:      float = BETA_I) -> UBSolution:
    """Solve the 5-equation simultaneous system explicitly.

    Strategy: E4∧E5 → r;  E3 → g_eff;  E1,E2 → F_UBi, F_UBi_i;
    χ_UB closed form → g_universal.
    """
    if M_star_kg <= 0 or M_body_kg <= 0 or R_body_m <= 0:
        raise ValueError("masses and radius must be > 0")
    if L_star_W <= 0:
        raise ValueError("L_star must be > 0")
    if T_target_K <= 0:
        raise ValueError("T_target must be > 0")
    if not (0.0 <= albedo < 1.0):
        raise ValueError("albedo must be in [0, 1)")

    # E4 ∧ E5
    r = habitable_radius_radiative(L_star_W, albedo, T_target_K)

    # E3
    g_eff = G_SI * M_star_kg / (r * r)

    # body
    V_body   = (4.0 / 3.0) * PI * R_body_m**3
    rho_body = M_body_kg / V_body

    # E1, E2
    F_UBi   = (RHO_UA / (C_LIGHT * C_LIGHT)) * V_body * g_eff
    F_UBi_i = beta_i * (RHO_SCm / (C_LIGHT * C_LIGHT)) * V_body * g_eff * (1.0 - F_TRZ)

    # χ_UB & Universal Gravity
    K = 10.0 - 9.0 * beta_i / 10.0
    chi = K * RHO_SCm / (C_LIGHT * C_LIGHT * rho_body)
    delta_g = g_eff * chi          # ~ 1e-58 m/s^2 for terrestrial; stored explicitly
    g_U = g_eff - delta_g          # numerically == g_eff in float64; algebraically distinct

    F_grav = M_body_kg * g_eff

    # habitable test: T_eq in [273.15, 373.15] K
    hab = (273.15 <= T_target_K <= 373.15)

    eqns = [
        f"E1: F_UBi   = (ρ_UA/c²)·V·g_eff           = {F_UBi:.6e} N",
        f"E2: F_UBi_i = β_i·(ρ_SCm/c²)·V·g·(1-F_TRZ) = {F_UBi_i:.6e} N",
        f"E3: g_eff   = G·M_*/r²                    = {g_eff:.6e} m/s² @ r={r/AU_M:.4f} AU",
        f"E4: σ·T⁴   = L_*·(1-A)/(16π·r²)           → r = {r:.6e} m",
        f"E5: T_eq    = T_target                    = {T_target_K:.2f} K",
        f"-- χ_UB = K_UB·ρ_SCm/(c²·ρ_body) = {chi:.6e}   (r-independent)",
        f"-- g_U  = g_eff·(1-χ_UB)         = {g_U:.6e} m/s²",
    ]
    return UBSolution(
        r_hz_m=r, T_eq_K=T_target_K,
        F_UBi_N=F_UBi, F_UBi_i_N=F_UBi_i,
        g_eff_m_s2=g_eff, g_universal_m_s2=g_U,
        delta_g_aether_m_s2=delta_g,
        chi_UB=chi, K_UB=K, F_grav_N=F_grav,
        rho_body_kg_m3=rho_body, r_hz_AU=r / AU_M,
        habitable=hab, primary_equations=eqns,
    )


# ---------------------------------------------------------------------------
# SECTION 4 — CALCULATOR CLASS (CP3 interface)
# ---------------------------------------------------------------------------
class UniversalBuoyancySimultaneousSolver:
    """CP3-compatible calculator.

    Solves the Universal-Buoyancy ⇌ Universal-Gravity simultaneous system
    (5 equations, 5 unknowns) and returns the Aether-derived effective
    gravity g_U at the radiative habitable-zone radius.
    """
    cp4_id = 432
    audit_session = 288
    description = (
        "Universal Buoyancy simultaneous equation solver: F_U_Bi vs F_U_Bi_i "
        "(Aether UA vacuum, Greek φύσις counter-balance) closed against "
        "Newton emergent gravity + Stefan-Boltzmann radiative habitable zone."
    )

    def compute(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        sol = solve_universal_buoyancy(
            M_star_kg  = float(dataset["M_star_kg"]),
            L_star_W   = float(dataset["L_star_W"]),
            M_body_kg  = float(dataset["M_body_kg"]),
            R_body_m   = float(dataset["R_body_m"]),
            albedo     = float(dataset["albedo"]),
            T_target_K = float(dataset.get("T_target_K", 288.0)),
            beta_i     = float(dataset.get("beta_i", BETA_I)),
        )

        # Available equations the solver can also expose
        available = [
            "F_U_Bi(r,t_n)  = ρ_UA·(4π/3)·r·c²·cos(πt_n)            (QCalcGeom v2)",
            "F_U_Bi_i(r,t_n)= -β_i·G·M²/r²·cos(πt_n)·orbit          (QCalcGeom v2)",
            "F_U      = ΣUg_i + Um − F_U_Bi + F_U_Bi_i              (canonical UQFF)",
            "g_U(r)   = g_eff(r)·(1 − χ_UB)                         (this module)",
            "χ_UB     = K_UB·ρ_SCm/(c²·ρ_body)                      (closed identity)",
            "K_UB     = 10 − 9·β_i/10                               (G7 + G9 + F_TRZ closure)",
        ]

        # Simulation set: scan a few interpretation points
        sim_set = [
            {"case": "0.5·T_target", "r_AU": habitable_radius_radiative(
                sol.r_hz_AU and dataset["L_star_W"], dataset["albedo"],
                0.5 * dataset.get("T_target_K", 288.0)) / AU_M},
            {"case": "T_target",     "r_AU": sol.r_hz_AU},
            {"case": "2·T_target",   "r_AU": habitable_radius_radiative(
                dataset["L_star_W"], dataset["albedo"],
                2.0 * dataset.get("T_target_K", 288.0)) / AU_M},
        ]

        return {
            "primary_equations": sol.primary_equations,
            "available_equations": available,
            "simulation_set": sim_set,
            "headline": {
                "r_hz_m":             sol.r_hz_m,
                "r_hz_AU":            sol.r_hz_AU,
                "T_eq_K":             sol.T_eq_K,
                "g_eff_m_s2":         sol.g_eff_m_s2,
                "g_universal_m_s2":   sol.g_universal_m_s2,
                "delta_g_aether_m_s2": sol.delta_g_aether_m_s2,
                "chi_UB":             sol.chi_UB,
                "K_UB":               sol.K_UB,
                "F_UBi_N":            sol.F_UBi_N,
                "F_UBi_i_N":          sol.F_UBi_i_N,
                "F_grav_N":           sol.F_grav_N,
                "rho_body_kg_m3":     sol.rho_body_kg_m3,
                "habitable":          sol.habitable,
            },
        }


SESSION_288_CALCULATORS: Dict[str, type] = {
    "UniversalBuoyancySimultaneousSolver": UniversalBuoyancySimultaneousSolver,
}


# ---------------------------------------------------------------------------
# SECTION 5 — SMOKE TESTS
# ---------------------------------------------------------------------------
def _run_tests() -> int:
    bar = "=" * 72
    print(bar)
    print("Session 288 — Universal Buoyancy simultaneous solver smoke tests")
    print(bar)
    n = 0

    def assert_(name: str, cond: bool, detail: str = "") -> None:
        nonlocal n
        if cond:
            n += 1
            print(f"  [PASS] {name}  {detail}")
        else:
            print(f"  [FAIL] {name}  {detail}")

    # ---- Symbolic identities -------------------------------------------------
    sysm = symbolic_system()

    # T-1: symbolic system has 5 equations
    assert_("T-1 symbolic system has 5 equations",
            len(sysm["equations"]) == 5,
            f"n={len(sysm['equations'])}")

    # T-2: χ_UB canonical reduces to K_UB·ρ_SCm/(c²·ρ_body)
    K_sym = sp.simplify(sysm["K_UB_symbolic"])
    K_num = K_sym.subs({sp.Symbol("beta_i", positive=True): sp.Rational(6029, 10000)})
    K_val = float(sp.nsimplify(K_num, rational=False))
    assert_("T-2 K_UB symbolic == 10 − 9·β_i/10",
            abs(K_val - K_UB) < 1e-12,
            f"sym={K_val:.6f}  num={K_UB:.6f}")

    # T-3: χ_UB is r-independent (no r symbol in simplified form)
    chi = sp.simplify(sysm["chi_UB_canonical"])
    r_sym = sysm["symbols"]["r"]
    assert_("T-3 χ_UB canonical is r-independent",
            r_sym not in chi.free_symbols,
            f"free={sorted(s.name for s in chi.free_symbols)}")

    # T-4: ρ_UA = 10·ρ_SCm closure
    assert_("T-4 ρ_UA / ρ_SCm = 10 (G7 companion)",
            abs(RHO_UA / RHO_SCm - 10.0) < 1e-12,
            f"ratio={RHO_UA/RHO_SCm:.12f}")

    # T-5: K_UB numeric value
    assert_("T-5 K_UB = 9.45739",
            abs(K_UB - 9.45739) < 1e-5,
            f"K_UB={K_UB:.6f}")

    # ---- Earth (Sun-Earth) ---------------------------------------------------
    earth = solve_universal_buoyancy(
        M_star_kg=M_SUN, L_star_W=L_SUN,
        M_body_kg=5.972e24, R_body_m=6.371e6, albedo=0.30,
        T_target_K=288.0,
    )
    # T-6: Earth HZ radius ≈ 1.0 AU (within albedo-corrected ~0.85 AU)
    assert_("T-6 Earth HZ radius in [0.7, 1.1] AU",
            0.70 <= earth.r_hz_AU <= 1.10,
            f"r_hz={earth.r_hz_AU:.4f} AU")

    # T-7: g_eff at Earth orbit ≈ G·M_sun/AU² ≈ 5.93e-3 m/s²
    g_expected = G_SI * M_SUN / (earth.r_hz_m ** 2)
    assert_("T-7 g_eff = G·M_sun/r²",
            abs(earth.g_eff_m_s2 - g_expected) / g_expected < 1e-12,
            f"g={earth.g_eff_m_s2:.6e}")

    # T-8: χ_UB for Earth (ρ≈5515 kg/m³)
    chi_earth = chi_UB(5515.0)
    assert_("T-8 χ_UB(Earth) ≈ 1.35e-56 (Aether under-buoys rocky matter)",
            5e-57 < chi_earth < 5e-56,
            f"χ={chi_earth:.3e}")

    # T-9: F_UBi > F_UBi_i (Aether-UA outward dominates SCm collapse)
    assert_("T-9 F_UBi > F_UBi_i (UA dominates SCm)",
            earth.F_UBi_N > earth.F_UBi_i_N,
            f"F_UBi/F_UBi_i = {earth.F_UBi_N/earth.F_UBi_i_N:.3f}")

    # T-10: F_UBi/F_UBi_i = ρ_UA / (β_i·ρ_SCm·(1−F_TRZ)) = 10/(β_i·0.9) ≈ 18.42
    ratio_expected = 10.0 / (BETA_I * 0.9)
    ratio_actual   = earth.F_UBi_N / earth.F_UBi_i_N
    assert_("T-10 F_UBi/F_UBi_i = 10/(β·9/10) ≈ 18.42",
            abs(ratio_actual - ratio_expected) / ratio_expected < 1e-10,
            f"ratio={ratio_actual:.6f}  expected={ratio_expected:.6f}")

    # T-11: habitable flag True for T_target=288 K
    assert_("T-11 habitable=True at T_target=288 K",
            earth.habitable, f"hab={earth.habitable}")

    # ---- Mars ---------------------------------------------------------------
    mars = solve_universal_buoyancy(
        M_star_kg=M_SUN, L_star_W=L_SUN,
        M_body_kg=6.4171e23, R_body_m=3.3895e6, albedo=0.25,
        T_target_K=210.0,   # Mars T_eq target
    )
    # T-12: Mars HZ radius ≈ 1.5 AU
    assert_("T-12 Mars HZ radius in [1.3, 1.7] AU for T=210 K",
            1.30 <= mars.r_hz_AU <= 1.70,
            f"r={mars.r_hz_AU:.4f} AU")
    # T-13: Mars not habitable at T=210 K
    assert_("T-13 habitable=False at T=210 K (Mars cold)",
            not mars.habitable, f"hab={mars.habitable}")

    # ---- Proxima Centauri b ------------------------------------------------
    proxb = solve_universal_buoyancy(
        M_star_kg=0.1221 * M_SUN, L_star_W=0.0017 * L_SUN,
        M_body_kg=1.27 * 5.972e24, R_body_m=1.08 * 6.371e6, albedo=0.30,
        T_target_K=234.0,
    )
    # T-14: Proxima b HZ ≈ 0.05 AU
    assert_("T-14 Proxima b HZ radius in [0.03, 0.08] AU",
            0.03 <= proxb.r_hz_AU <= 0.08,
            f"r={proxb.r_hz_AU:.4f} AU")

    # ---- K2-18b -------------------------------------------------------------
    k2_18 = solve_universal_buoyancy(
        M_star_kg=0.495 * M_SUN, L_star_W=0.023 * L_SUN,
        M_body_kg=8.63 * 5.972e24, R_body_m=2.61 * 6.371e6, albedo=0.30,
        T_target_K=255.0,
    )
    # T-15: K2-18b HZ ≈ 0.14 AU (super-Earth)
    assert_("T-15 K2-18b HZ radius in [0.1, 0.25] AU",
            0.10 <= k2_18.r_hz_AU <= 0.25,
            f"r={k2_18.r_hz_AU:.4f} AU")

    # ---- Universal Gravity correction ---------------------------------------
    # T-16: explicit Δg = g_eff·χ_UB stored (survives float64 underflow)
    ok16 = True
    for nm, s in [("Earth", earth), ("Mars", mars), ("Prox b", proxb), ("K2-18b", k2_18)]:
        expect = s.g_eff_m_s2 * s.chi_UB
        if expect <= 0 or abs(s.delta_g_aether_m_s2 - expect) / expect > 1e-12:
            assert_(f"T-16/{nm} Δg_aether == g_eff·χ_UB", False,
                    f"Δg={s.delta_g_aether_m_s2:.3e} expect={expect:.3e}")
            ok16 = False
            break
    if ok16:
        assert_("T-16 Δg_aether = g_eff·χ_UB for all bodies", True, "4/4")

    # T-17: Aether buoyancy correction is strictly positive (lifts mass)
    assert_("T-17 Δg_aether > 0 (Aether buoyancy lifts mass)",
            earth.delta_g_aether_m_s2 > 0.0,
            f"Δg={earth.delta_g_aether_m_s2:.3e} m/s²")

    # ---- Validation guards --------------------------------------------------
    raised = 0
    for bad in [
        {"M_star_kg": -1, "L_star_W": L_SUN, "M_body_kg": 1, "R_body_m": 1, "albedo": 0.3},
        {"M_star_kg": M_SUN, "L_star_W": -1, "M_body_kg": 1, "R_body_m": 1, "albedo": 0.3},
        {"M_star_kg": M_SUN, "L_star_W": L_SUN, "M_body_kg": -1, "R_body_m": 1, "albedo": 0.3},
        {"M_star_kg": M_SUN, "L_star_W": L_SUN, "M_body_kg": 1, "R_body_m": -1, "albedo": 0.3},
        {"M_star_kg": M_SUN, "L_star_W": L_SUN, "M_body_kg": 1, "R_body_m": 1, "albedo": 1.5},
        {"M_star_kg": M_SUN, "L_star_W": L_SUN, "M_body_kg": 1, "R_body_m": 1, "albedo": 0.3, "T_target_K": -1},
    ]:
        try:
            solve_universal_buoyancy(**bad)
        except ValueError:
            raised += 1
    assert_("T-18 invalid inputs raise ValueError", raised == 6, f"raised={raised}/6")

    # T-19 chi_UB rejects ρ<=0
    try:
        chi_UB(-1)
        assert_("T-19 chi_UB(<=0) raises ValueError", False, "no raise")
    except ValueError:
        assert_("T-19 chi_UB(<=0) raises ValueError", True, "")

    # T-20 calculator integration via compute()
    calc = UniversalBuoyancySimultaneousSolver()
    out = calc.compute({
        "M_star_kg": M_SUN, "L_star_W": L_SUN,
        "M_body_kg": 5.972e24, "R_body_m": 6.371e6, "albedo": 0.30,
        "T_target_K": 288.0,
    })
    keys_required = {"primary_equations", "available_equations", "simulation_set", "headline"}
    assert_("T-20 compute() returns required keys",
            keys_required.issubset(out.keys()),
            f"missing={keys_required - out.keys()}")

    # T-21 calculator registry
    assert_("T-21 SESSION_288_CALCULATORS exposes class",
            "UniversalBuoyancySimultaneousSolver" in SESSION_288_CALCULATORS, "")

    print("-" * 72)
    print(f"  RESULT: {n}/21 tests passed")
    print(bar)
    return n


if __name__ == "__main__":
    n = _run_tests()
    print()
    print("Universal Buoyancy simultaneous solver — headline (Sun-Earth, 288 K):")
    sol = solve_universal_buoyancy(
        M_star_kg=M_SUN, L_star_W=L_SUN,
        M_body_kg=5.972e24, R_body_m=6.371e6, albedo=0.30,
        T_target_K=288.0,
    )
    print(f"  r_hz        = {sol.r_hz_AU:.6f} AU   ({sol.r_hz_m:.6e} m)")
    print(f"  T_eq        = {sol.T_eq_K:.2f} K")
    print(f"  g_eff       = {sol.g_eff_m_s2:.6e} m/s^2     (Newton emergent)")
    print(f"  g_universal = {sol.g_universal_m_s2:.6e} m/s^2     (after Aether buoyancy)")
    print(f"  Δg_aether   = {sol.delta_g_aether_m_s2:.6e} m/s^2     (= g_eff · χ_UB)")
    print(f"  chi_UB      = {sol.chi_UB:.6e}                 (r-independent identity)")
    print(f"  K_UB        = {sol.K_UB:.6f}                   (10 - 9*beta_i/10)")
    print(f"  F_UBi       = {sol.F_UBi_N:.6e} N")
    print(f"  F_UBi_i     = {sol.F_UBi_i_N:.6e} N")
    print(f"  F_grav      = {sol.F_grav_N:.6e} N")
    print(f"  habitable   = {sol.habitable}")
    assert n == 21, f"expected 21/21, got {n}/21"
