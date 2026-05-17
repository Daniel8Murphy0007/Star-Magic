"""
_session282_alma_molecular_gas.py
=================================
Session 282 — Closes UQFF_CALIBRATION_AUDIT Gap #8.

ALMAMolecularGasCalculator: standard radiative-transfer / LTE conversion
of ALMA-observed integrated line intensities for CO, CS, HCN into
physical molecular-gas parameters (column density, H2 mass, dense-gas
mass, virial mass, brightness temperature).

Pure calculator — no hardcoded sources. Receives dataset from source2.cpp
or direct call; emits primary_equations / available_equations /
simulation_set per CondensedPhysics3.py MANDATORY ARCHITECTURE RULES.

Physics references
------------------
- Mangum & Shirley (2015), PASP 127, 266 — molecular column density formulae
- Bolatto, Wolfire & Leroy (2013), ARA&A 51, 207 — CO-to-H2 conversion factor
- Gao & Solomon (2004), ApJ 606, 271 — HCN dense-gas mass conversion
- Solomon et al. (1987), ApJ 319, 730 — virial mass for molecular clouds
"""

from __future__ import annotations
import math
from typing import Dict, Any, Optional

# ---------------------------------------------------------------------------
# Physical constants (SI)
# ---------------------------------------------------------------------------
K_B = 1.380649e-23          # J/K
H_PLANCK = 6.62607015e-34   # J s
C_LIGHT = 2.99792458e8      # m/s
M_SUN = 1.98892e30          # kg
PC = 3.0857e16              # m

# ---------------------------------------------------------------------------
# Reference transitions — frequencies (GHz), A-coefficients (s^-1),
# upper-state energies E_u/k (K), critical densities (cm^-3),
# partition-function approximation constants (Q ~ kT/hB at T >> T_rot)
# Values from CDMS / JPL / Mangum-Shirley Table 1.
# ---------------------------------------------------------------------------
TRANSITIONS = {
    "CO_10":  {"nu_GHz": 115.2712018, "A_ul": 7.203e-8,  "E_u_K": 5.5321, "g_u": 3.0, "B0_K": 2.766,  "n_crit_cm3": 2.2e3},
    "CO_21":  {"nu_GHz": 230.5380000, "A_ul": 6.910e-7,  "E_u_K": 16.596, "g_u": 5.0, "B0_K": 2.766,  "n_crit_cm3": 1.1e4},
    "CS_21":  {"nu_GHz": 97.98095330, "A_ul": 1.679e-5,  "E_u_K": 7.054,  "g_u": 5.0, "B0_K": 2.351,  "n_crit_cm3": 3.0e5},
    "HCN_10": {"nu_GHz": 88.63160230, "A_ul": 2.407e-5,  "E_u_K": 4.2536, "g_u": 9.0, "B0_K": 2.1284, "n_crit_cm3": 2.6e6},
}

# X-factors (CO -> H2; HCN -> dense H2) in M_sun / (K km s^-1 pc^2)
X_CO_MW = 4.3       # Bolatto+ 2013 Milky-Way average
X_HCN   = 10.0      # Gao & Solomon 2004 dense-gas

# Virial coefficient (Solomon+ 1987): M_vir [M_sun] = 1040 * sigma_v^2 [km/s] * R [pc]
K_VIR = 1040.0

# ---------------------------------------------------------------------------
# Core radiative-transfer helpers
# ---------------------------------------------------------------------------
def _hnu_kT(nu_Hz: float, T: float) -> float:
    return H_PLANCK * nu_Hz / (K_B * T)


def brightness_temperature(I_nu_Jy_per_sr: float, nu_Hz: float) -> float:
    """Rayleigh-Jeans T_B [K] from intensity [Jy/sr]."""
    # 1 Jy = 1e-26 W/m^2/Hz
    I_SI = I_nu_Jy_per_sr * 1e-26
    return (C_LIGHT ** 2 * I_SI) / (2.0 * K_B * nu_Hz ** 2)


def excitation_temperature_from_ratio(line_ratio: float, transition_hi: str,
                                      transition_lo: str) -> float:
    """T_ex [K] from integrated brightness-temperature line ratio (same species).

    Canonical optically-thin LTE form in K km/s units (e.g. Sakamoto+ 1994,
    Leroy+ 2009 for CO ladders):
        I_hi / I_lo = (nu_hi/nu_lo)^2 · exp(-(E_u_hi - E_u_lo)/T_ex)
    Solve:
        T_ex = (E_u_hi - E_u_lo) / ln[(nu_hi/nu_lo)^2 / R]
    Returns +inf if R exceeds the (nu_hi/nu_lo)^2 ceiling (super-thermal /
    optically thick / masing regime).
    """
    if transition_hi not in TRANSITIONS or transition_lo not in TRANSITIONS:
        raise ValueError(f"Unknown transition(s): {transition_hi}, {transition_lo}")
    thi = TRANSITIONS[transition_hi]
    tlo = TRANSITIONS[transition_lo]
    dE = thi["E_u_K"] - tlo["E_u_K"]
    if dE <= 0:
        raise ValueError("transition_hi must have higher E_u than transition_lo")
    pref = (thi["nu_GHz"] / tlo["nu_GHz"]) ** 2
    arg = pref / max(line_ratio, 1e-30)
    if arg <= 1.0:
        return float("inf")  # super-thermal / inverted
    return dE / math.log(arg)


def column_density_LTE(integrated_intensity_K_kms: float, transition: str,
                       T_ex: float) -> float:
    """LTE column density [cm^-2] from integrated brightness ∫T_B dv [K km/s].

    Mangum & Shirley (2015) Eq. 80, optically-thin limit, linear-rotor partition Q ≈ kT/(hB) + 1/3:
        N_tot = (8π k ν² / h c³ A_ul g_u) · Q(T_ex) · exp(E_u/T_ex)
                · (∫T_B dv) / (1 - exp(-hν/kT_ex))

    The numeric prefactor for CO(1-0), CS(2-1), HCN(1-0) reproduces
    canonical values (~10^15 N_X per K km/s at T_ex ≈ 30 K).
    """
    if transition not in TRANSITIONS:
        raise ValueError(f"Unknown transition: {transition}")
    tr = TRANSITIONS[transition]
    nu_Hz = tr["nu_GHz"] * 1e9
    A_ul = tr["A_ul"]
    g_u = tr["g_u"]
    E_u = tr["E_u_K"]
    B0 = tr["B0_K"]

    # Partition function (linear rotor, T_ex >> B0)
    Q = (T_ex / B0) + (1.0 / 3.0)

    # Prefactor in cgs to match Mangum-Shirley (returns cm^-2 per K km s^-1)
    # 8π k ν² / (h c³ A_ul g_u) — convert: use SI then translate K km/s -> SI (K m/s)
    prefactor_SI = (8.0 * math.pi * K_B * nu_Hz ** 2) / (H_PLANCK * (C_LIGHT ** 3) * A_ul * g_u)
    # ∫T_B dv: K km/s -> K m/s
    intI_SI = integrated_intensity_K_kms * 1.0e3
    # 1/(1 - exp(-hnu/kT))
    x = _hnu_kT(nu_Hz, T_ex)
    rj_corr = 1.0 / (1.0 - math.exp(-x)) if x < 50 else 1.0
    boltz = math.exp(E_u / T_ex)

    N_m2 = prefactor_SI * Q * boltz * rj_corr * intI_SI  # column in m^-2
    return N_m2 * 1.0e-4  # convert to cm^-2


def h2_mass_from_CO(L_CO_K_kms_pc2: float, alpha_CO: float = X_CO_MW) -> float:
    """M(H2) [M_sun] from CO luminosity L'_CO [K km/s pc^2] using X-factor."""
    return alpha_CO * L_CO_K_kms_pc2


def dense_gas_mass_from_HCN(L_HCN_K_kms_pc2: float, alpha_HCN: float = X_HCN) -> float:
    """M_dense [M_sun] from HCN(1-0) luminosity via Gao & Solomon 2004."""
    return alpha_HCN * L_HCN_K_kms_pc2


def virial_mass(sigma_v_kms: float, radius_pc: float) -> float:
    """Virial mass [M_sun] from line-of-sight velocity dispersion and cloud radius."""
    return K_VIR * (sigma_v_kms ** 2) * radius_pc


def line_luminosity_K_kms_pc2(integrated_intensity_K_kms: float, source_size_arcsec: float,
                              distance_pc: float) -> float:
    """L' [K km/s pc^2] = I_integrated · area_pc2.

    Source assumed top-hat with angular diameter `source_size_arcsec`.
    """
    theta_rad = source_size_arcsec * (math.pi / (180.0 * 3600.0))
    diameter_pc = theta_rad * distance_pc
    area_pc2 = math.pi * (diameter_pc / 2.0) ** 2
    return integrated_intensity_K_kms * area_pc2


# ---------------------------------------------------------------------------
# Calculator class (CondensedPhysics3.py-compliant interface)
# ---------------------------------------------------------------------------
class ALMAMolecularGasCalculator:
    """ALMA molecular-gas physical-parameter calculator.

    Accepted dataset keys (all optional unless noted):
        transition          : str — one of 'CO_10', 'CO_21', 'CS_21', 'HCN_10' (required)
        integrated_intensity_K_kms : float — ∫T_B dv (required)
        T_ex_K              : float — excitation temperature (default 30 K)
        line_ratio_hi_lo    : float — used with two-transition mode to derive T_ex
        transition_hi/lo    : str   — companion transitions for ratio mode
        source_size_arcsec  : float — angular size (for luminosity)
        distance_pc         : float — source distance (for luminosity / mass)
        sigma_v_kms         : float — velocity dispersion (for virial mass)
        radius_pc           : float — cloud radius (for virial mass)
        alpha_CO / alpha_HCN: float — override X-factor
        label               : str   — diagnostic label
    """

    def compute(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        transition = dataset.get("transition")
        if not transition or transition not in TRANSITIONS:
            raise ValueError(
                f"dataset['transition'] required; valid: {list(TRANSITIONS.keys())}")

        intI = float(dataset.get("integrated_intensity_K_kms", 0.0))
        if intI <= 0:
            raise ValueError("integrated_intensity_K_kms must be > 0")

        # Excitation temperature: either supplied, derived from ratio, or default
        T_ex: float = float(dataset.get("T_ex_K", 30.0))
        ratio = dataset.get("line_ratio_hi_lo")
        if ratio is not None and "transition_hi" in dataset and "transition_lo" in dataset:
            T_ex = excitation_temperature_from_ratio(
                float(ratio),
                str(dataset["transition_hi"]),
                str(dataset["transition_lo"]),
            )

        # Column density (always computable)
        N_col_cm2 = column_density_LTE(intI, transition, T_ex)

        # Brightness temperature peak (if continuum intensity supplied)
        T_B: Optional[float] = None
        if "I_nu_Jy_per_sr" in dataset:
            T_B = brightness_temperature(
                float(dataset["I_nu_Jy_per_sr"]),
                TRANSITIONS[transition]["nu_GHz"] * 1e9,
            )

        # Luminosity & mass (require source size + distance)
        L_prime: Optional[float] = None
        M_H2: Optional[float] = None
        M_dense: Optional[float] = None
        if "source_size_arcsec" in dataset and "distance_pc" in dataset:
            L_prime = line_luminosity_K_kms_pc2(
                intI,
                float(dataset["source_size_arcsec"]),
                float(dataset["distance_pc"]),
            )
            if transition.startswith("CO"):
                alpha = float(dataset.get("alpha_CO", X_CO_MW))
                M_H2 = h2_mass_from_CO(L_prime, alpha)
            elif transition == "HCN_10":
                alpha = float(dataset.get("alpha_HCN", X_HCN))
                M_dense = dense_gas_mass_from_HCN(L_prime, alpha)

        # Virial mass
        M_vir: Optional[float] = None
        if "sigma_v_kms" in dataset and "radius_pc" in dataset:
            M_vir = virial_mass(float(dataset["sigma_v_kms"]),
                                float(dataset["radius_pc"]))

        # Density regime (compare to n_crit)
        n_crit = TRANSITIONS[transition]["n_crit_cm3"]

        primary = [
            f"transition = {transition} (nu = {TRANSITIONS[transition]['nu_GHz']:.3f} GHz)",
            f"T_ex = {T_ex:.2f} K",
            f"N_col(LTE) = {N_col_cm2:.3e} cm^-2",
        ]
        if T_B is not None:
            primary.append(f"T_B(peak) = {T_B:.2f} K")
        if L_prime is not None:
            primary.append(f"L' = {L_prime:.3e} K km s^-1 pc^2")
        if M_H2 is not None:
            primary.append(f"M(H2) = alpha_CO · L' = {M_H2:.3e} M_sun")
        if M_dense is not None:
            primary.append(f"M_dense = alpha_HCN · L' = {M_dense:.3e} M_sun")
        if M_vir is not None:
            primary.append(f"M_vir = 1040·sigma_v^2·R = {M_vir:.3e} M_sun")

        available = [
            "N_col_LTE(intI, transition, T_ex) — Mangum-Shirley 2015 Eq. 80",
            "T_ex from line ratio (hi/lo, same species) — Boltzmann LTE",
            "M(H2) = alpha_CO · L'_CO  [Bolatto+ 2013, alpha=4.3 MW]",
            "M_dense = alpha_HCN · L'_HCN  [Gao & Solomon 2004, alpha=10]",
            "M_vir = 1040 · sigma_v^2 · R  [Solomon+ 1987]",
            "T_B = c^2 I_nu / (2 k nu^2)  [Rayleigh-Jeans]",
            f"Density regime: n_crit({transition}) = {n_crit:.1e} cm^-3",
        ]
        simulation = [
            {"step": "1. Estimate T_ex from line ratio or assume LTE",
             "needs": ["line_ratio_hi_lo", "transitions_hi_lo"]},
            {"step": "2. Compute N_col from ∫T_B dv",
             "function": "column_density_LTE"},
            {"step": "3. Convert L' to mass via X-factor",
             "needs": ["source_size_arcsec", "distance_pc"]},
            {"step": "4. Cross-check with virial mass",
             "needs": ["sigma_v_kms", "radius_pc"]},
        ]

        result: Dict[str, Any] = {
            "primary_equations": primary,
            "available_equations": available,
            "simulation_set": simulation,
            "transition": transition,
            "nu_GHz": TRANSITIONS[transition]["nu_GHz"],
            "T_ex_K": T_ex,
            "N_col_cm2": N_col_cm2,
            "n_crit_cm3": n_crit,
            "label": dataset.get("label", transition),
        }
        if T_B is not None:
            result["T_B_K"] = T_B
        if L_prime is not None:
            result["L_prime_K_kms_pc2"] = L_prime
        if M_H2 is not None:
            result["M_H2_Msun"] = M_H2
        if M_dense is not None:
            result["M_dense_Msun"] = M_dense
        if M_vir is not None:
            result["M_vir_Msun"] = M_vir
        return result


SESSION_282_CALCULATORS = {
    "ALMAMolecularGasCalculator": ALMAMolecularGasCalculator,
}

# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------
def _run_tests() -> int:
    pass_count = 0
    fail_count = 0

    def check(name: str, cond: bool, detail: str = "") -> None:
        nonlocal pass_count, fail_count
        if cond:
            pass_count += 1
            print(f"  [PASS] {name}  {detail}")
        else:
            fail_count += 1
            print(f"  [FAIL] {name}  {detail}")

    print("Session 282 — ALMA molecular-gas calculator smoke tests")
    print("-" * 70)

    calc = ALMAMolecularGasCalculator()

    # T-1: CO(1-0) Orion-like dense core
    r1 = calc.compute({
        "transition": "CO_10",
        "integrated_intensity_K_kms": 100.0,
        "T_ex_K": 30.0,
        "label": "Orion_core_proxy",
    })
    check("T-1 CO(1-0) N_col positive & finite",
          math.isfinite(r1["N_col_cm2"]) and r1["N_col_cm2"] > 0,
          f"N={r1['N_col_cm2']:.3e}")

    # T-2: order of magnitude — CO(1-0) at T=30K, intI=100 K km/s -> N ~ 10^18 cm^-2
    check("T-2 CO(1-0) N_col ~ 10^17–10^19 cm^-2 (canonical range)",
          1e17 <= r1["N_col_cm2"] <= 1e19,
          f"N={r1['N_col_cm2']:.3e}")

    # T-3: temperature scaling — higher T_ex -> larger N (for fixed intI in optically thin)
    r2 = calc.compute({
        "transition": "CO_10",
        "integrated_intensity_K_kms": 100.0,
        "T_ex_K": 60.0,
    })
    check("T-3 N_col(T=60) > N_col(T=30) (Q growth dominates)",
          r2["N_col_cm2"] > r1["N_col_cm2"],
          f"ratio={r2['N_col_cm2']/r1['N_col_cm2']:.3f}")

    # T-4: intensity linearity — double intI -> double N
    r3 = calc.compute({
        "transition": "CO_10",
        "integrated_intensity_K_kms": 200.0,
        "T_ex_K": 30.0,
    })
    check("T-4 N_col linear in integrated intensity",
          abs(r3["N_col_cm2"] / r1["N_col_cm2"] - 2.0) < 1e-9,
          f"ratio={r3['N_col_cm2']/r1['N_col_cm2']:.6f}")

    # T-5: T_ex from line ratio CO(2-1)/CO(1-0)=2.5 (typical warm GMC, R21~0.6 in T units;
    # ratio in K km/s units thermalizes to (nu_hi/nu_lo)^2 = 4 at high T).
    Tex_inferred = excitation_temperature_from_ratio(2.5, "CO_21", "CO_10")
    check("T-5 T_ex from CO(2-1)/(1-0)=2.5 in physical GMC range (10–100 K)",
          10.0 < Tex_inferred < 100.0,
          f"T_ex={Tex_inferred:.2f} K")

    # T-6: HCN(1-0) — dense-gas tracer; n_crit ~ 2.6e6
    r_hcn = calc.compute({
        "transition": "HCN_10",
        "integrated_intensity_K_kms": 10.0,
        "T_ex_K": 30.0,
    })
    check("T-6 HCN(1-0) n_crit matches dense-gas regime",
          r_hcn["n_crit_cm3"] >= 1e6,
          f"n_crit={r_hcn['n_crit_cm3']:.1e}")

    # T-7: CS(2-1) — dense-gas tracer (~3e5)
    r_cs = calc.compute({
        "transition": "CS_21",
        "integrated_intensity_K_kms": 5.0,
        "T_ex_K": 20.0,
    })
    check("T-7 CS(2-1) n_crit in 1e5–1e6 range",
          1e5 <= r_cs["n_crit_cm3"] <= 1e6,
          f"n_crit={r_cs['n_crit_cm3']:.1e}")

    # T-8: virial mass Solomon+ 1987 — sigma=2 km/s, R=1 pc -> M_vir = 4160 M_sun
    r_vir = calc.compute({
        "transition": "CO_10",
        "integrated_intensity_K_kms": 50.0,
        "T_ex_K": 25.0,
        "sigma_v_kms": 2.0,
        "radius_pc": 1.0,
    })
    check("T-8 M_vir(sigma=2,R=1) = 4160 M_sun",
          abs(r_vir["M_vir_Msun"] - 4160.0) < 1e-6,
          f"M_vir={r_vir['M_vir_Msun']:.1f}")

    # T-9: H2 mass via X-factor — L'=1e6 K km/s pc^2 -> M_H2 = 4.3e6 M_sun (MW)
    r_mass = calc.compute({
        "transition": "CO_10",
        "integrated_intensity_K_kms": 50.0,
        "T_ex_K": 25.0,
        "source_size_arcsec": 30.0,
        "distance_pc": 5e6,    # ~5 Mpc
    })
    expected_L = line_luminosity_K_kms_pc2(50.0, 30.0, 5e6)
    check("T-9 L' matches helper computation",
          abs(r_mass["L_prime_K_kms_pc2"] - expected_L) < 1e-6,
          f"L'={r_mass['L_prime_K_kms_pc2']:.3e}")
    check("T-10 M_H2 = alpha_CO · L'",
          abs(r_mass["M_H2_Msun"] - X_CO_MW * expected_L) < 1e-3,
          f"M_H2={r_mass['M_H2_Msun']:.3e}")

    # T-11: dense-gas mass via HCN
    r_dense = calc.compute({
        "transition": "HCN_10",
        "integrated_intensity_K_kms": 5.0,
        "T_ex_K": 30.0,
        "source_size_arcsec": 10.0,
        "distance_pc": 5e6,
    })
    check("T-11 M_dense from HCN populated",
          "M_dense_Msun" in r_dense and r_dense["M_dense_Msun"] > 0,
          f"M_dense={r_dense.get('M_dense_Msun', 0):.3e}")

    # T-12: brightness temperature (RJ) sanity
    T_B = brightness_temperature(1e6, 115e9)  # 1 MJy/sr at 115 GHz
    check("T-12 T_B finite, positive for 1 MJy/sr @ 115 GHz",
          math.isfinite(T_B) and T_B > 0,
          f"T_B={T_B:.3e} K")

    # T-13: error handling — bad transition
    try:
        calc.compute({"transition": "NOPE", "integrated_intensity_K_kms": 1.0})
        check("T-13 unknown transition raises", False, "no exception")
    except ValueError:
        check("T-13 unknown transition raises", True)

    # T-14: error handling — zero intensity
    try:
        calc.compute({"transition": "CO_10", "integrated_intensity_K_kms": 0.0})
        check("T-14 zero intensity raises", False, "no exception")
    except ValueError:
        check("T-14 zero intensity raises", True)

    # T-15: architecture compliance — three required output keys
    check("T-15 output has primary/available/simulation keys",
          all(k in r1 for k in
              ("primary_equations", "available_equations", "simulation_set")))

    # T-16: registration dict
    check("T-16 SESSION_282_CALCULATORS exports calculator",
          "ALMAMolecularGasCalculator" in SESSION_282_CALCULATORS)

    print("-" * 70)
    print(f"Results: {pass_count} PASS, {fail_count} FAIL")
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    import sys
    sys.exit(_run_tests())
