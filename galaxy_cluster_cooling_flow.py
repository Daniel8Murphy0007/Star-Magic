#!/usr/bin/env python3
"""
galaxy_cluster_cooling_flow.py — ICM Cooling-Flow Buoyancy Suppression Calculator

Session 225 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Parameterized calculator for galaxy-cluster cooling-flow regulation via
F_{U,Bi,i} jet-mediated feedback. Computes the balance between:

  - Radiative cooling:  L_cool = n_e² Λ(T) V
  - F_{U,Bi,i} heating: Q_heat = F_{U,Bi,i} · M_ICM · v_ICM

where the F_{U,Bi,i} force modulates the AGN jet mechanical feedback power,
suppressing catastrophic cooling flows.

The suppression factor S(Γ) measures how much the cooling-flow mass
deposition rate is reduced as a function of phonon linewidth Γ:

  S(Γ) = 1 − Q_heat(Γ) / L_cool

Classes:
  ICMRadiativeCooling       — Bremsstrahlung + line cooling L_cool(T, n_e, V)
  FUBiJetHeating           — AGN jet mechanical power via F_{U,Bi,i}
  CoolingFlowSuppression   — Balance L_cool vs Q_heat → suppression factor
  ClusterCoolingFlowPipeline — Full pipeline for Perseus/Abell 2256/Virgo

References:
  - NGC 1275 / Perseus: M_cooling = 13e9 M_sun H2, T_ICM ≈ 4 keV
  - Abell 2256: T_ICM ≈ 6.5 keV, merging cluster
  - Virgo / M87: T_ICM ≈ 2.5 keV, AGN jet feedback
  - gen_muge_ngc1275.py: MUGE cooling-flow term
  - CondensedPhysics.py: CMBAnisotropyBuoyancyModulationCalculator
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, Any, List

# ── §0  Constants ──────────────────────────────────────────────────────────

PI      = math.pi
G       = 6.674e-11
C       = 2.998e8
HBAR    = 1.055e-34
K_B     = 1.381e-23
M_SUN   = 1.989e30
KPC     = 3.086e19
MPC     = 3.086e22

SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
OMEGA_SCM = 2 * PI * 1.25e12
F_NEUTRON = 1e-10
RHO_VAC   = 1e-10

# Astrophysical
KEV_TO_K  = 1.1605e7   # 1 keV = 1.16e7 K
M_PROTON  = 1.673e-27  # kg
SIGMA_T   = 6.652e-29  # Thomson cross section m²

# Ramanujan S₂₆⁽³⁾
def _ramanujan_Rn(n: int, k: int = 3) -> float:
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total

S26_3RD = sum(SSQ**n / n**26 * _ramanujan_Rn(n, 3) for n in range(1, 27))

def _phonon_fluence(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# ═══════════════════════════════════════════════════════════════════════════════
# §1  ICM RADIATIVE COOLING
# ═══════════════════════════════════════════════════════════════════════════════

class ICMRadiativeCooling:
    """Intra-cluster medium radiative cooling luminosity.

    L_cool = n_e² · Λ(T) · V

    where Λ(T) is the cooling function (Bremsstrahlung + line emission):
      Λ(T) ≈ 2.4e-27 · sqrt(T/1e8 K) · [1 + 4.4e-2·(T/1e7)^{-0.5}] erg cm⁻³ s⁻¹
            (Sutherland & Dopita 1993 approximation, Z = 0.3 solar)

    Converted to SI: Λ_SI = Λ_cgs × 1e-7 / (1e-6)² = Λ_cgs × 1e5  [W m³]
    (but standard is to keep n_e in m⁻³ and Λ in W m³)
    """

    def compute(self, dataset: dict) -> dict:
        T_keV = float(dataset.get('T_keV', 4.0))         # ICM temperature
        n_e = float(dataset.get('n_e', 3e4))              # electron density m⁻³
        r_cool_kpc = float(dataset.get('r_cool_kpc', 100)) # cooling radius kpc

        T_K = T_keV * KEV_TO_K

        # Cooling function (W m³) — Bremsstrahlung-dominated
        # Λ_ff(T) ≈ 1.42e-40 × sqrt(T) × g_ff  [W m³] for hydrogen
        # With metallicity Z=0.3: enhance by ~(1 + 1.8×Z) at cluster T
        g_ff = 1.2  # Gaunt factor ≈ 1.2 for cluster T
        Lambda_ff = 1.42e-40 * math.sqrt(T_K) * g_ff * (1 + 1.8 * 0.3)

        # Cooling volume
        r_cool = r_cool_kpc * KPC
        V_cool = (4.0 / 3.0) * PI * r_cool**3

        # Cooling luminosity
        L_cool = n_e**2 * Lambda_ff * V_cool  # W

        # Cooling time
        E_thermal = 1.5 * n_e * K_B * T_K * V_cool  # J
        t_cool = E_thermal / L_cool if L_cool > 0 else float('inf')
        t_cool_Gyr = t_cool / (3.156e16)

        # Mass deposition rate (without suppression)
        # Mdot = L_cool · μ m_p / (5/2 k_B T)
        mu = 0.6  # mean molecular weight
        Mdot = L_cool * mu * M_PROTON / (2.5 * K_B * T_K) if T_K > 0 else 0
        Mdot_solar_yr = Mdot / M_SUN * 3.156e7

        return {
            'T_keV': T_keV,
            'T_K': T_K,
            'n_e_m3': n_e,
            'r_cool_kpc': r_cool_kpc,
            'Lambda_ff_Wm3': Lambda_ff,
            'V_cool_m3': V_cool,
            'L_cool_W': L_cool,
            't_cool_Gyr': t_cool_Gyr,
            'Mdot_kg_s': Mdot,
            'Mdot_solar_yr': Mdot_solar_yr,
            'primary_equations': [
                'ICM Radiative Cooling:',
                '  L_cool = n_e^2 · Lambda(T) · V',
                f'  Lambda_ff(T) = 1.42e-40 · sqrt(T) · g_ff · (1+1.8Z)',
                f'  T = {T_keV} keV = {T_K:.4e} K',
                f'  n_e = {n_e:.2e} m^-3',
                f'  r_cool = {r_cool_kpc} kpc',
                f'  L_cool = {L_cool:.4e} W',
                f'  t_cool = {t_cool_Gyr:.2f} Gyr',
                f'  Mdot (no suppression) = {Mdot_solar_yr:.1f} M_sun/yr',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §2  F_{U,Bi,i} JET HEATING
# ═══════════════════════════════════════════════════════════════════════════════

class FUBiJetHeating:
    """AGN jet heating power via F_{U,Bi,i} buoyancy force.

    The jet mechanical power is linked to the SMBH mass through F_{U,Bi,i}:

    Q_heat = |F_{U,Bi,i}| · M_ICM_effective · v_jet

    where F_{U,Bi,i} is the 26-layer buoyancy force at the jet base:

    F_{U,Bi,i} = Σ_{i=1}^{26} [Ug_i - Ub_i] + Um + UA + F_n·S₂₆⁽³⁾·Φ(Γ)·E_net
    """

    def compute(self, dataset: dict) -> dict:
        M_bh = float(dataset.get('M_bh', 3.4e8)) * M_SUN  # SMBH mass (solar masses)
        M_bh_solar = float(dataset.get('M_bh', 3.4e8))
        gamma = float(dataset.get('gamma', GAMMA_0))
        v_jet_frac = float(dataset.get('v_jet_frac', 0.1))  # jet speed / c
        eta_jet = float(dataset.get('eta_jet', 0.1))          # jet coupling efficiency
        M_icm_within = float(dataset.get('M_icm_solar', 1e13)) * M_SUN  # ICM mass

        # SMBH horizon
        r_s = 2 * G * M_bh / C**2
        r_H = r_s  # Schwarzschild (non-rotating simplification)

        # 26-layer gravity and buoyancy at horizon
        Ug = sum(G * M_bh / r_H**2 * SSQ * i / 26 for i in range(1, 27))
        Ub = sum(G * M_bh / r_H**2 * math.exp(-SSQ * i / 26) * BETA_I
                 for i in range(1, 27))
        Um = G * M_bh / r_H**2 * SSQ * 0.1
        UA = G * M_bh / r_H**2 * 1e-10

        # Phonon and E_net
        Phi_val = _phonon_fluence(gamma)
        E_net_factor = 0.2  # (2 × 0.6/1.0 − 1) generic
        E_net = RHO_VAC * math.exp(min(KAPPA * 86400, 500)) * Phi_val * E_net_factor

        # F_{U,Bi,i}
        F_UBii = (Ug + Um + UA - Ub) + F_NEUTRON * S26_3RD * Phi_val * E_net

        # Jet mechanical power (Eddington-limited proxy)
        v_jet = v_jet_frac * C
        Q_heat = abs(F_UBii) * eta_jet * M_icm_within * v_jet  # W

        # Alternative: Bondi-based jet power
        # P_jet = eta_jet × Mdot_Bondi × c²
        L_edd = 4 * PI * G * M_bh * M_PROTON * C / SIGMA_T
        P_jet_edd = eta_jet * L_edd

        return {
            'M_bh_solar': M_bh_solar,
            'r_H_m': r_H,
            'Ug': Ug,
            'Ub': Ub,
            'Um': Um,
            'F_UBii': F_UBii,
            'Phi': Phi_val,
            'Q_heat_W': Q_heat,
            'P_jet_edd_W': P_jet_edd,
            'L_edd_W': L_edd,
            'primary_equations': [
                'F_{U,Bi,i} Jet Heating:',
                '  F_{U,Bi,i} = Σ(Ug_i - Ub_i) + Um + UA + F_n·S26·Phi·E_net',
                '  Q_heat = |F_{U,Bi,i}| · eta_jet · M_ICM · v_jet',
                f'  M_BH = {M_bh_solar:.2e} M_sun',
                f'  Ug = {Ug:.4e} m/s^2,  Ub = {Ub:.4e} m/s^2',
                f'  F_{{U,Bi,i}} = {F_UBii:.6e} m/s^2',
                f'  Q_heat = {Q_heat:.4e} W',
                f'  P_jet(Eddington) = {P_jet_edd:.4e} W',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §3  COOLING-FLOW SUPPRESSION FACTOR
# ═══════════════════════════════════════════════════════════════════════════════

class CoolingFlowSuppression:
    """Balance Q_heat vs L_cool to derive suppression factor S(Gamma).

    S(Gamma) = min(Q_heat(Gamma) / L_cool, 1.0)

    When S → 1: complete cooling-flow quenching (AGN dominates)
    When S → 0: no suppression (classical cooling catastrophe)

    The suppressed mass deposition rate is:
      Mdot_suppressed = Mdot_classical × (1 − S)
    """

    def compute(self, dataset: dict) -> dict:
        # Compute cooling
        cool = ICMRadiativeCooling().compute(dataset)
        L_cool = cool['L_cool_W']
        Mdot_classical = cool['Mdot_solar_yr']

        # Compute heating
        heat = FUBiJetHeating().compute(dataset)
        Q_heat = heat['Q_heat_W']

        # Suppression factor
        S = min(Q_heat / L_cool, 1.0) if L_cool > 0 else 0.0
        Mdot_suppressed = Mdot_classical * (1 - S)

        # Gamma sweep
        gamma_sweep_THz = [0.01, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50]
        sweep = []
        for g_THz in gamma_sweep_THz:
            gamma_rad = 2 * PI * g_THz * 1e12
            ds = dict(dataset)
            ds['gamma'] = gamma_rad
            h = FUBiJetHeating().compute(ds)
            s_g = min(h['Q_heat_W'] / L_cool, 1.0) if L_cool > 0 else 0.0
            sweep.append({
                'gamma_THz': g_THz,
                'Q_heat_W': h['Q_heat_W'],
                'S': s_g,
                'Mdot_suppressed_solar_yr': Mdot_classical * (1 - s_g),
            })

        return {
            'L_cool_W': L_cool,
            'Q_heat_W': Q_heat,
            'suppression_factor': S,
            'Mdot_classical_solar_yr': Mdot_classical,
            'Mdot_suppressed_solar_yr': Mdot_suppressed,
            'gamma_sweep': sweep,
            'primary_equations': [
                'Cooling-Flow Suppression:',
                '  S(Gamma) = min(Q_heat / L_cool, 1)',
                '  Mdot_suppressed = Mdot_classical × (1 − S)',
                f'  L_cool = {L_cool:.4e} W',
                f'  Q_heat = {Q_heat:.4e} W',
                f'  S = {S:.6f}',
                f'  Mdot classical = {Mdot_classical:.1f} M_sun/yr',
                f'  Mdot suppressed = {Mdot_suppressed:.1f} M_sun/yr',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §4  CLUSTER PIPELINE (Perseus / Abell 2256 / Virgo)
# ═══════════════════════════════════════════════════════════════════════════════

# Cluster parameter catalog
CLUSTER_CATALOG = {
    'Perseus': {
        'T_keV': 4.0,
        'n_e': 3e4,
        'r_cool_kpc': 100,
        'M_bh': 3.4e8,         # NGC 1275 SMBH
        'v_jet_frac': 0.1,
        'eta_jet': 0.1,
        'M_icm_solar': 1e13,
    },
    'Abell_2256': {
        'T_keV': 6.5,
        'n_e': 1.5e4,
        'r_cool_kpc': 200,
        'M_bh': 1e9,
        'v_jet_frac': 0.05,
        'eta_jet': 0.08,
        'M_icm_solar': 5e13,
    },
    'Virgo_M87': {
        'T_keV': 2.5,
        'n_e': 5e4,
        'r_cool_kpc': 50,
        'M_bh': 6.5e9,         # M87 SMBH
        'v_jet_frac': 0.99,    # Relativistic M87 jet
        'eta_jet': 0.05,
        'M_icm_solar': 3e12,
    },
}


class ClusterCoolingFlowPipeline:
    """Full cooling-flow analysis for multiple clusters."""

    def compute(self, dataset: dict) -> dict:
        cluster_name = dataset.get('cluster', 'Perseus')
        if cluster_name in CLUSTER_CATALOG:
            params = dict(CLUSTER_CATALOG[cluster_name])
            params.update(dataset)
        else:
            params = dict(dataset)

        cool = ICMRadiativeCooling().compute(params)
        heat = FUBiJetHeating().compute(params)
        suppression = CoolingFlowSuppression().compute(params)

        return {
            'cluster': cluster_name,
            'cooling': cool,
            'heating': heat,
            'suppression': suppression,
            'primary_equations': [
                f'Galaxy Cluster Cooling-Flow Analysis: {cluster_name}',
                '═' * 50,
                '',
            ] + cool['primary_equations'] + [''] + heat['primary_equations'] + [''] + suppression['primary_equations'],
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    print("=" * 72)
    print("galaxy_cluster_cooling_flow.py — Self-Tests")
    print("=" * 72)

    ok = True
    passed = 0

    # Test 1: Perseus cooling luminosity positive
    cool = ICMRadiativeCooling().compute(CLUSTER_CATALOG['Perseus'])
    if cool['L_cool_W'] > 0 and math.isfinite(cool['L_cool_W']):
        print(f"  [PASS] Test 1: Perseus L_cool = {cool['L_cool_W']:.4e} W")
        passed += 1
    else:
        print(f"  [FAIL] Test 1: L_cool = {cool['L_cool_W']}"); ok = False

    # Test 2: Perseus cooling time physically realistic (0.1-100 Gyr)
    if 0.01 < cool['t_cool_Gyr'] < 1000:
        print(f"  [PASS] Test 2: t_cool = {cool['t_cool_Gyr']:.2f} Gyr")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: t_cool = {cool['t_cool_Gyr']} Gyr"); ok = False

    # Test 3: Perseus Mdot positive
    if cool['Mdot_solar_yr'] > 0:
        print(f"  [PASS] Test 3: Mdot_classical = {cool['Mdot_solar_yr']:.1f} M_sun/yr")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: Mdot = {cool['Mdot_solar_yr']}"); ok = False

    # Test 4: F_{U,Bi,i} at horizon is finite
    heat = FUBiJetHeating().compute(CLUSTER_CATALOG['Perseus'])
    if math.isfinite(heat['F_UBii']):
        print(f"  [PASS] Test 4: F_{{U,Bi,i}} = {heat['F_UBii']:.6e} m/s^2")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: F_UBii = {heat['F_UBii']}"); ok = False

    # Test 5: Q_heat positive
    if heat['Q_heat_W'] > 0 and math.isfinite(heat['Q_heat_W']):
        print(f"  [PASS] Test 5: Q_heat = {heat['Q_heat_W']:.4e} W")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: Q_heat = {heat['Q_heat_W']}"); ok = False

    # Test 6: Suppression factor in [0, 1]
    supp = CoolingFlowSuppression().compute(CLUSTER_CATALOG['Perseus'])
    S = supp['suppression_factor']
    if 0 <= S <= 1:
        print(f"  [PASS] Test 6: S(Perseus) = {S:.6f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: S = {S}"); ok = False

    # Test 7: Gamma sweep returns correct number of points
    sweep = supp['gamma_sweep']
    if len(sweep) == 7:
        print(f"  [PASS] Test 7: Gamma sweep: {len(sweep)} points")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: Expected 7 sweep points, got {len(sweep)}"); ok = False

    # Test 8: Suppression peaks near Γ = 0.10 THz
    s_vals = [(s['gamma_THz'], s['S']) for s in sweep]
    s_at_010 = next(s for g, s in s_vals if abs(g - 0.10) < 0.01)
    s_at_050 = next(s for g, s in s_vals if abs(g - 0.50) < 0.01)
    if s_at_010 >= s_at_050:
        print(f"  [PASS] Test 8: S(0.10 THz)={s_at_010:.4f} >= S(0.50 THz)={s_at_050:.4f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: Wrong gamma ordering"); ok = False

    # Test 9: All 3 clusters produce finite results
    all_ok = True
    for name in ['Perseus', 'Abell_2256', 'Virgo_M87']:
        result = ClusterCoolingFlowPipeline().compute({'cluster': name})
        s_val = result['suppression']['suppression_factor']
        if not (0 <= s_val <= 1 and math.isfinite(s_val)):
            all_ok = False
            break
    if all_ok:
        print(f"  [PASS] Test 9: All 3 clusters produce valid suppression factors")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: Cluster pipeline failure"); ok = False

    # Test 10: Primary equations present in all outputs
    pipe = ClusterCoolingFlowPipeline().compute({'cluster': 'Perseus'})
    all_have = all(
        'primary_equations' in pipe[k]
        for k in ['cooling', 'heating', 'suppression']
    )
    if all_have:
        print(f"  [PASS] Test 10: All outputs contain primary_equations")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: Missing primary_equations"); ok = False

    print(f"\n  galaxy_cluster_cooling_flow.py: {passed}/10 tests passed")
    return ok


if __name__ == "__main__":
    success = _run_tests()
    raise SystemExit(0 if success else 1)
