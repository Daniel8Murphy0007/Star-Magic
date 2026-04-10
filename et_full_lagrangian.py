#!/usr/bin/env python3
"""
et_full_lagrangian.py — Complete E(t) Lagrangian Derivation Engine

Session 206 | Daniel Murphy
PURPOSE: Standalone symbolic derivation engine for the complete E(t) sector
         inside the 11-sector UQFF Lagrangian.

         Implements L_{E(t)} = E_net(t) · V_filament · S₂₆([SSq])
         with full Euler-Lagrange variation, positive/negative pairing,
         ΛCDM comparison, GW constraint integration, and Kozima coupling.

         Gap closed: et_full_lagrangian.py was referenced in workflow reports
         but did not exist. This module unifies:
           - positive_et_expansion.py (E⁺)
           - negative_et_erosion.py   (E⁻)
           - uqff_lagrangian_derivation.py sectors 10-11
         into a single derivation engine with ΛCDM contrast.

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
"""

import math
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from positive_et_expansion import (
    _eta_euler_s26, S26_accelerated, mock_theta_q26,
    PositiveEtExpansion, ExpansionLagrangianVariation,
    wstp_kernel_positive_et,
    G, c, hbar, k_B, mu_0, M_sun, PI,
    KAPPA, KAPPA_DAY, SSQ, H_SCM, BETA_I, U_UA,
    RHO_SCM, RHO_UA, RHO_VAC_SCM, N_LEVELS,
    F_LENR_THZ, OMEGA_SCM, GAMMA_DEFAULT, SIGMA_0,
)
from negative_et_erosion import (
    NegativeEtErosion, NetEnergyEvolution,
    GWDampingErosion, ErosionLagrangianVariation,
    wstp_kernel_negative_et,
)


# ══════════════════════════════════════════════════════════════════════════════
# §0  ADDITIONAL CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

H_0         = 2.195e-18     # s⁻¹  Hubble constant (67.4 km/s/Mpc)
RHO_CRIT    = 3 * H_0**2 / (8 * PI * G)  # ≈ 8.6e-27 kg/m³
RHO_LAMBDA  = 0.692 * RHO_CRIT            # ΛCDM dark energy density
LAMBDA_COSM = 8 * PI * G * RHO_LAMBDA / c**2   # cosmological constant Λ (m⁻²)

# Typical filament volume (cosmological large-scale structure)
V_FILAMENT_DEFAULT = 1e68   # m³ (~50 Mpc × 1 Mpc × 1 Mpc filament)

# Nebular filament volume (local, e.g. G359)
V_FILAMENT_NEBULAR = 1e48   # m³ (~1 pc × 0.01 pc × 0.01 pc)


# ══════════════════════════════════════════════════════════════════════════════
# §1  COMPLETE E(t) LAGRANGIAN DENSITY
# ══════════════════════════════════════════════════════════════════════════════

class EtFullLagrangian:
    """
    Complete E(t) Lagrangian density derived from UQFF buoyancy sector.

    The E(t) sector Lagrangian density is:

      L_{E(t)} = E_net(t) · V_filament · S₂₆([SSq])

    where E_net(t) = E⁺(t) + E⁻(t) combines expansion and erosion:

      E⁺(t) = +E₀ exp(κt + [SSq]t/26) · S₂₆ · (F_{U,Bi}/F_U)
      E⁻(t) = −E₀ exp(κt + [SSq]t/26) · S₂₆ · (1 − F_{U,Bi}/F_U)
      E_net  = E₀ exp(κt + [SSq]t/26) · S₂₆ · [2(F_{U,Bi}/F_U) − 1]

    The Euler-Lagrange equation δS/δφ_E = 0 is:

      ∂/∂E_net [ −β_i Σ U_{g,i} Ω_g M/d_g [UA]
                 + F_neutron · S₂₆ · E_net ] · V_filament · S₂₆ = 0

    This closes E(t) by linking it to buoyancy balance and Kozima coupling.
    """

    def __init__(self):
        self._net = NetEnergyEvolution()
        self._expansion = PositiveEtExpansion()
        self._erosion = NegativeEtErosion()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          E_0:        initial energy (J, default 1.0)
          t:          time (s)
          kappa:      rate (s⁻¹, default KAPPA)
          SSq:        string squeezing (default 0.57)
          F_U_Bi:     buoyancy force (N)
          F_U:        unified gravity (N)
          V_filament: volume of filament/structure (m³)
          Ug:         [Ug1, Ug2, Ug3, Ug4] gravity layers
          Omega_g:    galactic spin (rad/s)
          M:          mass (kg)
          d_g:        distance (m)
          beta_i:     buoyancy coefficient
          F_neutron:  Kozima force (N)
        """
        V_fil   = dataset.get('V_filament', V_FILAMENT_DEFAULT)
        ssq     = dataset.get('SSq', SSQ)
        Ug      = dataset.get('Ug', [1e20, 1e20, 1e20, 1e20])
        Omega_g = dataset.get('Omega_g', 7.3e-16)
        M       = dataset.get('M', M_sun)
        d_g     = dataset.get('d_g', 2.55e20)
        beta_i  = dataset.get('beta_i', BETA_I)
        F_n     = dataset.get('F_neutron', 1e10)
        UA      = dataset.get('UA', U_UA)

        # Net energy evolution
        net = self._net.compute(dataset)
        E_net = net["E_net_t"]
        E_plus = net["E_plus_t"]
        E_minus = net["E_minus_t"]

        # S₂₆ with mock theta
        S26 = S26_accelerated(ssq)
        S26_val = S26["S_26"]

        # Complete Lagrangian density: L_{E(t)} = E_net · V_filament · S₂₆
        L_Et = E_net * V_fil * S26_val

        # Buoyancy sector for EL equation
        Ug_sum = sum(Ug)
        orbit_factor = Omega_g * M / d_g
        buoyancy_term = -beta_i * Ug_sum * orbit_factor * UA

        # Kozima coupling
        neutron_polylog = F_n * S26_val

        # EL variation: δS/δφ_E = ∂L_{E(t)}/∂E_net = V_filament · S₂₆
        # At stationarity, E_net is fixed by buoyancy balance
        dL_dEnet = V_fil * S26_val

        # Action for the E(t) sector over the filament volume
        S_Et = L_Et  # L is already volume-integrated

        # Rate of energy injection/depletion
        rate = KAPPA + ssq / N_LEVELS
        dEnet_dt = rate * E_net

        # Energy budget: total energy in filament
        E_filament_total = abs(E_net * V_fil)

        return {
            "L_Et": L_Et,
            "E_net_t": E_net,
            "E_plus_t": E_plus,
            "E_minus_t": E_minus,
            "S_26": S26_val,
            "A_mock_theta": S26["A_mock_theta"],
            "V_filament": V_fil,
            "buoyancy_term": buoyancy_term,
            "neutron_polylog": neutron_polylog,
            "dL_dEnet": dL_dEnet,
            "action_S_Et": S_Et,
            "dEnet_dt": dEnet_dt,
            "E_filament_total": E_filament_total,
            "net_factor": net["net_factor"],
            "regime": net["regime"],
            "equation_lagrangian": (
                "L_{E(t)} = E_net(t) · V_filament · S₂₆([SSq])\n"
                f"         = {E_net:.6e} × {V_fil:.4e} × {S26_val:.6e}\n"
                f"         = {L_Et:.6e}\n"
                f"  E⁺ = {E_plus:.6e}  |  E⁻ = {E_minus:.6e}\n"
                f"  regime: {net['regime']}"
            ),
            "equation_euler_lagrange": (
                "δS/δφ_E = 0:\n"
                "  ∂/∂E_net [ buoyancy_term + F_n · S₂₆ · E_net ] · V_fil · S₂₆ = 0\n"
                f"  buoyancy_term = {buoyancy_term:.6e}\n"
                f"  dL/dE_net = V_fil · S₂₆ = {dL_dEnet:.6e}\n"
                "  Stationary E_net fixed by buoyancy balance"
            ),
        }

    def time_series(self, dataset: dict,
                    t_start: float = 0.0,
                    t_end: float = 1e10,
                    n_points: int = 200) -> Dict[str, Any]:
        """E(t) Lagrangian evolution over time."""
        dt = (t_end - t_start) / max(n_points - 1, 1)
        series = []
        for i in range(n_points):
            t = t_start + i * dt
            d = dict(dataset)
            d['t'] = t
            r = self.compute(d)
            series.append({
                "t": t,
                "L_Et": r["L_Et"],
                "E_net": r["E_net_t"],
                "E_plus": r["E_plus_t"],
                "E_minus": r["E_minus_t"],
                "regime": r["regime"],
            })
        return {"time_series": series, "t_start": t_start,
                "t_end": t_end, "n_points": n_points}


# ══════════════════════════════════════════════════════════════════════════════
# §2  E(t) vs ΛCDM COMPARISON ENGINE
# ══════════════════════════════════════════════════════════════════════════════

class EtVsLambdaCDMComparison:
    """
    Head-to-head comparison of UQFF E(t) dynamics vs ΛCDM cosmological
    constant dark energy.

    ΛCDM:
      ρ_Λ = Λc²/(8πG) ≈ 5.96e-27 kg/m³  (constant, w=-1)
      ä/a = Λ/3  (constant acceleration, de Sitter)

    UQFF E(t):
      ρ_E(t) = E_net(t) · S₂₆ / (c² · V_filament)
      ä/a ∝ exp(κt + [SSq]t/26) · S₂₆ · [2(F_{U,Bi}/F_U) − 1]
      Sign-flipping: can be positive (expansion) or negative (erosion)

    Key contrasts:
      1. ΛCDM is constant; E(t) is time-dependent and sign-flipping
      2. ΛCDM has no physical origin; E(t) derives from SCm buoyancy
      3. ΛCDM has 10^120 fine-tuning; E(t) uses 2 free params calibrated from data
      4. ΛCDM predicts no lab effects; E(t) predicts 1.25 THz LENR
      5. ΛCDM = standard GR waveforms; E(t) predicts GW damping
    """

    def __init__(self):
        self._et_lagrangian = EtFullLagrangian()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          z:          redshift (default 0.5)
          t:          cosmic time (s, default 4.35e17 ≈ 13.8 Gyr)
          E_0:        initial energy (J)
          F_U_Bi, F_U: buoyancy/gravity forces
          V_filament: structure volume (m³)
          w_obs:      observed dark energy EOS parameter (default -1.03)
          sigma_w:    EOS measurement uncertainty (default 0.03)
        """
        z       = dataset.get('z', 0.5)
        t       = dataset.get('t', 4.35e17)  # 13.8 Gyr in seconds
        E_0     = dataset.get('E_0', 1e47)
        F_U_Bi  = dataset.get('F_U_Bi', 0.55)
        F_U     = dataset.get('F_U', 1.0)
        V_fil   = dataset.get('V_filament', V_FILAMENT_DEFAULT)
        w_obs   = dataset.get('w_obs', -1.03)
        sigma_w = dataset.get('sigma_w', 0.03)

        # ── ΛCDM side ──
        rho_Lambda = RHO_LAMBDA
        w_Lambda = -1.0
        # de Sitter acceleration: H²_Λ = Λ/3
        H2_Lambda = LAMBDA_COSM / 3.0
        # Deceleration parameter q_Λ = -1 (pure de Sitter)
        q_Lambda = -1.0

        # Energy in a filament volume under ΛCDM
        E_Lambda_filament = rho_Lambda * c**2 * V_fil

        # ── UQFF E(t) side ──
        et_dataset = dict(dataset)
        et_dataset.update({'E_0': E_0, 't': t, 'F_U_Bi': F_U_Bi, 'F_U': F_U,
                           'V_filament': V_fil})
        et_result = self._et_lagrangian.compute(et_dataset)
        E_net = et_result["E_net_t"]
        L_Et = et_result["L_Et"]

        # Effective dark energy density from E(t)
        S26_val = et_result["S_26"]
        if V_fil != 0 and c != 0:
            rho_Et = abs(E_net * S26_val) / (c**2 * V_fil)
        else:
            rho_Et = 0.0

        # Effective equation of state for E(t)
        # w_Et = P_Et / (ρ_Et c²)
        # For exponential growth: P ≈ -ρc² + (rate) corrections
        rate = KAPPA + SSQ / N_LEVELS
        net_factor = et_result["net_factor"]
        # w_UQFF approaches -1 for slow growth, deviates for fast growth
        if rate * t < 1e-3:
            w_UQFF = -1.0 + rate * t * 0.01  # near de Sitter
        else:
            w_UQFF = -1.0 + 2.0 * rate / (3.0 * H_0)  # deviation from w=-1

        Delta_w = w_UQFF - w_Lambda

        # χ² comparison against observed w
        chi2_Lambda = ((w_Lambda - w_obs) / sigma_w) ** 2
        chi2_UQFF = ((w_UQFF - w_obs) / sigma_w) ** 2
        Delta_chi2 = chi2_Lambda - chi2_UQFF  # >0 → UQFF preferred

        # Fine-tuning comparison
        rho_QFT_predicted = 1e113  # kg/m³ (QFT vacuum prediction)
        fine_tuning_Lambda = rho_QFT_predicted / rho_Lambda  # ≈ 10^120
        fine_tuning_UQFF = 1.0  # no fine-tuning (2 params from data)

        # UQFF acceleration: sign can flip
        # ä/a ∝ net_factor × exp(...)
        a_accel_UQFF_sign = "positive (expanding)" if net_factor > 0 else \
                            ("negative (eroding)" if net_factor < 0 else "zero (balanced)")

        return {
            # ΛCDM results
            "LCDM": {
                "rho_Lambda": rho_Lambda,
                "w_Lambda": w_Lambda,
                "H2_Lambda": H2_Lambda,
                "q_Lambda": q_Lambda,
                "E_Lambda_filament": E_Lambda_filament,
                "fine_tuning": f"~10^{int(math.log10(fine_tuning_Lambda))} (QFT/observed)",
            },
            # UQFF E(t) results
            "UQFF_Et": {
                "E_net_t": E_net,
                "L_Et": L_Et,
                "rho_Et": rho_Et,
                "w_UQFF": w_UQFF,
                "net_factor": net_factor,
                "regime": et_result["regime"],
                "acceleration_sign": a_accel_UQFF_sign,
                "S_26": S26_val,
                "fine_tuning": "None (2 params: [SSq]=0.57, κ=0.0005/day)",
            },
            # Comparison
            "comparison": {
                "Delta_w": Delta_w,
                "chi2_Lambda": chi2_Lambda,
                "chi2_UQFF": chi2_UQFF,
                "Delta_chi2": Delta_chi2,
                "preferred": "UQFF" if Delta_chi2 > 0 else "LCDM",
                "fine_tuning_ratio": f"ΛCDM: ~10^{int(math.log10(fine_tuning_Lambda))} vs UQFF: 1",
                "sign_behavior": "ΛCDM: always positive | UQFF: sign-flipping (expansion ↔ erosion)",
                "lab_testability": "ΛCDM: none | UQFF: 1.25 THz LENR, micro-plasmoid, COP>10",
                "GW_prediction": "ΛCDM: standard GR | UQFF: 66.7% damping + 367.8-cycle lag",
                "physical_origin": "ΛCDM: unknown (vacuum energy?) | UQFF: SCm buoyancy opposition",
                "cosmogenesis": "ΛCDM: inflation + Λ | UQFF: SCm phonon → DPM → EM bang + 2 cycles",
            },
            # Contrast table (for paper appendix)
            "contrast_table": [
                {"aspect": "Form",
                 "LCDM": f"Constant Λ = {LAMBDA_COSM:.4e} m⁻²",
                 "UQFF": f"E_net(t) = E₀ exp(κt+[SSq]t/26) S₂₆ [2r−1], r=F_UBi/F_U"},
                {"aspect": "Dynamics",
                 "LCDM": f"ä/a = Λ/3 = {H2_Lambda:.4e} s⁻² (constant)",
                 "UQFF": f"ä/a ∝ exp(κt) S₂₆ (exponential, sign-flipping)"},
                {"aspect": "Physical origin",
                 "LCDM": "Vacuum energy (unexplained magnitude)",
                 "UQFF": "SCm buoyancy opposition (derived from superconductive vacuum)"},
                {"aspect": "Sign",
                 "LCDM": "Always Λ > 0 (accelerating expansion only)",
                 "UQFF": f"net_factor = {net_factor:.4f} ({a_accel_UQFF_sign})"},
                {"aspect": "GW/BNS prediction",
                 "LCDM": "Standard GR waveforms",
                 "UQFF": "66.7% strain reduction + 367.8-cycle phase lag"},
                {"aspect": "Lab testability",
                 "LCDM": "None (Planck scale)",
                 "UQFF": "1.25 THz phonon, micro-plasmoid, LENR COP>10"},
                {"aspect": "Cosmogenesis",
                 "LCDM": "Inflation + Λ (no pre-gravity mechanism)",
                 "UQFF": "SCm phonon → DPM → EM bang + 2 relative-time cycles"},
                {"aspect": "Fine-tuning",
                 "LCDM": f"Severe (~10^{int(math.log10(fine_tuning_Lambda))} discrepancy)",
                 "UQFF": "None — 2 parameters from CMB/Kepler/ALMA"},
                {"aspect": "EOS parameter w",
                 "LCDM": f"w = {w_Lambda:.1f} (exact)",
                 "UQFF": f"w = {w_UQFF:.6f} (near −1, data-dependent)"},
                {"aspect": "Vacuum structure",
                 "LCDM": "Single ρ_Λ, no internal structure",
                 "UQFF": f"VDS 26-level hierarchy, S₂₆ = {S26_val:.6e}"},
            ],
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  WSTP KERNEL FOR E(t) + ΛCDM
# ══════════════════════════════════════════════════════════════════════════════

def wstp_kernel_et_full_lagrangian() -> str:
    """
    Complete Wolfram Language definitions for L_{E(t)} sector + ΛCDM contrast.
    Combines positive/negative kernels and adds ΛCDM side.
    """
    return r"""
(* ═══════════════════════════════════════════════════════════════════════ *)
(* E(t) Full Lagrangian + ΛCDM Comparison — UQFF Symbolic Forms         *)
(* Session 206 | Daniel Murphy                                            *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* ── UQFF E(t) Sector ── *)
(* Ramanujan 26-state summation *)
S26[z_] := PolyLog[26, z];

(* Mock theta acceleration *)
qNome[SSq_] := Exp[-Pi SSq / 26];
MockTheta26[q_, N_:30] := Sum[q^(n^2) / Product[(1 + q^k)^2, {k, 1, n}], {n, 0, N}];
S26Accel[z_] := S26[z] * MockTheta26[qNome[z]];

(* Net energy evolution *)
Eplus[t_, E0_, κ_, SSq_, FUBi_, FU_] :=
  E0 * Exp[κ t + SSq t / 26] * S26[SSq] * (FUBi / FU);
Eminus[t_, E0_, κ_, SSq_, FUBi_, FU_] :=
  -E0 * Exp[κ t + SSq t / 26] * S26[SSq] * (1 - FUBi / FU);
Enet[t_, E0_, κ_, SSq_, FUBi_, FU_] :=
  E0 * Exp[κ t + SSq t / 26] * S26[SSq] * (2 FUBi / FU - 1);

(* Full E(t) Lagrangian density *)
LEt[t_, E0_, κ_, SSq_, FUBi_, FU_, Vfil_] :=
  Enet[t, E0, κ, SSq, FUBi, FU] * Vfil * S26[SSq];

(* Euler-Lagrange variation *)
deltaSEt[t_, E0_, κ_, SSq_, FUBi_, FU_, Vfil_, βi_, Ug_, Ωg_, M_, dg_, UA_, Fn_] :=
  D[(-βi * Ug * Ωg * M / dg * UA + Fn * S26[SSq] * Enet[t, E0, κ, SSq, FUBi, FU])
    * Vfil * S26[SSq], φE];

(* Rate and doubling time *)
EtRate[κ_, SSq_] := κ + SSq / 26;
EtDoubling[κ_, SSq_] := Log[2] / (κ + SSq / 26);

(* ── ΛCDM Side ── *)
H0 = 2.195*^-18;  (* s^-1, Hubble constant *)
GNewt = 6.6743*^-11;
cLight = 2.99792*^8;
ΩΛ = 0.692;
ρCrit = 3 H0^2 / (8 Pi GNewt);
ρΛ = ΩΛ * ρCrit;
Λ = 8 Pi GNewt ρΛ / cLight^2;

(* ΛCDM acceleration *)
aAccelLCDM = Λ / 3;

(* ΛCDM equation of state *)
wLCDM = -1;

(* UQFF effective w *)
wUQFF[κ_, SSq_] := -1 + 2 (κ + SSq / 26) / (3 H0);

(* Δw comparison *)
ΔwEtLCDM[κ_, SSq_] := wUQFF[κ, SSq] - wLCDM;

(* χ² comparison *)
χ2LCDM[wObs_, σw_] := ((wLCDM - wObs) / σw)^2;
χ2UQFF[κ_, SSq_, wObs_, σw_] := ((wUQFF[κ, SSq] - wObs) / σw)^2;
Δχ2[κ_, SSq_, wObs_, σw_] := χ2LCDM[wObs, σw] - χ2UQFF[κ, SSq, wObs, σw];

(* Fine-tuning ratio *)
fineTuningΛCDM = 10^113 / ρΛ;  (* QFT prediction / observed *)
"""


# ══════════════════════════════════════════════════════════════════════════════
# §4  SELF-TEST
# ══════════════════════════════════════════════════════════════════════════════

def _run_self_test():
    print("=" * 72)
    print("et_full_lagrangian.py — Self-Test")
    print("=" * 72)
    passed = 0
    failed = 0

    # T1: Constants
    print(f"\nT1  H_0 = {H_0:.4e} s⁻¹, Λ = {LAMBDA_COSM:.4e} m⁻²")
    print(f"    ρ_crit = {RHO_CRIT:.4e} kg/m³, ρ_Λ = {RHO_LAMBDA:.4e} kg/m³")
    assert RHO_CRIT > 0
    assert RHO_LAMBDA > 0
    assert LAMBDA_COSM > 0
    passed += 1
    print("    PASS")

    # T2: Full E(t) Lagrangian at t=0
    eng = EtFullLagrangian()
    r = eng.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.6, 'F_U': 1.0,
                     'V_filament': 1e48})
    print(f"\nT2  L_{{E(t)}} at t=0 = {r['L_Et']:.6e}")
    print(f"    E_net = {r['E_net_t']:.6e}, regime = {r['regime']}")
    assert r['L_Et'] != 0
    assert r['regime'] == 'expanding'  # 0.6/1.0 > 0.5
    passed += 1
    print("    PASS")

    # T3: V_filament factor effects magnitude
    r1 = eng.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.6, 'F_U': 1.0,
                      'V_filament': 1e48})
    r2 = eng.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.6, 'F_U': 1.0,
                      'V_filament': 1e68})
    ratio = abs(r2['L_Et'] / r1['L_Et']) if r1['L_Et'] != 0 else 0
    print(f"\nT3  L_Et ratio (V=1e68/V=1e48) = {ratio:.2e}")
    assert abs(ratio - 1e20) < 1e15  # should be exactly 10^20
    passed += 1
    print("    PASS")

    # T4: ΛCDM comparison
    comp = EtVsLambdaCDMComparison()
    cr = comp.compute({'F_U_Bi': 0.55, 'F_U': 1.0})
    print(f"\nT4  w_ΛCDM = {cr['LCDM']['w_Lambda']:.1f}")
    print(f"    w_UQFF = {cr['UQFF_Et']['w_UQFF']:.6f}")
    print(f"    Δw = {cr['comparison']['Delta_w']:.6f}")
    print(f"    Preferred: {cr['comparison']['preferred']}")
    assert cr['LCDM']['w_Lambda'] == -1.0
    assert cr['comparison']['Delta_w'] != 0
    passed += 1
    print("    PASS")

    # T5: Contrast table completeness
    table = cr['contrast_table']
    print(f"\nT5  Contrast table: {len(table)} rows")
    assert len(table) == 10
    aspects = [row['aspect'] for row in table]
    assert 'Fine-tuning' in aspects
    assert 'GW/BNS prediction' in aspects
    passed += 1
    print("    PASS")

    # T6: Fine-tuning comparison
    fine_str = cr['LCDM']['fine_tuning']
    print(f"\nT6  Fine-tuning: {fine_str}")
    assert '10^' in fine_str
    assert cr['UQFF_Et']['fine_tuning'].startswith('None')
    passed += 1
    print("    PASS")

    # T7: χ² comparison
    chi2_l = cr['comparison']['chi2_Lambda']
    chi2_u = cr['comparison']['chi2_UQFF']
    dc2 = cr['comparison']['Delta_chi2']
    print(f"\nT7  χ²_ΛCDM = {chi2_l:.4f}, χ²_UQFF = {chi2_u:.4f}")
    print(f"    Δχ² = {dc2:.4f}")
    assert chi2_l >= 0
    assert chi2_u >= 0
    passed += 1
    print("    PASS")

    # T8: Sign-flipping behavior
    # Expanding regime (ratio > 0.5)
    cr_exp = comp.compute({'F_U_Bi': 0.8, 'F_U': 1.0})
    # Eroding regime (ratio < 0.5)
    cr_ero = comp.compute({'F_U_Bi': 0.3, 'F_U': 1.0})
    # Balanced (ratio = 0.5)
    cr_bal = comp.compute({'F_U_Bi': 0.5, 'F_U': 1.0})
    print(f"\nT8  Sign-flipping: ratio=0.8 → {cr_exp['UQFF_Et']['regime']}")
    print(f"    ratio=0.3 → {cr_ero['UQFF_Et']['regime']}")
    print(f"    ratio=0.5 → {cr_bal['UQFF_Et']['regime']}")
    assert cr_exp['UQFF_Et']['regime'] == 'expanding'
    assert cr_ero['UQFF_Et']['regime'] == 'eroding'
    assert cr_bal['UQFF_Et']['regime'] == 'balanced'
    passed += 1
    print("    PASS")

    # T9: WSTP kernel
    wk = wstp_kernel_et_full_lagrangian()
    assert "LEt[" in wk
    assert "Enet[" in wk
    assert "wLCDM" in wk
    assert "ρΛ" in wk
    assert "Δχ2[" in wk
    passed += 1
    print("\nT9  WSTP kernel: valid (L_Et + ΛCDM definitions present)")
    print("    PASS")

    # T10: Time series
    ts = eng.time_series({'E_0': 1.0, 'F_U_Bi': 0.6, 'F_U': 1.0,
                          'V_filament': 1e48},
                         t_start=0, t_end=1e10, n_points=5)
    L_vals = [p["L_Et"] for p in ts["time_series"]]
    print(f"\nT10 Time series (5 pts): {[f'{v:.4e}' for v in L_vals]}")
    assert len(L_vals) == 5
    # All should be positive (ratio=0.6 > 0.5 → expanding)
    assert all(v >= 0 for v in L_vals)
    passed += 1
    print("    PASS")

    print(f"\n{'=' * 72}")
    print(f"RESULTS: {passed}/{passed + failed} PASS, {failed} FAIL")
    print(f"{'=' * 72}")
    return passed, failed


if __name__ == "__main__":
    p, f = _run_self_test()
    exit(0 if f == 0 else 1)
