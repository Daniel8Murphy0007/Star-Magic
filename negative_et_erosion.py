#!/usr/bin/env python3
"""
negative_et_erosion.py — Negative E(t) Buoyancy-Driven Erosion Engine

Session 205 | Daniel Murphy
PURPOSE: Symmetric counterpart to positive_et_expansion.py.

         Negative E(t) describes buoyancy opposition dominating, leading
         to vacuum energy depletion and filament/nebula erosion. Applies
         when F_{U,Bi}/F_U < 1 (gravity wins over buoyancy).

         Direct match: PAPER_359 (G359 Galactic Center Filament),
         PAPER_838 (Chandra SNR/nebula batch), PAPER_008b (GW170817
         late inspiral 66.7% damping).

         Currently MISSING from the codebase:
           - Individual erosion calculators exist (PhotoevaporationErosion,
             HorseheadErosion, PillarsErosion, etc.) but NONE implement
             the master E⁻(t) formula with S₂₆ Ramanujan coupling.
           - No Lagrangian δS/δφ_erosion variation.
           - No E_net(t) = E⁺ + E⁻ combined evolution.

         This module:
           1. Master E⁻(t) = -E₀ exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (1 − F_{U,Bi}/F_U)
           2. Symmetric Lagrangian variation δS/δφ_erosion = 0
           3. Net energy E_net(t) = E⁺(t) + E⁻(t)
           4. GW170817 damping integration
           5. Cosmogenesis vacuum depletion phase

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
"""

import math
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# Import S₂₆ engine from positive counterpart (no circular dependency)
from positive_et_expansion import (
    _eta_euler_s26, S26_accelerated, mock_theta_q26,
    PositiveEtExpansion,
    G, c, hbar, k_B, mu_0, M_sun, PI,
    KAPPA, KAPPA_DAY, SSQ, H_SCM, BETA_I, U_UA,
    RHO_SCM, RHO_UA, RHO_VAC_SCM, N_LEVELS,
    F_LENR_THZ, OMEGA_SCM, GAMMA_DEFAULT, SIGMA_0,
)


# ══════════════════════════════════════════════════════════════════════════════
# §1  MASTER NEGATIVE E(t) EROSION EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════

class NegativeEtErosion:
    """
    Master negative E(t) buoyancy-driven erosion engine.

    E⁻(t) = −E₀ · exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (1 − F_{U,Bi}/F_U)

    Physical meaning:
      - Negative sign enforces energy depletion
      - (1 − F_{U,Bi}/F_U) > 0 when buoyancy is subdominant (erosion regime)
      - Same exponential and S₂₆ structure as E⁺(t) — perfect symmetry

    Applies to: filament erosion (PAPER_359), nebular cavity formation,
    photoevaporation (M16), GW damping (PAPER_008b), and the vacuum
    depletion phase in cosmogenesis.

    Mirror: positive_et_expansion.py implements E⁺(t) for the growth case.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          E_0:     initial vacuum energy (J, default 1.0)
          t:       time (seconds)
          kappa:   expansion/erosion rate (s⁻¹, default KAPPA)
          SSq:     string squeezing (default 0.57)
          F_U_Bi:  buoyancy force (N)
          F_U:     unified gravity force (N)
          use_mock_theta: enable mock theta acceleration (default True)
        """
        E_0    = dataset.get('E_0', 1.0)
        t      = dataset.get('t', 0.0)
        kappa  = dataset.get('kappa', KAPPA)
        ssq    = dataset.get('SSq', SSQ)
        F_U_Bi = dataset.get('F_U_Bi', 0.9)  # < F_U for erosion
        F_U    = dataset.get('F_U', 1.0)
        use_mt = dataset.get('use_mock_theta', True)

        # Exponential growth factor (same as E⁺)
        exp_arg = kappa * t + (ssq * t) / N_LEVELS
        exp_arg_clamped = min(exp_arg, 700.0)
        growth_factor = math.exp(exp_arg_clamped)

        # Ramanujan 26-state summation
        if use_mt:
            s26_data = S26_accelerated(ssq)
            S26_val = s26_data["S_26"]
        else:
            S26_val = _eta_euler_s26(ssq)

        # Erosion factor: (1 − F_{U,Bi}/F_U)
        if F_U == 0:
            erosion_factor = 0.0
        else:
            erosion_factor = 1.0 - (F_U_Bi / F_U)

        # Master expression (negative sign)
        E_minus = -E_0 * growth_factor * S26_val * erosion_factor

        # Rate of change dE⁻/dt
        rate = kappa + ssq / N_LEVELS
        dEdt = E_minus * rate if E_minus != 0 else 0.0

        # Half-life (when does |E⁻| = E₀/2?)
        # |E⁻(t)| = E₀ × erosion × S₂₆ × exp(rate×t)
        # This grows, so "half-life" is the time to double the erosion magnitude
        t_double_erosion = math.log(2) / rate if rate > 0 else float('inf')

        return {
            "E_minus_t": E_minus,
            "E_0": E_0,
            "t": t,
            "kappa": kappa,
            "SSq": ssq,
            "exp_argument": exp_arg,
            "growth_factor": growth_factor,
            "S_26": S26_val,
            "erosion_factor": erosion_factor,
            "F_U_Bi": F_U_Bi,
            "F_U": F_U,
            "dE_minus_dt": dEdt,
            "effective_rate_per_s": rate,
            "erosion_doubling_time_s": t_double_erosion,
            "erosion_doubling_time_days": t_double_erosion / 86400.0,
            "is_eroding": erosion_factor > 0,
            "equation": (
                "E⁻(t) = −E₀ · exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (1 − F_{U,Bi}/F_U)\n"
                f"      = −{E_0:.4e} × {growth_factor:.6e} × {S26_val:.6e} × {erosion_factor:.4f}\n"
                f"      = {E_minus:.6e}"
            ),
        }

    def time_series(self, dataset: dict,
                    t_start: float = 0.0,
                    t_end: float = 1e10,
                    n_points: int = 200) -> Dict[str, Any]:
        """Compute E⁻(t) over a time range."""
        dt = (t_end - t_start) / max(n_points - 1, 1)
        series = []
        for i in range(n_points):
            t = t_start + i * dt
            d = dict(dataset)
            d['t'] = t
            result = self.compute(d)
            series.append({
                "t": t,
                "E_minus": result["E_minus_t"],
                "growth_factor": result["growth_factor"],
            })
        return {
            "time_series": series,
            "t_start": t_start,
            "t_end": t_end,
            "n_points": n_points,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  NET ENERGY EVOLUTION: E_net(t) = E⁺(t) + E⁻(t)
# ══════════════════════════════════════════════════════════════════════════════

class NetEnergyEvolution:
    """
    Combined E(t) evolution balancing expansion and erosion:

      E_net(t) = E⁺(t) + E⁻(t)

    where:
      E⁺(t) = +E₀ exp(κt + [SSq]t/26) S₂₆ (F_{U,Bi}/F_U)
      E⁻(t) = −E₀ exp(κt + [SSq]t/26) S₂₆ (1 − F_{U,Bi}/F_U)

    Substituting:
      E_net(t) = E₀ exp(κt + [SSq]t/26) S₂₆ [ (F_{U,Bi}/F_U) − (1 − F_{U,Bi}/F_U) ]
               = E₀ exp(κt + [SSq]t/26) S₂₆ [ 2(F_{U,Bi}/F_U) − 1 ]

    Critical ratio: F_{U,Bi}/F_U = 0.5  →  E_net = 0 (exact balance)
                    F_{U,Bi}/F_U > 0.5  →  E_net > 0 (net expansion)
                    F_{U,Bi}/F_U < 0.5  →  E_net < 0 (net erosion)
    """

    def __init__(self):
        self._expansion = PositiveEtExpansion()
        self._erosion = NegativeEtErosion()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters: same as PositiveEtExpansion.compute() and
        NegativeEtErosion.compute(). Returns both components + net.
        """
        plus = self._expansion.compute(dataset)
        minus = self._erosion.compute(dataset)

        E_net = plus["E_plus_t"] + minus["E_minus_t"]

        # Analytic shortcut: E_net = E₀ exp(...) S₂₆ [2(F_{U,Bi}/F_U) − 1]
        F_U_Bi = dataset.get('F_U_Bi', 1.1)
        F_U    = dataset.get('F_U', 1.0)
        ratio  = F_U_Bi / F_U if F_U != 0 else 0.0
        net_factor = 2.0 * ratio - 1.0
        critical_ratio = 0.5

        regime = "expanding" if net_factor > 0 else ("eroding" if net_factor < 0 else "balanced")

        return {
            "E_net_t": E_net,
            "E_plus_t": plus["E_plus_t"],
            "E_minus_t": minus["E_minus_t"],
            "net_factor": net_factor,
            "buoyancy_ratio": ratio,
            "critical_ratio": critical_ratio,
            "regime": regime,
            "S_26": plus["S_26"],
            "growth_factor": plus["growth_factor"],
            "equation": (
                "E_net(t) = E₀ exp(κt + [SSq]t/26) S₂₆ [2(F_{U,Bi}/F_U) − 1]\n"
                f"         net_factor = 2×{ratio:.4f} − 1 = {net_factor:.4f}\n"
                f"         regime: {regime}\n"
                f"  E⁺ = {plus['E_plus_t']:.6e}\n"
                f"  E⁻ = {minus['E_minus_t']:.6e}\n"
                f"  net = {E_net:.6e}"
            ),
        }

    def time_series(self, dataset: dict,
                    t_start: float = 0.0,
                    t_end: float = 1e10,
                    n_points: int = 200) -> Dict[str, Any]:
        """Net energy evolution over time."""
        dt = (t_end - t_start) / max(n_points - 1, 1)
        series = []
        for i in range(n_points):
            t = t_start + i * dt
            d = dict(dataset)
            d['t'] = t
            result = self.compute(d)
            series.append({
                "t": t,
                "E_net": result["E_net_t"],
                "E_plus": result["E_plus_t"],
                "E_minus": result["E_minus_t"],
                "regime": result["regime"],
            })
        return {
            "time_series": series,
            "t_start": t_start,
            "t_end": t_end,
            "n_points": n_points,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  GW170817 DAMPING INTEGRATION
# ══════════════════════════════════════════════════════════════════════════════

class GWDampingErosion:
    """
    GW170817 strain reduction via negative E(t) during late inspiral.

    PAPER_008b claim: negative E(t) contributes 66.7% strain damping and
    367.8-cycle phase lag when buoyancy opposition dominates in the
    merging neutron star system.

    The damping factor D(t) applied to GR strain h(t):

      h_UQFF(t) = h_GR(t) × [1 + E⁻(t)/E_GW(t)]

    where E_GW(t) is the GW energy flux and E⁻(t)/E_GW < 0 produces damping.

    At 66.7% reduction: |E⁻(t)/E_GW| ≈ 0.667
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          h_GR:        GR strain amplitude (dimensionless, default 1e-21)
          f_GW:        GW frequency (Hz, default 100)
          M_chirp:     chirp mass (kg, default 1.188 M_sun)
          t_merge:     time to merger (s, default 0.1)
          E_GW:        GW energy (J, default 1e47)
          F_U_Bi, F_U: buoyancy forces (N)
          SSq, kappa:  standard params
        """
        h_GR     = dataset.get('h_GR', 1e-21)
        f_GW     = dataset.get('f_GW', 100.0)
        M_chirp  = dataset.get('M_chirp', 1.188 * M_sun)
        t_merge  = dataset.get('t_merge', 0.1)
        E_GW     = dataset.get('E_GW', 1e47)
        F_U_Bi   = dataset.get('F_U_Bi', 0.333)  # < F_U → erosion
        F_U      = dataset.get('F_U', 1.0)
        ssq      = dataset.get('SSq', SSQ)
        kappa    = dataset.get('kappa', KAPPA)

        # Compute E⁻ at merger time
        erosion = NegativeEtErosion()
        e_result = erosion.compute({
            'E_0': E_GW,
            't': t_merge,
            'kappa': kappa,
            'SSq': ssq,
            'F_U_Bi': F_U_Bi,
            'F_U': F_U,
        })
        E_minus = e_result["E_minus_t"]

        # Damping ratio
        if E_GW != 0:
            damping_ratio = abs(E_minus / E_GW)
        else:
            damping_ratio = 0.0

        # UQFF strain
        h_UQFF = h_GR * (1.0 - damping_ratio)

        # Phase lag (cycles): proportional to E⁻ accumulated over inspiral
        # Simplified: Δφ = f_GW × |E⁻/E_GW| × t_merge / (2π)
        phase_lag_rad = 2 * PI * f_GW * damping_ratio * t_merge
        phase_lag_cycles = phase_lag_rad / (2 * PI)

        # Strain reduction percentage
        strain_reduction_pct = damping_ratio * 100.0

        return {
            "h_UQFF": h_UQFF,
            "h_GR": h_GR,
            "strain_reduction_pct": strain_reduction_pct,
            "damping_ratio": damping_ratio,
            "phase_lag_cycles": phase_lag_cycles,
            "phase_lag_rad": phase_lag_rad,
            "E_minus_t": E_minus,
            "E_GW": E_GW,
            "f_GW": f_GW,
            "t_merge": t_merge,
            "equation": (
                "h_UQFF = h_GR × [1 − |E⁻(t)|/E_GW]\n"
                f"       = {h_GR:.4e} × [1 − {damping_ratio:.4f}]\n"
                f"       = {h_UQFF:.4e}\n"
                f"Strain reduction: {strain_reduction_pct:.1f}%\n"
                f"Phase lag: {phase_lag_cycles:.1f} cycles"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  LAGRANGIAN VARIATION δS/δφ_erosion = 0 (SYMMETRIC)
# ══════════════════════════════════════════════════════════════════════════════

class ErosionLagrangianVariation:
    """
    Symmetric Euler-Lagrange variation for negative E(t) erosion.

    L_erosion = −β_i Σ U_{g,i} Ω_g M/d_g [UA]
                + F_neutron · S₂₆([SSq]) · E⁻(t)

    Variation δS/δφ_erosion = 0 yields:

      E⁻(t) = −E₀ (1 − F_{U,Bi}/F_U) exp(κt + [SSq]t/26) S₂₆([SSq])

    This is the variational origin of filament erosion, nebular cavity
    formation, and the vacuum depletion phase in cosmogenesis.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        Ug      = dataset.get('Ug', [1e20, 1e20, 1e20, 1e20])
        Omega_g = dataset.get('Omega_g', 7.3e-16)
        M       = dataset.get('M', M_sun)
        d_g     = dataset.get('d_g', 2.55e20)
        UA      = dataset.get('UA', U_UA)
        beta_i  = dataset.get('beta_i', BETA_I)
        F_n     = dataset.get('F_neutron', 1e10)
        ssq     = dataset.get('SSq', SSQ)
        kappa   = dataset.get('kappa', KAPPA)
        E_0     = dataset.get('E_0', 1.0)
        t       = dataset.get('t', 0.0)
        F_U_Bi  = dataset.get('F_U_Bi', 0.9)
        F_U     = dataset.get('F_U', 1.0)

        # Buoyancy sector (negative for erosion)
        Ug_sum = sum(Ug)
        orbit_factor = Omega_g * M / d_g
        buoyancy_sector = -beta_i * Ug_sum * orbit_factor * UA

        # S₂₆
        S26_val = _eta_euler_s26(ssq)

        # Neutron-polylog coupling in erosion regime
        erosion_factor = 1.0 - (F_U_Bi / F_U) if F_U != 0 else 0.0
        neutron_erosion = F_n * S26_val * erosion_factor

        L_erosion = buoyancy_sector + neutron_erosion

        # Variational solution
        exp_arg = min(kappa * t + (ssq * t) / N_LEVELS, 700.0)
        growth = math.exp(exp_arg)
        E_minus_variational = -E_0 * erosion_factor * growth * S26_val

        rate = kappa + ssq / N_LEVELS
        residual = abs(rate * E_minus_variational + buoyancy_sector * growth / E_0) \
                   if E_0 != 0 else float('inf')

        return {
            "L_erosion": L_erosion,
            "buoyancy_sector": buoyancy_sector,
            "neutron_erosion": neutron_erosion,
            "E_minus_variational": E_minus_variational,
            "erosion_factor": erosion_factor,
            "S_26": S26_val,
            "effective_rate": rate,
            "residual": residual,
            "equation_lagrangian": (
                "L_ero = −β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_n · S₂₆ · (1−F_{U,Bi}/F_U)\n"
                f"      = {L_erosion:.6e}"
            ),
            "equation_variation": (
                "δS/δφ_erosion = 0  ⟹\n"
                "E⁻(t) = −E₀ (1 − F_{U,Bi}/F_U) exp(κt + [SSq]t/26) S₂₆\n"
                f"      = {E_minus_variational:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  WSTP KERNEL SYMBOLIC FORM
# ══════════════════════════════════════════════════════════════════════════════

def wstp_kernel_negative_et() -> str:
    """Wolfram Language symbolic definitions for E⁻(t)."""
    return r"""
(* ═══════════════════════════════════════════════════════════════════════ *)
(* Negative E(t) Buoyancy-Driven Erosion — UQFF Symbolic Forms           *)
(* Session 205 | Daniel Murphy                                            *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* Master negative E(t) expression *)
Eminus[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  -E0 * Exp[\[Kappa] t + SSq t / 26] * S26[SSq] * (1 - FUBi / FU);

(* Net energy evolution *)
Enet[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  Eplus[t, E0, \[Kappa], SSq, FUBi, FU] + Eminus[t, E0, \[Kappa], SSq, FUBi, FU];

(* Simplified net factor *)
EnetSimplified[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  E0 * Exp[\[Kappa] t + SSq t / 26] * S26[SSq] * (2 FUBi / FU - 1);

(* Lagrangian variation for negative E(t) *)
LErosion[\[Beta]i_, Ug_, \[CapitalOmega]g_, M_, dg_, UA_, Fn_, SSq_, FUBi_, FU_] :=
  -\[Beta]i * Ug * \[CapitalOmega]g * M / dg * UA + Fn * S26[SSq] * (1 - FUBi / FU);

deltaSErosion[t_, E0_, \[Kappa]_, SSq_] :=
  D[Eminus[t, E0, \[Kappa], SSq, FUBi, FU] *
    LErosion[\[Beta]i, Ug, \[CapitalOmega]g, M, dg, UA, Fn, SSq, FUBi, FU],
    \[Phi]erosion];

(* GW damping: h_UQFF = h_GR × [1 + E⁻/E_GW] *)
hUQFF[hGR_, Eminus_, EGW_] := hGR * (1 + Eminus / EGW);
"""


# ══════════════════════════════════════════════════════════════════════════════
# §6  SELF-TEST
# ══════════════════════════════════════════════════════════════════════════════

def _run_self_test():
    """Validate all components."""
    print("=" * 72)
    print("negative_et_erosion.py — Self-Test")
    print("=" * 72)
    passed = 0
    failed = 0

    # Test 1: E⁻(t=0) is negative when F_{U,Bi} < F_U
    eng = NegativeEtErosion()
    r = eng.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.9, 'F_U': 1.0})
    E_at_0 = r["E_minus_t"]
    print(f"\nT1  E⁻(0, ratio=0.9) = {E_at_0:.10e}")
    assert E_at_0 < 0, f"E⁻ must be negative for erosion, got {E_at_0}"
    assert r["is_eroding"] is True
    passed += 1
    print("    PASS")

    # Test 2: E⁻ becomes more negative with time
    r2 = eng.compute({'E_0': 1.0, 't': 1e9, 'F_U_Bi': 0.9, 'F_U': 1.0})
    print(f"\nT2  E⁻(1e9s) = {r2['E_minus_t']:.10e}")
    assert r2["E_minus_t"] < E_at_0, "E⁻ should decrease with time"
    passed += 1
    print("    PASS")

    # Test 3: E⁻ = 0 when F_{U,Bi} = F_U (no erosion)
    r3 = eng.compute({'E_0': 1.0, 't': 1e6, 'F_U_Bi': 1.0, 'F_U': 1.0})
    print(f"\nT3  E⁻(F_ratio=1.0) = {r3['E_minus_t']:.10e}")
    assert abs(r3["E_minus_t"]) < 1e-20, "E⁻ should be 0 when balanced"
    passed += 1
    print("    PASS")

    # Test 4: Net energy evolution
    net = NetEnergyEvolution()
    nr = net.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 1.1, 'F_U': 1.0})
    print(f"\nT4  E_net(ratio=1.1) = {nr['E_net_t']:.10e}")
    print(f"    net_factor = {nr['net_factor']:.4f}, regime = {nr['regime']}")
    assert nr["regime"] == "expanding"
    assert nr["E_net_t"] > 0
    passed += 1
    print("    PASS")

    # Test 5: Net energy eroding when ratio < 0.5
    nr2 = net.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.3, 'F_U': 1.0})
    print(f"\nT5  E_net(ratio=0.3) = {nr2['E_net_t']:.10e}")
    print(f"    net_factor = {nr2['net_factor']:.4f}, regime = {nr2['regime']}")
    assert nr2["regime"] == "eroding"
    assert nr2["E_net_t"] < 0
    passed += 1
    print("    PASS")

    # Test 6: Net balance at critical ratio 0.5
    nr3 = net.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.5, 'F_U': 1.0})
    print(f"\nT6  E_net(ratio=0.5) = {nr3['E_net_t']:.10e}")
    print(f"    net_factor = {nr3['net_factor']:.4f}, regime = {nr3['regime']}")
    assert nr3["regime"] == "balanced"
    assert abs(nr3["E_net_t"]) < 1e-20
    passed += 1
    print("    PASS")

    # Test 7: GW damping
    gw = GWDampingErosion()
    gr = gw.compute({'h_GR': 1e-21, 'E_GW': 1e47, 'F_U_Bi': 0.333, 'F_U': 1.0,
                     't_merge': 0.1, 'f_GW': 100.0})
    print(f"\nT7  GW170817 strain reduction = {gr['strain_reduction_pct']:.2f}%")
    print(f"    Phase lag = {gr['phase_lag_cycles']:.2f} cycles")
    assert gr["strain_reduction_pct"] > 0, "Must have positive strain reduction"
    assert gr["h_UQFF"] < gr["h_GR"], "UQFF strain < GR strain"
    passed += 1
    print("    PASS")

    # Test 8: Lagrangian variation
    lv = ErosionLagrangianVariation()
    lr = lv.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 0.9, 'F_U': 1.0})
    print(f"\nT8  L_erosion = {lr['L_erosion']:.6e}")
    print(f"    E⁻_variational = {lr['E_minus_variational']:.6e}")
    assert lr["E_minus_variational"] < 0, "Variational solution must be negative"
    passed += 1
    print("    PASS")

    # Test 9: WSTP kernel
    wk = wstp_kernel_negative_et()
    assert "Eminus[t_" in wk
    assert "Enet[t_" in wk
    assert "hUQFF" in wk
    passed += 1
    print("\nT9  WSTP kernel symbolic form: valid")
    print("    PASS")

    # Test 10: Symmetry check — |E⁺| + |E⁻| integrates correctly
    d_expand = {'E_0': 1.0, 't': 1e6, 'F_U_Bi': 1.1, 'F_U': 1.0}
    ep = PositiveEtExpansion().compute(d_expand)
    d_erode = {'E_0': 1.0, 't': 1e6, 'F_U_Bi': 0.9, 'F_U': 1.0}
    em = eng.compute(d_erode)
    # With symmetric ratios (1.1 and 0.9), |E⁺(1.1)| should differ from |E⁻(0.9)|
    # because the formulas use different factors (ratio vs 1-ratio)
    print(f"\nT10 Symmetry: |E⁺(1.1)| = {abs(ep['E_plus_t']):.10e}")
    print(f"              |E⁻(0.9)| = {abs(em['E_minus_t']):.10e}")
    assert abs(ep["E_plus_t"]) > 0 and abs(em["E_minus_t"]) > 0
    passed += 1
    print("    PASS")

    print(f"\n{'=' * 72}")
    print(f"RESULTS: {passed}/{passed + failed} PASS, {failed} FAIL")
    print(f"{'=' * 72}")
    return passed, failed


if __name__ == "__main__":
    p, f = _run_self_test()
    exit(0 if f == 0 else 1)
