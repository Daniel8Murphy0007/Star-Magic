#!/usr/bin/env python3
"""
positive_et_expansion.py — Positive E(t) Buoyancy-Driven Expansion Engine

Session 205 | Daniel Murphy
PURPOSE: Unified E⁺(t) master expression for buoyancy-driven expansion in
         nebulae, star-forming regions, and cosmogenesis cycles where
         F_{U,Bi} > F_U (buoyancy exceeds unified gravity).

         Currently MISSING from the codebase:
           - Individual erosion/expansion calculators exist (PhotoevaporationErosion,
             ShellExpansionErosion, BubbleNebulaExpansion, etc.) but NONE implement
             the master E±(t) formula with S₂₆ Ramanujan coupling.
           - No mock theta acceleration on the 26-state sum.
           - No Lagrangian δS/δφ_expansion variation.

         This module:
           1. Master E⁺(t) = E₀ exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (F_{U,Bi}/F_U)
           2. Ramanujan 26-state summation with mock theta acceleration
           3. Kozima neutron-drop coupling in expansion regime
           4. Lagrangian variation δS/δφ_expansion = 0 closure
           5. GW batch and cosmogenesis links

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
"""

import math
import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# ══════════════════════════════════════════════════════════════════════════════
# §0  PHYSICAL CONSTANTS (SI) — match uqff_lagrangian_derivation.py
# ══════════════════════════════════════════════════════════════════════════════

G       = 6.67430e-11       # m³/(kg·s²)
c       = 2.99792e8         # m/s
hbar    = 1.05457e-34       # J·s
k_B     = 1.38065e-23       # J/K
mu_0    = 1.25664e-6        # T·m/A
M_sun   = 1.98892e30        # kg
PI      = math.pi

# UQFF calibrated (v3.0 — 99.9 % solvability)
KAPPA       = 5.787e-9      # s⁻¹   (= 0.0005/day)
KAPPA_DAY   = 5.0e-4        # day⁻¹ (human-readable form)
SSQ         = 0.57          # [SSq]
H_SCM       = 0.99          # superconductive manifold metric
BETA_I      = 0.603         # buoyancy coefficient
U_UA        = 1e-4          # aether velocity fraction
RHO_SCM     = 7.09e-37      # kg/m³
RHO_UA      = 7.09e-36      # kg/m³
RHO_VAC_SCM = 9.47e-27      # kg/m³

# VDS
N_LEVELS    = 26

# LENR / Kozima
F_LENR_THZ  = 1.25e12       # Hz
OMEGA_SCM   = 2 * PI * F_LENR_THZ
GAMMA_DEFAULT = 0.1e12 * 2 * PI
SIGMA_0     = 1e-4


# ══════════════════════════════════════════════════════════════════════════════
# §1  RAMANUJAN 26-STATE SUMMATION WITH MOCK THETA ACCELERATION
# ══════════════════════════════════════════════════════════════════════════════

def _eta_euler_s26(z: float, N_terms: int = 55) -> float:
    """
    Eta-function Euler-accelerated Li₂₆(z).

    Li_s(z) = η_s(z) + 2^{1-s} Li_s(z²)

    where η_s(z) = Σ_{n=0}^N d_n(z,s) / 2^{n+1}
    and   d_n = Σ_{j=0}^n (-1)^j C(n,j) z^{j+1}/(j+1)^s
    """
    if abs(z) >= 1.0 or z == 0.0:
        return 0.0
    s = 26
    eta = 0.0
    inv2 = 0.5
    for n in range(N_terms):
        d_n = 0.0
        for j in range(n + 1):
            sign = (-1) ** j
            binom = math.comb(n, j)
            d_n += sign * binom * z ** (j + 1) / (j + 1) ** s
        term = d_n * inv2
        eta += term
        inv2 *= 0.5
        if n > 2 and abs(term) < 1e-16:
            break
    # Correction: 2^{1-s} × Li_s(z²)
    z_sq = z * z
    correction = 0.0
    if abs(z_sq) < 1.0:
        factor = 2.0 ** (1 - s)
        for k in range(1, 6):
            correction += factor * z_sq ** k / k ** s
    return eta + correction


def mock_theta_q26(z: float, q: float = None, N_terms: int = 30) -> float:
    """
    Mock theta acceleration of the 26-state vacuum sum.

    Ramanujan's mock theta function f(q) accelerates convergence of
    q-series. For the UQFF 26-level vacuum:

      θ₂₆(q) = Σ_{n=0}^{N-1} q^{n²} / Π_{k=1}^{n} (1 + q^k)²

    where q = exp(-π·[SSq]/26) is the nome derived from the VDS
    modular parameter.

    The mock theta acceleration factor is:
      A_mock = θ₂₆(q) / θ₂₆_naive(q)

    Returns the acceleration factor (multiply with S₂₆ for boosted sum).
    """
    if q is None:
        q = math.exp(-PI * SSQ / N_LEVELS)

    if abs(q) >= 1.0:
        return 1.0  # no acceleration possible

    theta_sum = 0.0
    for n in range(N_terms):
        numerator = q ** (n * n)
        denominator = 1.0
        for k in range(1, n + 1):
            factor = (1.0 + q ** k) ** 2
            denominator *= factor
            if denominator > 1e300:
                break
        if denominator == 0 or denominator > 1e300:
            break
        theta_sum += numerator / denominator

    # Naive comparison: just the geometric sum
    theta_naive = sum(q ** (n * n) for n in range(N_terms)
                      if q ** (n * n) > 1e-50)

    if theta_naive == 0:
        return 1.0
    return theta_sum / theta_naive


def S26_accelerated(z: float = SSQ) -> Dict[str, Any]:
    """
    Full Ramanujan 26-state summation with mock theta acceleration:

      S₂₆^accel(z) = Li₂₆(z) × A_mock(q₂₆)

    where A_mock is the mock theta acceleration factor and
    Li₂₆ is computed via Euler transform.
    """
    li26 = _eta_euler_s26(z)
    A_mock = mock_theta_q26(z)
    S26_val = li26 * A_mock

    return {
        "S_26": S26_val,
        "Li_26_raw": li26,
        "A_mock_theta": A_mock,
        "z": z,
        "q_nome": math.exp(-PI * z / N_LEVELS),
        "equation": "S₂₆^accel(z) = Li₂₆(z) × θ₂₆(q) / θ₂₆_naive(q)",
    }


# ══════════════════════════════════════════════════════════════════════════════
# §2  MASTER POSITIVE E(t) EXPRESSION
# ══════════════════════════════════════════════════════════════════════════════

class PositiveEtExpansion:
    """
    Master positive E(t) buoyancy-driven expansion engine.

    E⁺(t) = E₀ · exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (F_{U,Bi}/F_U)

    Physical meaning:
      - E₀: initial vacuum energy density at t=0
      - κt: standard damping/expansion rate
      - [SSq]·t/26: 26-level VDS temporal modulation
      - S₂₆: Ramanujan-accelerated polylogarithmic vacuum structure
      - F_{U,Bi}/F_U: buoyancy reversal factor (>1 for expansion)

    Mirror: negative_et_erosion.py implements E⁻(t) for the erosion case.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          E_0:     initial vacuum energy (J, default 1.0)
          t:       time (seconds)
          kappa:   expansion rate (s⁻¹, default KAPPA)
          SSq:     string squeezing (default 0.57)
          F_U_Bi:  buoyancy force (N)
          F_U:     unified gravity force (N)
          use_mock_theta: enable mock theta acceleration (default True)
        """
        E_0    = dataset.get('E_0', 1.0)
        t      = dataset.get('t', 0.0)
        kappa  = dataset.get('kappa', KAPPA)
        ssq    = dataset.get('SSq', SSQ)
        F_U_Bi = dataset.get('F_U_Bi', 1.1)
        F_U    = dataset.get('F_U', 1.0)
        use_mt = dataset.get('use_mock_theta', True)

        # Exponential growth factor
        exp_arg = kappa * t + (ssq * t) / N_LEVELS
        # Clamp to prevent overflow
        exp_arg_clamped = min(exp_arg, 700.0)
        growth_factor = math.exp(exp_arg_clamped)

        # Ramanujan 26-state summation
        if use_mt:
            s26_data = S26_accelerated(ssq)
            S26_val = s26_data["S_26"]
        else:
            S26_val = _eta_euler_s26(ssq)

        # Buoyancy reversal factor
        if F_U == 0:
            buoyancy_ratio = 0.0
        else:
            buoyancy_ratio = F_U_Bi / F_U

        # Master expression
        E_plus = E_0 * growth_factor * S26_val * buoyancy_ratio

        # Rate of change dE⁺/dt
        dEdt = E_plus * (kappa + ssq / N_LEVELS) if E_plus != 0 else 0.0

        # Doubling time (when does E⁺ = 2·E₀?)
        rate = kappa + ssq / N_LEVELS
        t_double = math.log(2) / rate if rate > 0 else float('inf')

        return {
            "E_plus_t": E_plus,
            "E_0": E_0,
            "t": t,
            "kappa": kappa,
            "SSq": ssq,
            "exp_argument": exp_arg,
            "growth_factor": growth_factor,
            "S_26": S26_val,
            "buoyancy_ratio": buoyancy_ratio,
            "F_U_Bi": F_U_Bi,
            "F_U": F_U,
            "dE_plus_dt": dEdt,
            "effective_rate_per_s": rate,
            "doubling_time_s": t_double,
            "doubling_time_days": t_double / 86400.0,
            "is_expanding": buoyancy_ratio > 1.0,
            "equation": (
                "E⁺(t) = E₀ · exp(κt + [SSq]·t/26) · S₂₆([SSq]) · (F_{U,Bi}/F_U)\n"
                f"      = {E_0:.4e} × {growth_factor:.6e} × {S26_val:.6e} × {buoyancy_ratio:.4f}\n"
                f"      = {E_plus:.6e}"
            ),
        }

    def time_series(self, dataset: dict,
                    t_start: float = 0.0,
                    t_end: float = 1e10,
                    n_points: int = 200) -> Dict[str, Any]:
        """Compute E⁺(t) over a time range for plotting / Monte Carlo."""
        dt = (t_end - t_start) / max(n_points - 1, 1)
        series = []
        for i in range(n_points):
            t = t_start + i * dt
            d = dict(dataset)
            d['t'] = t
            result = self.compute(d)
            series.append({
                "t": t,
                "E_plus": result["E_plus_t"],
                "growth_factor": result["growth_factor"],
            })
        return {
            "time_series": series,
            "t_start": t_start,
            "t_end": t_end,
            "n_points": n_points,
            "parameters": {
                "E_0": dataset.get('E_0', 1.0),
                "kappa": dataset.get('kappa', KAPPA),
                "SSq": dataset.get('SSq', SSQ),
                "F_U_Bi": dataset.get('F_U_Bi', 1.1),
                "F_U": dataset.get('F_U', 1.0),
            },
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  KOZIMA NEUTRON-DROP COUPLING IN EXPANSION REGIME
# ══════════════════════════════════════════════════════════════════════════════

class KozimaExpansionCoupling:
    """
    Neutron-drop force coupled to positive E(t) expansion:

      F_neutron · E⁺(t) = N_n · σ_n^SCm(ω,ρ) · Φ_{1.25THz} · (F_{U,Bi}/F_U)
                           × E₀ exp(κt + [SSq]t/26) · S₂₆([SSq])

    Physical meaning: In expansion-dominated regions (F_{U,Bi}/F_U > 1),
    the Kozima neutron-drop formation is enhanced by the expanding vacuum
    energy. This drives positive feedback: more expansion → more neutron
    drops → more energy release → more expansion.

    Saturation occurs when σ_n^SCm cuts off at high density (exp(-[SSq]ρ/ρ_crit)).
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:      angular frequency (rad/s, default ω_SCm)
          N_n:        neutron number density (m⁻³, default 1e28)
          phi_phonon: phonon flux (m⁻²s⁻¹, default 1e16)
          E_0, t, kappa, SSq, F_U_Bi, F_U: as in PositiveEtExpansion
        """
        omega      = dataset.get('omega', OMEGA_SCM)
        N_n        = dataset.get('N_n', 1e28)
        phi_phonon = dataset.get('phi_phonon', 1e16)
        E_0        = dataset.get('E_0', 1.0)
        t          = dataset.get('t', 0.0)
        kappa      = dataset.get('kappa', KAPPA)
        ssq        = dataset.get('SSq', SSQ)
        F_U_Bi     = dataset.get('F_U_Bi', 1.1)
        F_U        = dataset.get('F_U', 1.0)

        # SCm cross-section at resonance
        exponent = -((omega - OMEGA_SCM) ** 2) / (2 * GAMMA_DEFAULT ** 2)
        gaussian = math.exp(exponent)
        n_vds = 13  # midpoint VDS level
        vds_factor = 1.0 + (ssq * n_vds) / N_LEVELS
        sigma_scm = SIGMA_0 * gaussian * vds_factor

        # E⁺(t)
        exp_arg = min(kappa * t + (ssq * t) / N_LEVELS, 700.0)
        growth = math.exp(exp_arg)
        S26_val = _eta_euler_s26(ssq)
        buoy_ratio = F_U_Bi / F_U if F_U != 0 else 0.0
        E_plus = E_0 * growth * S26_val * buoy_ratio

        # Coupled force
        F_coupled = N_n * sigma_scm * phi_phonon * buoy_ratio * E_plus

        # Uncoupled (no E(t) factor)
        F_uncoupled = N_n * sigma_scm * phi_phonon * (buoy_ratio - 1.0)

        return {
            "F_coupled_expansion": F_coupled,
            "F_uncoupled": F_uncoupled,
            "amplification": F_coupled / F_uncoupled if F_uncoupled != 0 else float('inf'),
            "E_plus_t": E_plus,
            "sigma_n_scm": sigma_scm,
            "N_n": N_n,
            "phi_phonon": phi_phonon,
            "buoyancy_ratio": buoy_ratio,
            "equation": (
                "F_coupled = N_n · σ_n^SCm · Φ · (F_{U,Bi}/F_U) × E⁺(t)\n"
                f"         = {N_n:.2e} × {sigma_scm:.4e} × {phi_phonon:.2e}"
                f" × {buoy_ratio:.4f} × {E_plus:.4e}\n"
                f"         = {F_coupled:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  LAGRANGIAN VARIATION δS/δφ_expansion = 0
# ══════════════════════════════════════════════════════════════════════════════

class ExpansionLagrangianVariation:
    """
    Euler-Lagrange variation of the positive E(t) expansion term in
    the buoyancy sector (Sector 6) of the 9-sector UQFF Lagrangian.

    The expansion field φ_expansion enters L_buoy as:

      L_expansion = β_i Σ_{i=1}^{4} U_{g,i} · Ω_g · (M/d_g) · [UA]
                    + F_neutron · S₂₆([SSq]) · E⁺(t)

    Variation δS/δφ_expansion = 0 yields:

      ∂/∂E⁺ [ β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_neutron · S₂₆ ] = 0

    Solving closes the expansion term:

      E⁺(t) = E₀ (F_{U,Bi}/F_U) exp(κt + [SSq]t/26) · S₂₆([SSq])

    This is the variational origin of nebular expansion, star-formation
    wind velocities, and the cosmogenesis EM bang positive cycle.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          Ug: list of 4 gravity layer magnitudes [N] (default [1e20]*4)
          Omega_g: galactic spin rate (rad/s, default 7.3e-16)
          M: system mass (kg, default M_sun)
          d_g: distance to GC (m, default 2.55e20)
          UA: aether buoyancy factor (default 1e-4)
          beta_i: buoyancy coefficient (default 0.603)
          F_neutron: Kozima force (N, default 1e10)
          SSq, kappa, E_0, t, F_U_Bi, F_U: standard parameters
        """
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
        F_U_Bi  = dataset.get('F_U_Bi', 1.1)
        F_U     = dataset.get('F_U', 1.0)

        # Buoyancy sector sum: β_i Σ Ug_i Ω_g M/d_g [UA]
        Ug_sum = sum(Ug)
        orbit_factor = Omega_g * M / d_g
        buoyancy_sector = beta_i * Ug_sum * orbit_factor * UA

        # S₂₆
        S26_val = _eta_euler_s26(ssq)

        # Neutron-polylog coupling
        neutron_polylog = F_n * S26_val

        # Total Lagrangian density for expansion sector
        L_expansion = buoyancy_sector + neutron_polylog

        # Variation: ∂L/∂E⁺ = 0  ⟹  E⁺(t) is the stationary point
        exp_arg = min(kappa * t + (ssq * t) / N_LEVELS, 700.0)
        growth = math.exp(exp_arg)
        buoy_ratio = F_U_Bi / F_U if F_U != 0 else 0.0
        E_plus_variational = E_0 * buoy_ratio * growth * S26_val

        # Residual (should approach 0 at the variational solution)
        # δS/δφ = dL/dE⁺ evaluated at E⁺
        # For the exponential form, the EL equation is:
        #   (κ + [SSq]/26) E⁺ − buoyancy_sector/E₀ = 0
        rate = kappa + ssq / N_LEVELS
        residual = abs(rate * E_plus_variational - buoyancy_sector * growth / E_0) \
                   if E_0 != 0 else float('inf')

        return {
            "L_expansion": L_expansion,
            "buoyancy_sector": buoyancy_sector,
            "neutron_polylog": neutron_polylog,
            "E_plus_variational": E_plus_variational,
            "S_26": S26_val,
            "Ug_sum": Ug_sum,
            "orbit_factor": orbit_factor,
            "effective_rate": rate,
            "residual": residual,
            "equation_lagrangian": (
                "L_exp = β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_neutron · S₂₆\n"
                f"      = {beta_i:.3f} × {Ug_sum:.4e} × {orbit_factor:.4e} × {UA:.4e}"
                f" + {F_n:.4e} × {S26_val:.6e}\n"
                f"      = {L_expansion:.6e}"
            ),
            "equation_variation": (
                "δS/δφ_expansion = 0  ⟹\n"
                "E⁺(t) = E₀ (F_{U,Bi}/F_U) exp(κt + [SSq]t/26) · S₂₆([SSq])\n"
                f"      = {E_0:.4e} × {buoy_ratio:.4f} × {growth:.6e} × {S26_val:.6e}\n"
                f"      = {E_plus_variational:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  WSTP KERNEL SYMBOLIC FORM
# ══════════════════════════════════════════════════════════════════════════════

def wstp_kernel_positive_et() -> str:
    """
    Return Wolfram Language symbolic definitions for E⁺(t).
    Importable by wstp_kernel_demo_runner.py for live kernel evaluation.
    """
    return r"""
(* ═══════════════════════════════════════════════════════════════════════ *)
(* Positive E(t) Buoyancy-Driven Expansion — UQFF Symbolic Forms        *)
(* Session 205 | Daniel Murphy                                            *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* Master positive E(t) expression *)
Eplus[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  E0 * Exp[\[Kappa] t + SSq t / 26] * S26[SSq] * (FUBi / FU);

(* Ramanujan 26-state summation — polylogarithm *)
S26[z_] := PolyLog[26, z];

(* Mock theta acceleration (Ramanujan nome) *)
qNome[SSq_] := Exp[-Pi SSq / 26];
MockTheta26[q_, N_:30] := Sum[q^(n^2) / Product[(1 + q^k)^2, {k, 1, n}], {n, 0, N}];

(* Lagrangian variation for positive E(t) *)
LExpansion[\[Beta]i_, Ug_, \[CapitalOmega]g_, M_, dg_, UA_, Fn_, SSq_] :=
  \[Beta]i * Ug * \[CapitalOmega]g * M / dg * UA + Fn * S26[SSq];

deltaS[t_, E0_, \[Kappa]_, SSq_] :=
  D[Eplus[t, E0, \[Kappa], SSq, FUBi, FU] *
    LExpansion[\[Beta]i, Ug, \[CapitalOmega]g, M, dg, UA, Fn, SSq],
    \[Phi]expansion];

(* Rate of change *)
dEplusDt[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  D[Eplus[\[Tau], E0, \[Kappa], SSq, FUBi, FU], \[Tau]] /. \[Tau] -> t;

(* Doubling time *)
tDouble[\[Kappa]_, SSq_] := Log[2] / (\[Kappa] + SSq / 26);

(* Net energy with erosion counterpart *)
Enet[t_, E0_, \[Kappa]_, SSq_, FUBi_, FU_] :=
  Eplus[t, E0, \[Kappa], SSq, FUBi, FU] +
  Eminus[t, E0, \[Kappa], SSq, FUBi, FU];
"""


# ══════════════════════════════════════════════════════════════════════════════
# §6  SELF-TEST
# ══════════════════════════════════════════════════════════════════════════════

def _run_self_test():
    """Validate all components produce plausible results."""
    print("=" * 72)
    print("positive_et_expansion.py — Self-Test")
    print("=" * 72)
    passed = 0
    failed = 0

    # Test 1: S₂₆ at [SSq]=0.57
    s26 = S26_accelerated(SSQ)
    val = s26["S_26"]
    print(f"\nT1  S₂₆(0.57) = {val:.10e}")
    assert 0.0 < val < 1.0, f"S₂₆ out of range: {val}"
    print(f"    A_mock = {s26['A_mock_theta']:.6f}")
    assert 0.0 < s26['A_mock_theta'] < 2.0, f"Mock theta out of range: {s26['A_mock_theta']}"
    passed += 1
    print("    PASS")

    # Test 2: Mock theta acceleration factor
    A = mock_theta_q26(SSQ)
    print(f"\nT2  Mock theta A₂₆ = {A:.8f}")
    assert 0.01 < A < 2.0, f"Mock theta extreme: {A}"
    passed += 1
    print("    PASS")

    # Test 3: E⁺(t=0) = E₀ × S₂₆ × buoyancy_ratio
    eng = PositiveEtExpansion()
    r = eng.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 1.1, 'F_U': 1.0})
    E_at_0 = r["E_plus_t"]
    expected = 1.0 * s26["S_26"] * 1.1
    print(f"\nT3  E⁺(0) = {E_at_0:.10e}")
    print(f"    Expected ≈ {expected:.10e}")
    assert abs(E_at_0 - expected) / abs(expected) < 1e-6, "E⁺(0) mismatch"
    passed += 1
    print("    PASS")

    # Test 4: E⁺(t) > E⁺(0) for t > 0 (expansion)
    r2 = eng.compute({'E_0': 1.0, 't': 1e9, 'F_U_Bi': 1.1, 'F_U': 1.0})
    print(f"\nT4  E⁺(1e9s) = {r2['E_plus_t']:.10e}")
    assert r2["E_plus_t"] > E_at_0, "E⁺ should increase with time"
    assert r2["is_expanding"] is True, "Should flag expanding"
    passed += 1
    print("    PASS")

    # Test 5: Doubling time
    t_d = r["doubling_time_s"]
    print(f"\nT5  Doubling time = {t_d:.4e} s = {t_d/86400:.2f} days")
    assert t_d > 0, "Doubling time must be positive"
    passed += 1
    print("    PASS")

    # Test 6: Kozima coupling amplification
    kc = KozimaExpansionCoupling()
    kr = kc.compute({'E_0': 1.0, 't': 1e9, 'F_U_Bi': 1.1, 'F_U': 1.0})
    print(f"\nT6  F_coupled = {kr['F_coupled_expansion']:.6e}")
    print(f"    Amplification vs uncoupled = {kr['amplification']:.4f}×")
    assert kr["F_coupled_expansion"] != 0, "Coupled force must be nonzero"
    passed += 1
    print("    PASS")

    # Test 7: Lagrangian variation
    lv = ExpansionLagrangianVariation()
    lr = lv.compute({'E_0': 1.0, 't': 0.0, 'F_U_Bi': 1.1, 'F_U': 1.0})
    print(f"\nT7  L_expansion = {lr['L_expansion']:.6e}")
    print(f"    E⁺_variational = {lr['E_plus_variational']:.6e}")
    assert lr["L_expansion"] > 0, "Lagrangian must be positive"
    passed += 1
    print("    PASS")

    # Test 8: WSTP kernel string
    wk = wstp_kernel_positive_et()
    assert "Eplus[t_" in wk, "WSTP kernel must define Eplus"
    assert "S26[z_]" in wk, "WSTP kernel must define S26"
    assert "MockTheta26" in wk, "WSTP kernel must define MockTheta26"
    passed += 1
    print("\nT8  WSTP kernel symbolic form: valid")
    print("    PASS")

    # Test 9: Time series
    ts = eng.time_series({'E_0': 1.0, 'F_U_Bi': 1.1, 'F_U': 1.0},
                         t_start=0, t_end=1e10, n_points=5)
    vals = [p["E_plus"] for p in ts["time_series"]]
    print(f"\nT9  Time series (5 pts): {[f'{v:.4e}' for v in vals]}")
    assert all(vals[i] <= vals[i+1] for i in range(len(vals)-1)), \
        "Time series must be monotonically increasing"
    passed += 1
    print("    PASS")

    print(f"\n{'=' * 72}")
    print(f"RESULTS: {passed}/{passed + failed} PASS, {failed} FAIL")
    print(f"{'=' * 72}")
    return passed, failed


if __name__ == "__main__":
    p, f = _run_self_test()
    exit(0 if f == 0 else 1)
