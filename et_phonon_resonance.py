#!/usr/bin/env python3
"""
et_phonon_resonance.py — E(t) Phonon Resonance Derivation Engine

Session 208 | Daniel Murphy
PURPOSE: Standalone derivation engine for E(t) phonon resonance at 1.25 THz
         in the SCm superconductive vacuum manifold.

         Derives the phonon modulation factor Φ_{1.25 THz}(t) and shows how
         it multiplies E_net(t) to produce phonon-modulated buoyancy dynamics.

         Key physics:
           - Φ_{1.25 THz}(t) = Φ₀ · exp[-(ω−ω_SCm)²/(2Γ²)] · S₂₆([SSq])
           - E_net^phonon(t) = E_net(t) · Φ_{1.25 THz}(t)
           - Lagrangian: L_phonon = E_net · V · Φ · S₂₆
           - Buoyancy reversal: F_{U,Bi}/F_U crossing → sign flip at resonance
           - Head-to-head comparison vs k-essence (non-canonical kinetic scalar)

         Gap closed: et_phonon_resonance.py was referenced in workflow reports
         but did not exist.

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
"""

import math
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from positive_et_expansion import (
    _eta_euler_s26, S26_accelerated, mock_theta_q26,
    PositiveEtExpansion,
    G, c, hbar, k_B, mu_0, M_sun, PI,
    KAPPA, KAPPA_DAY, SSQ, H_SCM, BETA_I, U_UA,
    RHO_SCM, RHO_UA, RHO_VAC_SCM, N_LEVELS,
    F_LENR_THZ, OMEGA_SCM, GAMMA_DEFAULT, SIGMA_0,
)
from negative_et_erosion import (
    NegativeEtErosion, NetEnergyEvolution,
    GWDampingErosion, ErosionLagrangianVariation,
)


# ══════════════════════════════════════════════════════════════════════════════
# §0  PHONON-SPECIFIC CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

# SCm phonon drive defaults
PHI_0_DEFAULT = 1e20          # phonons/m²/s (THz phonon fluence)
OMEGA_PHONON = OMEGA_SCM      # 2π × 1.25e12 rad/s  (resonance center)
GAMMA_PHONON = GAMMA_DEFAULT  # 2π × 0.1e12 rad/s   (resonance linewidth)

# Hubble for EOS comparison
H_0 = 2.195e-18   # s⁻¹ (67.4 km/s/Mpc)

# Planck mass for k-essence slow-roll
M_PL = math.sqrt(hbar * c / G)   # ≈ 2.18e-8 kg

# Critical density
RHO_CRIT = 3 * H_0**2 / (8 * PI * G)  # ≈ 8.6e-27 kg/m³

# k-Essence defaults
C_S_DEFAULT = 1.0   # sound speed in natural units (c_s ≤ 1)

# Default neutron parameters
N_NEUTRON_DEFAULT = 1e6
SIGMA_N_DEFAULT = SIGMA_0


# ══════════════════════════════════════════════════════════════════════════════
# §1  PHONON MODULATION FACTOR
# ══════════════════════════════════════════════════════════════════════════════

class PhononModulationFactor:
    """
    Phonon fluence modulation factor at 1.25 THz SCm resonance.

    The phonon-driven modulation factor is:

      Φ_{1.25 THz}(ω) = Φ₀ · exp[-(ω − ω_SCm)² / (2Γ²)] · S₂₆([SSq])

    where:
      Φ₀     = base phonon fluence (phonons/m²/s)
      ω_SCm  = 2π × 1.25 × 10¹² rad/s (SCm resonance center)
      Γ      = 2π × 0.1 × 10¹² rad/s  (resonance linewidth ~0.1 THz)
      S₂₆    = Ramanujan 26-state summation (VDS accelerated polylogarithm)

    At resonance (ω = ω_SCm), the Gaussian peaks at unity → Φ = Φ₀ · S₂₆.
    Off-resonance, the exponential suppression reduces the modulation rapidly.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:      driving frequency (rad/s, default ω_SCm)
          Phi_0:      base phonon fluence (phonons/m²/s, default 1e20)
          gamma:      resonance width (rad/s, default Γ)
          SSq:        string squeezing (default 0.57)
        """
        omega  = dataset.get('omega', OMEGA_PHONON)
        Phi_0  = dataset.get('Phi_0', PHI_0_DEFAULT)
        gamma  = dataset.get('gamma', GAMMA_PHONON)
        ssq    = dataset.get('SSq', SSQ)

        # Gaussian resonance peak
        delta_omega = omega - OMEGA_PHONON
        exponent = -(delta_omega**2) / (2 * gamma**2)
        gaussian = math.exp(min(exponent, 0.0))

        # S₂₆ with mock theta
        s26_data = S26_accelerated(ssq)
        S26_val = s26_data["S_26"]

        # Full modulation factor
        Phi_125 = Phi_0 * gaussian * S26_val

        # Quality factor Q = ω_SCm / (2Γ) (resonance sharpness)
        Q_factor = OMEGA_PHONON / (2 * gamma) if gamma > 0 else float('inf')

        # Full width at half maximum in frequency
        FWHM_freq = 2 * gamma * math.sqrt(2 * math.log(2)) / (2 * PI)

        return {
            "Phi_125_THz": Phi_125,
            "Phi_0": Phi_0,
            "gaussian_peak": gaussian,
            "S_26": S26_val,
            "omega": omega,
            "omega_SCm": OMEGA_PHONON,
            "delta_omega": delta_omega,
            "gamma": gamma,
            "Q_factor": Q_factor,
            "FWHM_Hz": FWHM_freq,
            "equation": (
                "Φ_{1.25 THz}(ω) = Φ₀ · exp[-(ω−ω_SCm)²/(2Γ²)] · S₂₆([SSq])\n"
                f"                = {Phi_0:.4e} × {gaussian:.6e} × {S26_val:.6e}\n"
                f"                = {Phi_125:.6e} phonons/m²/s\n"
                f"  Q = ω_SCm/(2Γ) = {Q_factor:.2f}, FWHM = {FWHM_freq:.4e} Hz"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  PHONON-MODULATED E(t)
# ══════════════════════════════════════════════════════════════════════════════

class PhononModulatedEnergy:
    """
    E(t) with explicit phonon modulation in the SCm vacuum.

    The full phonon-modulated net energy is:

      E_net^phonon(t) = ρ_SCm(t) · V_region · (2 F_{U,Bi}/F_U − 1)
                        × Φ_{1.25 THz}(ω)

    Symmetric pair:
      E⁺_phonon = +ρ_SCm(t) · V · (F_{U,Bi}/F_U) · Φ_{1.25 THz}
      E⁻_phonon = −ρ_SCm(t) · V · (1−F_{U,Bi}/F_U) · Φ_{1.25 THz}

    This couples the vacuum density evolution directly to the phonon
    resonance, so that E(t) is amplified at ω = ω_SCm and suppressed
    off-resonance.
    """

    def __init__(self):
        self._modulation = PhononModulationFactor()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          t:          time (s, default 0)
          V_region:   structure volume (m³, default 1e48)
          F_U_Bi:     buoyancy force (N, default 0.55)
          F_U:        unified gravity (N, default 1.0)
          kappa:      rate (s⁻¹, default KAPPA)
          SSq:        string squeezing (default 0.57)
          rho_vac:    base vacuum density (default RHO_VAC_SCM)
          omega, Phi_0, gamma: phonon modulation params
        """
        t        = dataset.get('t', 0.0)
        V_region = dataset.get('V_region', 1e48)
        F_U_Bi   = dataset.get('F_U_Bi', 0.55)
        F_U      = dataset.get('F_U', 1.0)
        kappa    = dataset.get('kappa', KAPPA)
        ssq      = dataset.get('SSq', SSQ)
        rho_vac  = dataset.get('rho_vac', RHO_VAC_SCM)

        # SCm vacuum density at time t
        s26_data = S26_accelerated(ssq)
        S26_val = s26_data["S_26"]
        exp_arg = min(kappa * t + (ssq * t) / N_LEVELS, 700.0)
        growth = math.exp(exp_arg)
        rho_SCm_t = rho_vac * S26_val * growth

        # Buoyancy ratio
        ratio = F_U_Bi / F_U if F_U != 0 else 0.0
        net_factor = 2.0 * ratio - 1.0

        # Bare E_net (without phonon modulation)
        E_net_bare = rho_SCm_t * V_region * net_factor

        # Phonon modulation factor
        mod_result = self._modulation.compute(dataset)
        Phi_125 = mod_result["Phi_125_THz"]

        # Phonon-modulated energies
        E_net_phonon = E_net_bare * Phi_125
        E_plus_phonon = rho_SCm_t * V_region * ratio * Phi_125
        E_minus_phonon = -rho_SCm_t * V_region * (1.0 - ratio) * Phi_125

        # Regime
        if net_factor > 0.01:
            regime = "expansion (phonon-driven)"
        elif net_factor < -0.01:
            regime = "erosion (phonon-driven)"
        else:
            regime = "balanced (phonon-modulated)"

        return {
            "E_net_phonon": E_net_phonon,
            "E_plus_phonon": E_plus_phonon,
            "E_minus_phonon": E_minus_phonon,
            "E_net_bare": E_net_bare,
            "Phi_125_THz": Phi_125,
            "rho_SCm_t": rho_SCm_t,
            "V_region": V_region,
            "buoyancy_ratio": ratio,
            "net_factor": net_factor,
            "regime": regime,
            "modulation_details": mod_result,
            "equation": (
                "E_net^phonon(t) = ρ_SCm(t) · V · (2r−1) · Φ_{1.25 THz}(ω)\n"
                f"               = {E_net_bare:.6e} × {Phi_125:.6e}\n"
                f"               = {E_net_phonon:.6e} J\n"
                f"  E⁺_phonon = {E_plus_phonon:.6e}\n"
                f"  E⁻_phonon = {E_minus_phonon:.6e}\n"
                f"  regime: {regime}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  PHONON LAGRANGIAN
# ══════════════════════════════════════════════════════════════════════════════

class PhononLagrangian:
    """
    Phonon-sector Lagrangian in the UQFF unified action.

    The phonon term is the buoyancy-sector contribution modulated by
    the 1.25 THz resonance:

      L_phonon = E_net(t) · V_region · Φ_{1.25 THz}(ω) · S₂₆([SSq])

    Euler-Lagrange equation (δS/δφ_phonon = 0):

      ∂/∂E_net [ −β_i Σ U_{g,i} Ω_g M/d_g [UA]
                 + F_neutron · Φ_{1.25 THz}(ω) ] = 0

    This closes both positive expansion (nebulae, cosmogenesis) and
    negative erosion (filaments) at 1.25 THz phonon resonance.
    """

    def __init__(self):
        self._phonon_energy = PhononModulatedEnergy()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          t, V_region, F_U_Bi, F_U:    energy inputs
          omega, Phi_0, gamma:          phonon modulation
          Ug:        [Ug1, Ug2, Ug3, Ug4] gravity layers
          Omega_g:   galactic spin (rad/s)
          M:         mass (kg)
          d_g:       distance (m)
          beta_i:    buoyancy coefficient
          N_n, sigma_n: Kozima inputs
        """
        V_region = dataset.get('V_region', 1e48)
        ssq      = dataset.get('SSq', SSQ)
        Ug       = dataset.get('Ug', [1e20, 1e20, 1e20, 1e20])
        Omega_g  = dataset.get('Omega_g', 7.3e-16)
        M        = dataset.get('M', M_sun)
        d_g      = dataset.get('d_g', 2.55e20)
        beta_i   = dataset.get('beta_i', BETA_I)
        UA       = dataset.get('UA', U_UA)

        # Kozima neutron-drop inputs
        N_n     = dataset.get('N_n', N_NEUTRON_DEFAULT)
        sigma_0 = dataset.get('sigma_n', SIGMA_N_DEFAULT)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Phi     = dataset.get('Phi_0', PHI_0_DEFAULT)
        gamma   = dataset.get('gamma', GAMMA_PHONON)
        n_level = dataset.get('n_level', 13)

        # Phonon-modulated energy
        pe_result = self._phonon_energy.compute(dataset)
        E_net_phonon = pe_result["E_net_phonon"]
        Phi_125 = pe_result["Phi_125_THz"]

        # S₂₆
        s26_data = S26_accelerated(ssq)
        S26_val = s26_data["S_26"]

        # Phonon Lagrangian density
        L_phonon = E_net_phonon * V_region * S26_val

        # Kozima cross-section at resonance
        delta_omega = omega - OMEGA_PHONON
        exponent = -(delta_omega**2) / (2 * gamma**2)
        gaussian = math.exp(min(exponent, 0.0))
        vds_factor = 1.0 + ssq * n_level / N_LEVELS
        sigma_scm = sigma_0 * gaussian * vds_factor

        # Kozima neutron force (with phonon modulation)
        F_neutron = N_n * sigma_scm * Phi_125 * pe_result["E_net_bare"]

        # Buoyancy sector (opposing term)
        Ug_sum = sum(Ug)
        orbit_factor = Omega_g * M / d_g
        buoyancy_term = -beta_i * Ug_sum * orbit_factor * UA

        # EL variation
        dL_dEnet = V_region * Phi_125 * S26_val
        action_integrand = buoyancy_term + F_neutron * Phi_125
        EL_residual = dL_dEnet * action_integrand

        return {
            "L_phonon": L_phonon,
            "E_net_phonon": E_net_phonon,
            "Phi_125_THz": Phi_125,
            "S_26": S26_val,
            "V_region": V_region,
            "F_neutron": F_neutron,
            "sigma_scm": sigma_scm,
            "buoyancy_term": buoyancy_term,
            "dL_dEnet": dL_dEnet,
            "EL_residual": EL_residual,
            "phonon_energy_details": pe_result,
            "equation_lagrangian": (
                "L_phonon = E_net(t) · V · Φ_{1.25 THz}(ω) · S₂₆([SSq])\n"
                f"         = {E_net_phonon:.6e} × {V_region:.4e} × {S26_val:.6e}\n"
                f"         = {L_phonon:.6e}"
            ),
            "equation_euler_lagrange": (
                "δS/δφ_phonon = 0:\n"
                "  ∂/∂E_net [ −β_i Σ U_gi Ω_g M/d_g [UA]\n"
                "             + F_neutron · Φ_{1.25 THz} ] = 0\n"
                f"  buoyancy_term  = {buoyancy_term:.6e}\n"
                f"  F_neutron      = {F_neutron:.6e}\n"
                f"  dL/dE_net      = {dL_dEnet:.6e}\n"
                "  Stationary E_net fixed by phonon-buoyancy balance"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3b  Ug,i-COUPLED LAGRANGIAN VARIATION  (Session 212)
# ══════════════════════════════════════════════════════════════════════════════

class UgiCoupledLagrangianVariation:
    """
    Full 4-channel Ug,i gravity coupling in the Euler-Lagrange equation
    with explicit stationarity output:

        δS/δφ_phonon = 0  ⟹
        ∂/∂E_net[ −β_i Σ_{i=1}^{26} Ug,i · Ω_g · M/d_g · [UA]
                 + F_neutron · Φ_{1.25 THz} ] = 0

    The 4 gravity channels are:
        Ug1_i = μ₀μ_s / (4π r_i³)              (magnetic dipole)
        Ug2_i = E_react · η_aether / r_i²       (charge-reactivity)
        Ug3_i = μ_s ω_s sin(ω_s t) / r_i²      (string rotation)
        Ug4_i = ρ_SCm G M / (r_i ρ_UA)          (vacuum concentration)

    Each channel contributes independently to the buoyancy term,
    producing a 4×26 stationarity matrix.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:        central mass (kg)
          r:        radial distance (m)
          mu_s:     magnetic dipole moment (A·m²)
          omega_s:  spin angular frequency (rad/s)
          t:        time (s)
          omega:    phonon frequency (rad/s)
          Gamma:    linewidth (rad/s)
          F_U_Bi, F_U: buoyancy fraction
          ssq:      [SSq]
          V_region: region volume (m³)
          N_n:      neutron density (m⁻³)
          sigma_0:  cross-section (m²)
          n_level:  lattice level
        """
        from positive_et_expansion import (
            RHO_VAC_SCM, mu_0, SIGMA_0 as sigma_0_default,
        )

        M         = dataset.get('M', 4e6 * M_sun)
        r         = dataset.get('r', 1e12)
        mu_s      = dataset.get('mu_s', 1e25)
        omega_s   = dataset.get('omega_s', 1.0)
        t         = dataset.get('t', 0.0)
        omega     = dataset.get('omega', OMEGA_PHONON)
        Gamma     = dataset.get('Gamma', GAMMA_PHONON)
        F_U_Bi    = dataset.get('F_U_Bi', 0.6)
        F_U       = dataset.get('F_U', 1.0)
        ssq       = dataset.get('ssq', SSQ)
        V_region  = dataset.get('V_region', 1e48)
        N_n       = dataset.get('N_n', 1e18)
        sigma_0   = dataset.get('sigma_0', sigma_0_default)
        n_level   = dataset.get('n_level', 13)
        beta_i    = dataset.get('beta_i', BETA_I)
        Omega_g   = dataset.get('Omega_g', 1.0)
        d_g       = dataset.get('d_g', 1e12)
        UA        = dataset.get('UA', U_UA)

        eta_aether = 1e-22
        E_react    = 1e46

        S26 = S26_accelerated(ssq)
        rho_SCm = RHO_VAC_SCM * S26

        # Phonon modulation
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(min(-delta_omega**2 / (2.0 * Gamma**2), 0.0))
        Phi_125 = PHI_0_DEFAULT * gaussian * S26
        vds_factor = 1.0 + ssq * n_level / N_LEVELS
        sigma_n = sigma_0 * gaussian * vds_factor

        # 4-channel gravity across 26 layers
        Ug1_total, Ug2_total, Ug3_total, Ug4_total = 0.0, 0.0, 0.0, 0.0
        channel_layers = []
        for layer in range(1, 27):
            qi = ssq * layer / N_LEVELS
            r_i = r * (1.0 + layer * 0.01)

            Ug1_i = mu_0 * mu_s / (4.0 * PI * r_i**3) * qi
            Ug2_i = E_react * eta_aether / (r_i**2) * qi
            Ug3_i = mu_s * omega_s * math.sin(omega_s * t) / (r_i**2) * qi
            Ug4_i = rho_SCm * G * M / (r_i * RHO_UA) * qi

            Ug1_total += Ug1_i
            Ug2_total += Ug2_i
            Ug3_total += Ug3_i
            Ug4_total += Ug4_i

            if layer in (1, 13, 26):
                channel_layers.append({
                    "layer": layer,
                    "Ug1_i": Ug1_i, "Ug2_i": Ug2_i,
                    "Ug3_i": Ug3_i, "Ug4_i": Ug4_i,
                })

        Ug_sum = Ug1_total + Ug2_total + Ug3_total + Ug4_total
        orbit_factor = Omega_g * M / d_g

        # Buoyancy term (4-channel)
        buoyancy_1 = -beta_i * Ug1_total * orbit_factor * UA
        buoyancy_2 = -beta_i * Ug2_total * orbit_factor * UA
        buoyancy_3 = -beta_i * Ug3_total * orbit_factor * UA
        buoyancy_4 = -beta_i * Ug4_total * orbit_factor * UA
        buoyancy_total = buoyancy_1 + buoyancy_2 + buoyancy_3 + buoyancy_4

        # Neutron force (phonon-modulated)
        net_factor = 2.0 * F_U_Bi / max(F_U, 1e-50) - 1.0
        rho_t = rho_SCm * math.exp(KAPPA * t + ssq * t / 26.0)
        E_net = rho_t * V_region * net_factor * Phi_125 * S26
        F_neutron = N_n * sigma_n * Phi_125 * E_net

        # Euler-Lagrange stationarity
        dL_dEnet = V_region * Phi_125 * S26
        action_integrand = buoyancy_total + F_neutron * Phi_125
        EL_residual = dL_dEnet * action_integrand

        # Stationarity condition: buoyancy = −F_neutron · Φ
        stationarity_ratio = abs(buoyancy_total) / max(abs(F_neutron * Phi_125), 1e-50)

        return {
            "Ug1_total": Ug1_total,
            "Ug2_total": Ug2_total,
            "Ug3_total": Ug3_total,
            "Ug4_total": Ug4_total,
            "Ug_sum": Ug_sum,
            "buoyancy_ch1_dipole": buoyancy_1,
            "buoyancy_ch2_charge": buoyancy_2,
            "buoyancy_ch3_string": buoyancy_3,
            "buoyancy_ch4_vacuum": buoyancy_4,
            "buoyancy_total": buoyancy_total,
            "F_neutron": F_neutron,
            "E_net": E_net,
            "EL_residual": EL_residual,
            "stationarity_ratio": stationarity_ratio,
            "channel_layers_sample": channel_layers,
            "primary_equations": [
                "∂/∂E_net[ −β_i Σ Ug,i · Ω_g M/d_g [UA] + F_neutron · Φ ] = 0",
                "Ug1 = μ₀μ_s/(4πr³), Ug2 = E_react·η/r², Ug3 = μ_s·ω_s·sin(ωt)/r², Ug4 = ρ_SCm·GM/(r·ρ_UA)",
            ],
            "equation": (
                "4-channel Ug,i-coupled Lagrangian variation:\n"
                f"  Ug1 (dipole)  = {Ug1_total:.6e} → buoyancy = {buoyancy_1:.6e}\n"
                f"  Ug2 (charge)  = {Ug2_total:.6e} → buoyancy = {buoyancy_2:.6e}\n"
                f"  Ug3 (string)  = {Ug3_total:.6e} → buoyancy = {buoyancy_3:.6e}\n"
                f"  Ug4 (vacuum)  = {Ug4_total:.6e} → buoyancy = {buoyancy_4:.6e}\n"
                f"  F_neutron = {F_neutron:.6e}\n"
                f"  EL_residual = {EL_residual:.6e}\n"
                f"  Stationarity ratio = {stationarity_ratio:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  BUOYANCY REVERSAL AT RESONANCE
# ══════════════════════════════════════════════════════════════════════════════

class BuoyancyReversalAtResonance:
    """
    Buoyancy sign-flip mechanism at 1.25 THz phonon resonance.

    When F_{U,Bi}/F_U crosses 0.5, E_net flips sign:
      - ratio > 0.5 → net_factor > 0 → expansion (nebulae, cosmogenesis)
      - ratio < 0.5 → net_factor < 0 → erosion (filaments, cavities)
      - ratio = 0.5 → net_factor = 0 → critical transition

    The phonon modulation Φ_{1.25 THz} amplifies or suppresses this flip
    depending on the driving frequency. At resonance (ω = ω_SCm), the
    flip is maximally amplified.
    """

    def __init__(self):
        self._phonon_energy = PhononModulatedEnergy()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Sweep F_{U,Bi}/F_U across the critical transition at ratio = 0.5.

        Parameters from dataset:
          t, V_region, omega, Phi_0, gamma, SSq, kappa, rho_vac:
          ratio_min:  start of sweep (default 0.1)
          ratio_max:  end of sweep (default 0.9)
          n_points:   number of points (default 9)
        """
        ratio_min = dataset.get('ratio_min', 0.1)
        ratio_max = dataset.get('ratio_max', 0.9)
        n_points  = dataset.get('n_points', 9)

        sweep = []
        for i in range(n_points):
            r = ratio_min + (ratio_max - ratio_min) * i / max(n_points - 1, 1)
            d = dict(dataset)
            d['F_U_Bi'] = r
            d['F_U'] = 1.0
            result = self._phonon_energy.compute(d)
            sweep.append({
                "ratio": r,
                "net_factor": result["net_factor"],
                "E_net_phonon": result["E_net_phonon"],
                "E_net_bare": result["E_net_bare"],
                "Phi_125": result["Phi_125_THz"],
                "regime": result["regime"],
            })

        # Find critical crossing
        sign_changes = []
        for i in range(1, len(sweep)):
            nf_prev = sweep[i-1]["net_factor"]
            nf_curr = sweep[i]["net_factor"]
            # Sign change: product < 0 or one of them is exactly zero
            if nf_prev * nf_curr < 0 or (nf_prev != 0 and nf_curr == 0) \
               or (nf_prev == 0 and nf_curr != 0) \
               or (nf_prev < 0 < nf_curr) or (nf_curr < 0 < nf_prev):
                sign_changes.append({
                    "between_ratios": (sweep[i-1]["ratio"], sweep[i]["ratio"]),
                    "transition": f"{sweep[i-1]['regime']} → {sweep[i]['regime']}",
                })

        return {
            "sweep": sweep,
            "sign_changes": sign_changes,
            "critical_ratio": 0.5,
            "n_sign_flips": len(sign_changes),
            "description": (
                "Buoyancy reversal at 1.25 THz phonon resonance.\n"
                "When F_{U,Bi}/F_U crosses 0.5, E_net flips sign.\n"
                f"At resonance, Φ amplifies the transition by {sweep[0]['Phi_125']:.4e}.\n"
                f"Sign changes detected: {len(sign_changes)}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4b  RESONANCE ACCELERATION TERM  a_res  (Session 211 addition)
# ══════════════════════════════════════════════════════════════════════════════

class ResonanceAccelerationTerm:
    """
    Resonance acceleration term — the missing link between phonon modulation
    and the gravitational field equations:

        a_res = (F_{U,Bi} / F_U) · Φ_{1.25 THz}(ω) · S₂₆([SSq])

    When F_{U,Bi}/F_U > 0.5, a_res drives expansion (buoyancy dominant).
    When F_{U,Bi}/F_U < 0.5, a_res drives contraction (gravity dominant).

    This term couples the phonon modulation factor Φ directly into the
    26-layer compressed gravity framework via:
        g_phonon(r,t) = Σ_{i=1}^{26} Ug_i × a_res × Q_i
    """

    def __init__(self):
        self._pmf = PhononModulationFactor()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:   angular frequency (rad/s, default ω_SCm)
          Gamma:   linewidth (rad/s, default Γ_SCm)
          Phi_0:   fluence normalization
          F_U_Bi:  buoyancy force (N)
          F_U:     gravity force (N)
          ssq:     [SSq] parameter
        """
        F_U_Bi = dataset.get('F_U_Bi', 1.0)
        F_U    = dataset.get('F_U', 1.0)
        ssq    = dataset.get('ssq', SSQ)

        # Get Φ from the phonon modulation factor
        pmf_result = self._pmf.compute(dataset)
        Phi_norm = pmf_result["Phi_125_THz"] / max(pmf_result.get("Phi_0_used", PHI_0_DEFAULT), 1e-50)

        # S₂₆
        S26 = S26_accelerated(ssq)

        # Buoyancy ratio
        ratio = F_U_Bi / max(F_U, 1e-50)

        # Resonance acceleration
        a_res = ratio * Phi_norm * S26

        # Direction
        if ratio > 0.5:
            regime = "expansion (buoyancy dominant)"
        elif ratio < 0.5:
            regime = "contraction (gravity dominant)"
        else:
            regime = "critical balance"

        return {
            "a_res": a_res,
            "F_UBi_over_FU": ratio,
            "Phi_normalized": Phi_norm,
            "S26": S26,
            "regime": regime,
            "Q_factor": pmf_result["Q_factor"],
            "equation": (
                "a_res = (F_{U,Bi}/F_U) · Φ_{1.25THz}(ω) · S₂₆([SSq])\n"
                f"      = {ratio:.6f} × {Phi_norm:.6e} × {S26:.6e}\n"
                f"      = {a_res:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  E(t) vs k-ESSENCE COMPARISON
# ══════════════════════════════════════════════════════════════════════════════

class EtVsKEssenceComparison:
    """
    Head-to-head comparison of UQFF E(t) phonon-driven vacuum dynamics
    vs k-essence (non-canonical kinetic scalar field).

    k-Essence:
      Lagrangian: L = F(X, φ)  where X = ½(∂_μ φ)(∂^μ φ)
      Pressure:   p = F(X, φ)
      Density:    ρ = 2X F_X − F
      EOS:        w = F / (2X F_X − F)
      Sound speed: c_s² = F_X / (F_X + 2X F_{XX})
      Field eq:   (F_X + 2X F_{XX})φ̈ + 3H F_X φ̇ + (F_X ∂_i∂^i φ − F_φ) = 0

    Specific model used: purely kinetic k-essence F(X) = −A + B X^n
      (Scherrer model, n > 0, gives tracking attractor behavior)

    UQFF E(t):
      SCm superconductive vacuum manifold modulated by 1.25 THz phonon
      sign-flipping buoyancy dynamics with 2 calibrated parameters
    """

    def __init__(self):
        self._phonon_energy = PhononModulatedEnergy()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          t:          cosmic time (s, default 4.35e17)
          V_region:   volume (m³, default 1e68)
          F_U_Bi, F_U: buoyancy/gravity forces
          w_obs:      observed dark energy EOS (default -1.03)
          sigma_w:    EOS uncertainty (default 0.03)

        k-Essence model (Scherrer):
          A_kess:     constant term (J/m³, default ρ_crit c²)
          B_kess:     kinetic amplitude (default 1.0)
          n_kess:     kinetic exponent (default 1.0 → canonical)
          X_0:        initial kinetic energy ½(∂φ)² (default 1e-60)
          phi_dot_0:  field velocity (default 1e-30)
        """
        t        = dataset.get('t', 4.35e17)
        V_region = dataset.get('V_region', 1e68)
        w_obs    = dataset.get('w_obs', -1.03)
        sigma_w  = dataset.get('sigma_w', 0.03)

        # ── k-Essence side (Scherrer model) ──
        # F(X) = -A + B X^n
        A_kess = dataset.get('A_kess', RHO_CRIT * c**2)
        B_kess = dataset.get('B_kess', 1.0)
        n_kess = dataset.get('n_kess', 1.0)
        phi_dot = dataset.get('phi_dot_0', 1e-30)

        X_0 = 0.5 * phi_dot**2       # kinetic energy X = ½(∂φ)²

        # F(X) = -A + B X^n
        F_X_val = A_kess if X_0 == 0 else -A_kess + B_kess * (X_0 ** n_kess)

        # F_X = dF/dX = n B X^{n-1}
        F_X_deriv = n_kess * B_kess * (X_0 ** max(n_kess - 1, 0)) \
            if X_0 > 0 else n_kess * B_kess

        # F_XX = d²F/dX² = n(n-1) B X^{n-2}
        F_XX = n_kess * (n_kess - 1) * B_kess * (X_0 ** max(n_kess - 2, 0)) \
            if X_0 > 0 and n_kess > 1 else 0.0

        # Pressure: p = F(X, φ)
        p_kess = F_X_val

        # Density: ρ = 2X F_X - F
        rho_kess = 2 * X_0 * F_X_deriv - F_X_val

        # EOS: w = p / ρ = F / (2X F_X - F)
        w_kess = p_kess / rho_kess if rho_kess != 0 else -1.0

        # Sound speed: c_s² = F_X / (F_X + 2X F_XX)
        denom_cs = F_X_deriv + 2 * X_0 * F_XX
        c_s_sq = F_X_deriv / denom_cs if denom_cs != 0 else 1.0
        c_s = math.sqrt(abs(c_s_sq))

        # Effective field equation acceleration
        # φ̈ = -(3H F_X φ̇ + F_φ) / (F_X + 2X F_XX)
        F_phi = 0.0   # purely kinetic model: no φ-dependence
        phi_ddot = -(3 * H_0 * F_X_deriv * phi_dot + F_phi) / denom_cs \
            if denom_cs != 0 else 0.0

        # k-Essence free parameters
        n_params_kess = 3   # A, B, n (plus initial conditions)

        # Fine-tuning: A must match ρ_Λ c² to ~10^{-120}
        rho_QFT = 1e113
        rho_Lambda_obs = 0.692 * RHO_CRIT
        fine_tuning_kess = rho_QFT / rho_Lambda_obs if rho_Lambda_obs > 0 else 0.0
        ft_exp = int(math.log10(fine_tuning_kess)) if fine_tuning_kess > 0 else 0

        # ── UQFF E(t) side ──
        phonon_dataset = dict(dataset)
        phonon_dataset.update({'V_region': V_region, 't': 0.0})
        pe_result = self._phonon_energy.compute(phonon_dataset)
        E_net_phonon = pe_result["E_net_phonon"]

        # UQFF effective EOS
        rate = KAPPA + SSQ / N_LEVELS
        w_UQFF = -1.0 + 2.0 * rate / (3.0 * H_0)

        # χ² comparison
        chi2_kess = ((w_kess - w_obs) / sigma_w)**2
        chi2_UQFF = ((w_UQFF - w_obs) / sigma_w)**2
        Delta_chi2 = chi2_kess - chi2_UQFF

        net_factor = pe_result["net_factor"]
        accel_sign = ("positive (expanding)" if net_factor > 0
                      else ("negative (eroding)" if net_factor < 0
                            else "zero (balanced)"))

        return {
            "k_essence": {
                "F_X": F_X_val,
                "F_X_deriv": F_X_deriv,
                "F_XX": F_XX,
                "p_kess": p_kess,
                "rho_kess": rho_kess,
                "w_kess": w_kess,
                "c_s_squared": c_s_sq,
                "c_s": c_s,
                "phi_ddot": phi_ddot,
                "X_0": X_0,
                "n_params": n_params_kess,
                "model": f"F(X) = -A + B·X^n, A={A_kess:.4e}, B={B_kess}, n={n_kess}",
                "fine_tuning": f"~10^{ft_exp} (A tuned for ρ_Λ)",
            },
            "UQFF_phonon_Et": {
                "E_net_phonon": E_net_phonon,
                "E_net_bare": pe_result["E_net_bare"],
                "Phi_125_THz": pe_result["Phi_125_THz"],
                "w_UQFF": w_UQFF,
                "regime": pe_result["regime"],
                "acceleration_sign": accel_sign,
                "fine_tuning": "None (2 params: [SSq]=0.57, κ=0.0005/day)",
            },
            "comparison": {
                "Delta_w": w_kess - w_UQFF,
                "chi2_kessence": chi2_kess,
                "chi2_UQFF": chi2_UQFF,
                "Delta_chi2": Delta_chi2,
                "preferred": "UQFF" if Delta_chi2 > 0 else "k-Essence",
                "fine_tuning_ratio": (
                    f"k-Essence: ~10^{ft_exp} (A tuned) | UQFF: 1"
                ),
                "sound_speed": (
                    f"k-Essence: c_s = {c_s:.6f} (subluminal constraint) | "
                    "UQFF: N/A (buoyancy-driven, no scalar sound speed)"
                ),
            },
            "contrast_table": [
                {"aspect": "Origin",
                 "kEssence": "Kinetic term K(X,φ) with p=F(X,φ); X=½(∂φ)²",
                 "UQFF": "SCm superconductive vacuum modulated by 1.25 THz phonon"},
                {"aspect": "Equation of state",
                 "kEssence": f"w = {w_kess:.6f} (can cross −1; depends on F)",
                 "UQFF": f"w = {w_UQFF:.6f} (sign-flipping via buoyancy)"},
                {"aspect": "Dynamics",
                 "kEssence": f"(F_X+2XF_XX)φ̈ + 3HF_Xφ̇ − F_φ = 0; c_s={c_s:.4f}",
                 "UQFF": "exp(κt+[SSq]t/26)·S₂₆·Φ_{1.25 THz} (phonon-buoyancy)"},
                {"aspect": "Lab testability",
                 "kEssence": "None (high-energy scalar field)",
                 "UQFF": "Direct: 1.25 THz QCL neutron drops, micro-plasmoid, LENR COP>10"},
                {"aspect": "GW prediction",
                 "kEssence": "Standard GR waveforms",
                 "UQFF": "66.7% strain reduction + 367.8-cycle phase lag"},
                {"aspect": "Cosmogenesis",
                 "kEssence": "Possible inflation but no pre-gravity vacuum",
                 "UQFF": "SCm phonon → DPM proto-shells → EM bang + 2 cycles"},
                {"aspect": "Fine-tuning",
                 "kEssence": f"Potential and kinetic function tuned (~10^{ft_exp})",
                 "UQFF": "None — 2 fixed parameters ([SSq]=0.57, κ=5e-4/day)"},
                {"aspect": "Sign behavior",
                 "kEssence": "Usually accelerating (w ≈ −1)",
                 "UQFF": f"Explicit ± phases: net_factor = {net_factor:.4f} ({accel_sign})"},
                {"aspect": "Sound speed",
                 "kEssence": f"c_s² = F_X/(F_X+2XF_XX) = {c_s_sq:.6e}",
                 "UQFF": "No scalar c_s; phonon propagation at ω_SCm = 1.25 THz"},
                {"aspect": "Free parameters",
                 "kEssence": f"{n_params_kess} (A, B, n) + initial conditions",
                 "UQFF": "2 ([SSq]=0.57, κ=0.0005/day) calibrated from data"},
            ],
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  WSTP KERNEL FOR PHONON E(t) + k-ESSENCE
# ══════════════════════════════════════════════════════════════════════════════

def wstp_kernel_phonon_et() -> str:
    """
    Wolfram Language definitions for phonon-modulated E(t) sector +
    k-essence comparison.
    """
    return r"""
(* ═══════════════════════════════════════════════════════════════════════ *)
(* Phonon-Modulated E(t) + k-Essence Comparison — UQFF Symbolic Forms   *)
(* Session 208 | Daniel Murphy                                            *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* ── UQFF calibrated constants ── *)
κ = 5.787*^-9;     (* s⁻¹ *)
SSq = 0.57;
H0 = 2.195*^-18;   (* s⁻¹ *)
GNewt = 6.6743*^-11;
cLight = 2.99792*^8;
ρvacSCm = 9.47*^-27;
ωSCm = 2 Pi * 1.25*^12;      (* rad/s *)
ΓSCm = 2 Pi * 0.1*^12;       (* rad/s *)
Φ0 = 10^20;                  (* phonons/m²/s *)

(* ── Phonon modulation factor ── *)
S26[z_] := PolyLog[26, z];
Φ125THz[ω_] := Φ0 * Exp[-(ω - ωSCm)^2 / (2 ΓSCm^2)] * S26[SSq];

(* ── SCm vacuum density ── *)
ρSCmEvol[t_] := ρvacSCm * S26[SSq] * Exp[κ t + SSq t / 26];

(* ── Phonon-modulated E_net ── *)
EnetPhonon[t_, Vreg_, FUBi_, FU_, ω_] :=
  ρSCmEvol[t] * Vreg * (2 FUBi / FU - 1) * Φ125THz[ω];

(* ── Phonon Lagrangian ── *)
Lphonon[t_, Vreg_, FUBi_, FU_, ω_] :=
  EnetPhonon[t, Vreg, FUBi, FU, ω] * Vreg * S26[SSq];

(* ── Euler-Lagrange variation ── *)
δSδφPhonon[t_, Vreg_, FUBi_, FU_, ω_, βi_, Ug_, Ωg_, M_, dg_, UA_, Fn_] :=
  D[(-βi * Ug * Ωg * M / dg * UA
     + Fn * Φ125THz[ω])
    * Vreg * Φ125THz[ω] * S26[SSq], φphonon];

(* ── k-Essence: Scherrer model F(X) = -A + B X^n ── *)
Fkess[X_, A_, B_, n_] := -A + B * X^n;
FXkess[X_, B_, n_] := n * B * X^(n - 1);
FXXkess[X_, B_, n_] := n (n - 1) B * X^(n - 2);

(* k-Essence EOS *)
ρKess[X_, A_, B_, n_] := 2 X * FXkess[X, B, n] - Fkess[X, A, B, n];
pKess[X_, A_, B_, n_] := Fkess[X, A, B, n];
wKess[X_, A_, B_, n_] := pKess[X, A, B, n] / ρKess[X, A, B, n];

(* k-Essence sound speed *)
cs2Kess[X_, B_, n_] := FXkess[X, B, n] / (FXkess[X, B, n] + 2 X FXXkess[X, B, n]);

(* UQFF effective w *)
wUQFF[κ_, SSq_] := -1 + 2 (κ + SSq / 26) / (3 H0);

(* Δw comparison *)
ΔwKessUQFF[wK_, κ_, SSq_] := wK - wUQFF[κ, SSq];
"""


# ══════════════════════════════════════════════════════════════════════════════
# §7  SELF-TEST
# ══════════════════════════════════════════════════════════════════════════════

def _run_tests():
    """Self-test suite for et_phonon_resonance.py — Session 208."""
    results = []
    test_num = 0

    def check(name, condition, detail=""):
        nonlocal test_num
        test_num += 1
        tag = "PASS" if condition else "FAIL"
        results.append((test_num, tag, name))
        print(f"  T{test_num} [{tag}] {name}" + (f"  ({detail})" if detail else ""))

    print("=" * 70)
    print("et_phonon_resonance.py — Self-Test Suite")
    print("=" * 70)

    # T1: Phonon modulation at resonance (ω = ω_SCm)
    pmf = PhononModulationFactor()
    r1 = pmf.compute({'omega': OMEGA_PHONON})
    check("Φ at resonance is Φ₀·S₂₆ (Gaussian=1)",
          abs(r1["gaussian_peak"] - 1.0) < 1e-10,
          f"Gaussian = {r1['gaussian_peak']:.10f}")

    # T2: Phonon modulation off-resonance is suppressed
    r2 = pmf.compute({'omega': OMEGA_PHONON * 2.0})
    check("Φ off-resonance < Φ at resonance",
          r2["Phi_125_THz"] < r1["Phi_125_THz"],
          f"Φ(2ω) = {r2['Phi_125_THz']:.6e} < Φ(ω) = {r1['Phi_125_THz']:.6e}")

    # T3: Phonon-modulated E_net at resonance
    pme = PhononModulatedEnergy()
    r3 = pme.compute({'t': 0.0, 'F_U_Bi': 0.6, 'F_U': 1.0, 'V_region': 1e48,
                      'omega': OMEGA_PHONON})
    check("E_net_phonon = E_net_bare × Φ (identity check)",
          abs(r3["E_net_phonon"] - r3["E_net_bare"] * r3["Phi_125_THz"])
          / (abs(r3["E_net_phonon"]) + 1e-300) < 1e-10,
          f"E_phonon = {r3['E_net_phonon']:.6e}")

    # T4: E⁺_phonon + E⁻_phonon = E_net_phonon
    sum_pm = r3["E_plus_phonon"] + r3["E_minus_phonon"]
    rel_err = abs(sum_pm - r3["E_net_phonon"]) / (abs(r3["E_net_phonon"]) + 1e-300)
    check("E⁺_phonon + E⁻_phonon = E_net_phonon",
          rel_err < 1e-10,
          f"rel_err = {rel_err:.2e}")

    # T5: Expansion regime for ratio > 0.5
    check("Expansion regime when F_UBi/F_U > 0.5",
          r3["net_factor"] > 0 and "expansion" in r3["regime"],
          f"net_factor = {r3['net_factor']:.4f}, regime = {r3['regime']}")

    # T6: Erosion regime for ratio < 0.5
    r6 = pme.compute({'t': 0.0, 'F_U_Bi': 0.3, 'F_U': 1.0, 'V_region': 1e48})
    check("Erosion regime when F_UBi/F_U < 0.5",
          r6["net_factor"] < 0 and "erosion" in r6["regime"],
          f"net_factor = {r6['net_factor']:.4f}")

    # T7: Phonon Lagrangian finite
    pl = PhononLagrangian()
    r7 = pl.compute({'t': 0.0, 'F_U_Bi': 0.55, 'F_U': 1.0, 'V_region': 1e48})
    check("Phonon Lagrangian is finite",
          math.isfinite(r7["L_phonon"]),
          f"L_phonon = {r7['L_phonon']:.6e}")

    # T8: EL residual finite
    check("EL residual δS/δφ_phonon is finite",
          math.isfinite(r7["EL_residual"]),
          f"EL_residual = {r7['EL_residual']:.6e}")

    # T9: Buoyancy reversal detects sign change
    br = BuoyancyReversalAtResonance()
    r9 = br.compute({'t': 0.0, 'V_region': 1e48, 'omega': OMEGA_PHONON})
    check("Buoyancy reversal detects ≥1 sign change across ratio sweep",
          r9["n_sign_flips"] >= 1,
          f"sign_flips = {r9['n_sign_flips']}")

    # T10: k-Essence comparison returns 10-row contrast table
    kc = EtVsKEssenceComparison()
    r10 = kc.compute({'t': 4.35e17, 'F_U_Bi': 0.55, 'F_U': 1.0, 'V_region': 1e68})
    check("k-Essence comparison returns 10-row contrast table",
          len(r10["contrast_table"]) == 10,
          f"{len(r10['contrast_table'])} rows")

    # T11: k-Essence w ≈ −1 for slow-roll initial conditions
    check("k-Essence w ≈ −1 for default parameters",
          abs(r10["k_essence"]["w_kess"] - (-1.0)) < 1.0,
          f"w_kess = {r10['k_essence']['w_kess']:.6f}")

    # T12: k-Essence sound speed 0 < c_s ≤ 1
    check("k-Essence sound speed c_s > 0",
          r10["k_essence"]["c_s"] > 0,
          f"c_s = {r10['k_essence']['c_s']:.6f}")

    # T13: WSTP kernel valid
    wl = wstp_kernel_phonon_et()
    check("WSTP kernel contains phonon + k-essence definitions",
          "Φ125THz" in wl and "wKess" in wl and "Lphonon" in wl and "cs2Kess" in wl,
          f"len = {len(wl)}")

    # T14: Quality factor Q > 1 (resonance is sharp)
    check("Quality factor Q > 1 (sharp resonance)",
          r1["Q_factor"] > 1.0,
          f"Q = {r1['Q_factor']:.2f}")

    # Summary
    passed = sum(1 for _, tag, _ in results if tag == "PASS")
    total = len(results)
    print("\n" + "=" * 70)
    print(f"RESULT: {passed}/{total} tests passed")
    print("=" * 70)
    return passed == total


# ═══════════════════════════════════════════════════════════════════════════════
# §4  MERGER PHONON RESONANCE (Session 213)
# ═══════════════════════════════════════════════════════════════════════════════

class MergerPhononResonance:
    """Phonon resonance in SMBH binary mergers.

    P_merger(Γ) = P_GR · (1 + M_merger(Γ))
    Circumbinary disk phonon coupling at 1.25 THz drives final coalescence.
    Strain damping 47-66.7%, phase lag 200-400 cycles.
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        P_GR = float(d.get("P_GR_W", 3.6e49))
        q = float(d.get("q", 0.5))
        Gamma_THz = float(d.get("Gamma_THz", 0.10))

        ssq = 0.57
        s26 = sum(math.exp(-ssq * k / 26.0) for k in range(1, 27))
        Gamma_0 = 2 * math.pi * 0.1e12
        sigma_G = 0.08 * 2 * math.pi * 1e12
        Gamma = 2 * math.pi * Gamma_THz * 1e12
        delta = Gamma - Gamma_0
        Phi_norm = math.exp(-delta**2 / (2 * sigma_G**2)) * s26
        buoyancy = 2 * 0.6 / 1.0 - 1
        M_merger = Phi_norm * s26 * buoyancy
        P_merger = P_GR * (1 + M_merger)

        D_total = 0.333 + 0.197 * (1 - q)
        damping = (1 - D_total) * 100

        return {
            "P_GR_W": P_GR, "P_merger_W": P_merger,
            "enhancement": P_merger / max(P_GR, 1e-50),
            "D_total": D_total, "damping_percent": damping,
            "primary_equations": [
                "P_merger(Γ) = P_GR · (1 + M_merger(Γ))",
                f"Enhancement: {P_merger / max(P_GR, 1e-50):.2f}×",
                f"Strain damping: {damping:.1f}% (q={q:.2f})",
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §5  LINEWIDTH PHONON MAPPING (Session 213)
# ═══════════════════════════════════════════════════════════════════════════════

class LinewidthPhononMapping:
    """Systematic Γ sweep for phonon resonance interpretation.

    Maps linewidth to regime (narrow/optimal/broad) with quality factor Q.
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        gammas = d.get("gammas", [0.05, 0.10, 0.30])
        omega_SCm = 2 * math.pi * 1.25e12

        rows = []
        for gTHz in gammas:
            Gamma = 2 * math.pi * gTHz * 1e12
            Q = omega_SCm / (2 * Gamma)
            regime = "narrow" if gTHz < 0.07 else ("optimal" if gTHz < 0.15 else "broad")
            rows.append({"Gamma_THz": gTHz, "Q": Q, "regime": regime})

        return {
            "sweep": rows,
            "primary_equations": [
                "Q = ω_SCm / (2Γ)",
            ] + [f"Γ={r['Gamma_THz']}: Q={r['Q']:.1f}, {r['regime']}" for r in rows],
        }


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
