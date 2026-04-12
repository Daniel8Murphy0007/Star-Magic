#!/usr/bin/env python3
"""
bh_phonon_interaction.py — Black Hole Phonon Interaction Engine

Session 211 | Daniel Murphy
PURPOSE: Standalone physics engine for black hole-phonon interaction:
         ergosphere superradiance coupling, phonon-modified Hawking temperature,
         and QPO accretion disk phonon coupling.

         Key physics:
           - Superradiance: Γ_SR = Φ · (ω - m·Ω_H) · α_BH when ω < m·Ω_H
           - T_Hawking^phonon = T_H · (1 + Φ · S₂₆ · [SSq] / N_levels)
           - QPO frequency shift: f_QPO^phonon = f_QPO · (1 + M_jet(Γ) · Φ)
           - Ergosphere coupling: E_ergo = ℏ·ω · Φ · (r_ergo/r_H - 1)

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
"""

import math
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

from positive_et_expansion import (
    S26_accelerated,
    G, c, hbar, k_B, mu_0, M_sun, PI,
    KAPPA, SSQ, H_SCM, BETA_I, U_UA,
    RHO_SCM, RHO_UA, RHO_VAC_SCM, N_LEVELS,
    OMEGA_SCM, GAMMA_DEFAULT, SIGMA_0,
)

# ══════════════════════════════════════════════════════════════════════════════
# §0  BH-PHONON CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PHI_0_DEFAULT = 1e20        # phonons/m²/s (THz phonon fluence)
OMEGA_PHONON  = OMEGA_SCM   # 2π × 1.25e12 rad/s
GAMMA_PHONON  = GAMMA_DEFAULT


# ══════════════════════════════════════════════════════════════════════════════
# §1  PHONON-ERGOSPHERE SUPERRADIANCE
# ══════════════════════════════════════════════════════════════════════════════

class PhononErgosphereSuperradiance:
    """
    Phonon-modulated superradiance in a Kerr black hole ergosphere.

    Superradiance condition: ω < m·Ω_H where Ω_H is the horizon angular velocity.

    The phonon field modifies the superradiance amplification factor:
        Γ_SR = Φ_{1.25THz}(ω) · (m·Ω_H − ω) · α_BH

    where α_BH = r_S²/(r_ergo·r_H) is a geometric factor and
    r_ergo = r_S (equatorial), r_H = r_S/2·(1 + √(1−a²)) for spin a.

    Ergosphere phonon energy extraction:
        E_ergo = ℏ·ω_SCm · Φ · (r_ergo/r_H − 1) · S₂₆

    This is the phonon analog of the Penrose process.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:          BH mass (kg, default 10 M_sun)
          a_spin:     dimensionless spin parameter (default 0.9)
          m_mode:     azimuthal mode number (default 1)
          omega:      phonon frequency (rad/s, default ω_SCm)
          Gamma:      linewidth (rad/s, default Γ_SCm)
          Phi_0:      fluence normalization
          ssq:        [SSq] parameter
        """
        M       = dataset.get('M', 10.0 * M_sun)
        a_spin  = dataset.get('a_spin', 0.9)
        m_mode  = dataset.get('m_mode', 1)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0   = dataset.get('Phi_0', PHI_0_DEFAULT)
        ssq     = dataset.get('ssq', SSQ)

        # Schwarzschild radius
        r_S = 2.0 * G * M / c**2

        # Kerr horizon radius
        r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - a_spin**2, 0.0)))

        # Ergosphere radius (equatorial)
        r_ergo = r_S

        # Horizon angular velocity
        a_phys = a_spin * G * M / c  # physical spin parameter
        Omega_H = a_spin * c / (2.0 * r_H) if r_H > 0 else 0.0

        # S₂₆ and Φ
        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi = Phi_0 * gaussian * S26

        # Superradiance condition: ω < m·Ω_H
        omega_SR_max = m_mode * Omega_H
        is_superradiant = omega < omega_SR_max

        # Geometric factor
        alpha_BH = r_S**2 / (r_ergo * max(r_H, 1e-50))

        # Superradiance rate
        if is_superradiant:
            Gamma_SR = Phi * (omega_SR_max - omega) * alpha_BH
        else:
            Gamma_SR = 0.0

        # Ergosphere phonon energy extraction (Penrose analog)
        E_ergo = hbar * OMEGA_PHONON * Phi * (r_ergo / max(r_H, 1e-50) - 1.0) * S26

        # Amplification factor (per phonon cycle)
        Z_lm = Gamma_SR / max(omega, 1e-50) if is_superradiant else 0.0

        return {
            "Gamma_SR": Gamma_SR,
            "E_ergo": E_ergo,
            "Z_lm": Z_lm,
            "is_superradiant": is_superradiant,
            "omega_SR_max": omega_SR_max,
            "Omega_H": Omega_H,
            "r_S": r_S,
            "r_H": r_H,
            "r_ergo": r_ergo,
            "alpha_BH": alpha_BH,
            "Phi": Phi,
            "S26": S26,
            "a_spin": a_spin,
            "equation": (
                "Γ_SR = Φ · (m·Ω_H − ω) · α_BH\n"
                f"     = {Phi:.4e} × ({omega_SR_max:.4e} − {omega:.4e}) × {alpha_BH:.4e}\n"
                f"     = {Gamma_SR:.6e} rad/s\n\n"
                f"E_ergo = ℏω · Φ · (r_ergo/r_H − 1) · S₂₆ = {E_ergo:.6e} J"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  PHONON-MODIFIED HAWKING TEMPERATURE
# ══════════════════════════════════════════════════════════════════════════════

class PhononModifiedHawkingTemperature:
    """
    Phonon correction to the Hawking temperature of a black hole:

        T_H^phonon = T_H · (1 + Φ_{1.25THz} · S₂₆ · [SSq] / N_levels)

    where T_H = ℏc³ / (8π G M k_B) is the standard Hawking temperature.

    The phonon correction arises from the SCm vacuum manifold modifying
    the near-horizon modes. At 1.25 THz resonance, virtual phonon pairs
    in the vacuum contribute an additional spectral weight to the Hawking
    radiation.

    Also computes:
      - Modified BH lifetime: τ_BH^phonon = τ_BH / (1 + correction)³
      - Modified luminosity: L_H^phonon = L_H · (1 + correction)⁴
      - Page time modification
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:        BH mass (kg, default 10 M_sun)
          omega:    probe frequency (rad/s, default ω_SCm)
          Gamma:    linewidth (rad/s, default Γ_SCm)
          Phi_0:    fluence normalization
          ssq, n_levels: UQFF parameters
        """
        M       = dataset.get('M', 10.0 * M_sun)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0   = dataset.get('Phi_0', PHI_0_DEFAULT)
        ssq     = dataset.get('ssq', SSQ)
        n_levels = dataset.get('n_levels', N_LEVELS)

        S26 = S26_accelerated(ssq)

        # Standard Hawking temperature
        T_H = hbar * c**3 / (8.0 * PI * G * M * k_B)

        # Phonon modulation at resonance
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26  # normalized

        # Correction factor
        correction = Phi_norm * S26 * ssq / n_levels

        # Modified temperature
        T_H_phonon = T_H * (1.0 + correction)

        # Hawking luminosity: L_H = ℏc⁶ / (15360 π G² M²)
        L_H = hbar * c**6 / (15360.0 * PI * G**2 * M**2)
        L_H_phonon = L_H * (1.0 + correction)**4

        # BH evaporation time: τ = 5120 π G² M³ / (ℏ c⁴)
        tau_BH = 5120.0 * PI * G**2 * M**3 / (hbar * c**4)
        tau_BH_phonon = tau_BH / (1.0 + correction)**3

        # Page time (information release halfway): t_Page ≈ τ / 2
        t_Page = tau_BH_phonon / 2.0

        return {
            "T_H": T_H,
            "T_H_phonon": T_H_phonon,
            "correction": correction,
            "L_H": L_H,
            "L_H_phonon": L_H_phonon,
            "tau_BH": tau_BH,
            "tau_BH_phonon": tau_BH_phonon,
            "t_Page_phonon": t_Page,
            "Phi_norm": Phi_norm,
            "S26": S26,
            "equation": (
                "T_H^phonon = T_H · (1 + Φ · S₂₆ · [SSq] / N)\n"
                f"T_H     = {T_H:.6e} K\n"
                f"T_H^ph  = {T_H_phonon:.6e} K\n"
                f"ΔT/T    = {correction:.6e}\n"
                f"τ_BH    = {tau_BH:.6e} s → τ_BH^ph = {tau_BH_phonon:.6e} s"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  QPO ACCRETION DISK PHONON COUPLING
# ══════════════════════════════════════════════════════════════════════════════

class QPOAccretionDiskPhononCoupling:
    """
    Phonon coupling to quasi-periodic oscillations (QPOs) in
    accretion disks around black holes.

    The QPO frequency is shifted by the phonon modulation:
        f_QPO^phonon = f_QPO · (1 + M_jet(Γ) · Φ_{1.25THz})

    where M_jet(Γ) = 1 + A_jet · exp(-(Γ − Γ₀)²/(2σ_Γ²)) is the
    jet modulation factor defined in the jet phonon module.

    Also computes the phonon-modified ISCO frequency:
        f_ISCO^phonon = f_ISCO · (1 + Φ · S₂₆) where f_ISCO = c³/(6√6 π G M)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:         BH mass (kg)
          a_spin:    dimensionless spin
          f_QPO_Hz:  observed QPO frequency (Hz, default 300)
          A_jet:     jet modulation amplitude (default 1.5)
          sigma_Gamma_THz: jet modulation width in THz (default 0.08)
          omega, Gamma, Phi_0, ssq: phonon parameters
        """
        M       = dataset.get('M', 10.0 * M_sun)
        a_spin  = dataset.get('a_spin', 0.5)
        f_QPO   = dataset.get('f_QPO_Hz', 300.0)
        A_jet   = dataset.get('A_jet', 1.5)
        sigma_G = dataset.get('sigma_Gamma_THz', 0.08) * 2.0 * PI * 1e12
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0   = dataset.get('Phi_0', PHI_0_DEFAULT)
        ssq     = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)

        # Phonon modulation
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        # Jet modulation factor M_jet(Γ)
        Gamma_0 = GAMMA_PHONON
        M_jet = 1.0 + A_jet * math.exp(-(Gamma - Gamma_0)**2 / (2.0 * sigma_G**2))

        # Modified QPO frequency
        f_QPO_phonon = f_QPO * (1.0 + M_jet * Phi_norm)

        # ISCO frequency
        # For Schwarzschild: r_ISCO = 6GM/c², f_ISCO = c³/(6√6 π G M)
        f_ISCO = c**3 / (6.0 * math.sqrt(6.0) * PI * G * M)
        f_ISCO_phonon = f_ISCO * (1.0 + Phi_norm * S26)

        # QPO-to-ISCO ratio
        ratio = f_QPO / max(f_ISCO, 1e-50)
        ratio_phonon = f_QPO_phonon / max(f_ISCO_phonon, 1e-50)

        return {
            "f_QPO_Hz": f_QPO,
            "f_QPO_phonon_Hz": f_QPO_phonon,
            "Delta_f_QPO_Hz": f_QPO_phonon - f_QPO,
            "M_jet_Gamma": M_jet,
            "f_ISCO_Hz": f_ISCO,
            "f_ISCO_phonon_Hz": f_ISCO_phonon,
            "QPO_ISCO_ratio": ratio,
            "QPO_ISCO_ratio_phonon": ratio_phonon,
            "Phi_norm": Phi_norm,
            "S26": S26,
            "equation": (
                "f_QPO^phonon = f_QPO · (1 + M_jet(Γ) · Φ)\n"
                f"  = {f_QPO:.2f} × (1 + {M_jet:.4f} × {Phi_norm:.6e})\n"
                f"  = {f_QPO_phonon:.6f} Hz\n"
                f"f_ISCO^phonon = {f_ISCO_phonon:.6f} Hz"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  PHONON-MODIFIED BH ENTROPY
# ══════════════════════════════════════════════════════════════════════════════

class PhononModifiedBHEntropy:
    """
    Phonon correction to the Bekenstein-Hawking entropy:

        S_BH^phonon = S_BH · (1 + Φ · S₂₆ · [SSq] / N_levels)²

    where S_BH = A / (4 l_P²) = 4π G² M² / (ℏ c³).

    The phonon correction preserves the area law but modifies the
    effective Planck area via SCm vacuum corrections.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        M       = dataset.get('M', 10.0 * M_sun)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        ssq     = dataset.get('ssq', SSQ)
        n_levels = dataset.get('n_levels', N_LEVELS)

        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        correction = Phi_norm * S26 * ssq / n_levels

        # Planck length
        l_P = math.sqrt(hbar * G / c**3)

        # Schwarzschild radius and horizon area
        r_S = 2.0 * G * M / c**2
        A_horizon = 4.0 * PI * r_S**2

        # Bekenstein-Hawking entropy
        S_BH = A_horizon / (4.0 * l_P**2)
        S_BH_phonon = S_BH * (1.0 + correction)**2

        # Information capacity (bits)
        I_BH = S_BH / (k_B * math.log(2.0))
        I_BH_phonon = S_BH_phonon / (k_B * math.log(2.0))

        return {
            "S_BH": S_BH,
            "S_BH_phonon": S_BH_phonon,
            "Delta_S": S_BH_phonon - S_BH,
            "correction": correction,
            "I_BH_bits": I_BH,
            "I_BH_phonon_bits": I_BH_phonon,
            "A_horizon_m2": A_horizon,
            "r_S_m": r_S,
            "equation": (
                "S_BH^phonon = S_BH · (1 + Φ·S₂₆·[SSq]/N)²\n"
                f"S_BH     = {S_BH:.6e} J/K\n"
                f"S_BH^ph  = {S_BH_phonon:.6e} J/K\n"
                f"ΔS/S     = {(S_BH_phonon/S_BH - 1.0):.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  BLAZAR ERGOSPHERE PHONON COUPLING  (Session 212)
# ══════════════════════════════════════════════════════════════════════════════

class BlazarErgosphereCoupling:
    """
    Blazar-specific ergosphere phonon coupling for BL Lac / Mrk 421 type systems.

    Extends PhononErgosphereSuperradiance with relativistic Doppler boosting:

        E_ergo^blazar = (a/2) M c² · Φ · S₂₆ · δ_D
        δ_D = 1 / [Γ_bulk (1 − β cos θ)]     (Doppler factor)

    Superradiant condition in blazar frame: ω_obs = ω_SCm · δ_D < Ω_H.
    Jet-ergosphere coupling: P_ergo = E_ergo · Ω_H / (2π).
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:          BH mass (kg, default 3e8 M☉ — BL Lac type)
          a_spin:     spin (default 0.95)
          Gamma_bulk: bulk Lorentz factor (default 15)
          theta_obs:  observer angle (rad, default 0.05)
          omega:      phonon frequency (rad/s)
          Gamma:      phonon linewidth (rad/s)
          ssq:        [SSq]
        """
        M           = dataset.get('M', 3e8 * M_sun)
        a_spin      = dataset.get('a_spin', 0.95)
        Gamma_bulk  = dataset.get('Gamma_bulk', 15.0)
        theta_obs   = dataset.get('theta_obs', 0.05)
        omega       = dataset.get('omega', OMEGA_PHONON)
        Gamma       = dataset.get('Gamma', GAMMA_PHONON)
        ssq         = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)

        # Doppler factor
        beta = math.sqrt(1.0 - 1.0 / max(Gamma_bulk**2, 1.0 + 1e-50))
        delta_D = 1.0 / (Gamma_bulk * (1.0 - beta * math.cos(theta_obs)))

        # Horizon
        r_S = 2.0 * G * M / c**2
        r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - a_spin**2, 0.0)))
        Omega_H = a_spin * c / (2.0 * r_H)

        # Doppler-shifted phonon
        omega_obs = omega * delta_D
        is_superradiant = omega_obs < Omega_H

        # Phonon modulation
        delta_omega = omega_obs - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        # Ergosphere energy extraction
        E_ergo = (a_spin / 2.0) * M * c**2 * Phi_norm * S26 * delta_D
        P_ergo = E_ergo * Omega_H / (2.0 * PI)

        return {
            "delta_D": delta_D,
            "omega_obs_rad_s": omega_obs,
            "Omega_H_rad_s": Omega_H,
            "is_superradiant": is_superradiant,
            "E_ergo_J": E_ergo,
            "P_ergo_W": P_ergo,
            "Phi_norm": Phi_norm,
            "r_H_m": r_H,
            "equation": (
                "Blazar ergosphere phonon coupling:\n"
                f"  δ_D = {delta_D:.4f}\n"
                f"  ω_obs = {omega_obs:.6e} rad/s\n"
                f"  Ω_H = {Omega_H:.6e} rad/s\n"
                f"  Superradiant: {is_superradiant}\n"
                f"  E_ergo = {E_ergo:.6e} J\n"
                f"  P_ergo = {P_ergo:.6e} W"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate BH phonon interaction engine."""
    print("=" * 72)
    print("Black Hole Phonon Interaction Engine — Session 211+212")
    print("=" * 72)

    # §1 Ergosphere superradiance
    print("\n── §1 Phonon-Ergosphere Superradiance ──")
    sr = PhononErgosphereSuperradiance()
    result = sr.compute({'M': 10.0 * M_sun, 'a_spin': 0.9})
    print(f"  Γ_SR = {result['Gamma_SR']:.6e} rad/s")
    print(f"  E_ergo = {result['E_ergo']:.6e} J")
    print(f"  Superradiant: {result['is_superradiant']}")

    # §2 Hawking temperature
    print("\n── §2 Phonon-Modified Hawking Temperature ──")
    hawking = PhononModifiedHawkingTemperature()
    result = hawking.compute({'M': 10.0 * M_sun})
    print(f"  T_H = {result['T_H']:.6e} K")
    print(f"  T_H^phonon = {result['T_H_phonon']:.6e} K")
    print(f"  Correction: {result['correction']:.6e}")

    # §3 QPO coupling
    print("\n── §3 QPO Accretion Disk Phonon Coupling ──")
    qpo = QPOAccretionDiskPhononCoupling()
    result = qpo.compute({'M': 10.0 * M_sun, 'f_QPO_Hz': 300.0})
    print(f"  f_QPO = {result['f_QPO_Hz']:.2f} Hz → {result['f_QPO_phonon_Hz']:.6f} Hz")
    print(f"  M_jet(Γ) = {result['M_jet_Gamma']:.4f}")

    # §4 BH entropy
    print("\n── §4 Phonon-Modified BH Entropy ──")
    entropy = PhononModifiedBHEntropy()
    result = entropy.compute({'M': 10.0 * M_sun})
    print(f"  S_BH = {result['S_BH']:.6e} J/K")
    print(f"  S_BH^phonon = {result['S_BH_phonon']:.6e} J/K")

    print(f"\n{'=' * 72}")
    print("BH PHONON INTERACTION COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
