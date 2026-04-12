#!/usr/bin/env python3
"""
blazar_jet_phonon.py — Blazar Jet Phonon Coupling Engine

Session 212 | Daniel Murphy
PURPOSE: Standalone engine for phonon-modulated blazar jet physics,
         extending AGN jet models to relativistic blazar-specific systems.

         Key physics:
           - BlazarErgosphereResonance: BL Lac / Mrk 421 / PKS-type
             ergosphere phonon coupling with Doppler-boosted spectra
           - BlazarJetPowerGammaCurve: Γ-stepped jet power curves
             for blazar jets with beaming correction
           - BlazarMultiMessengerPhononCorrelation: VHE γ-ray and
             neutrino correlation via phonon-enhanced pair cascades

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
    RHO_SCM, RHO_UA, N_LEVELS,
    OMEGA_SCM, GAMMA_DEFAULT,
)


# ══════════════════════════════════════════════════════════════════════════════
# §0  BLAZAR-SPECIFIC CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PHI_0_DEFAULT    = 1e20
OMEGA_PHONON     = OMEGA_SCM
GAMMA_PHONON     = GAMMA_DEFAULT
A_JET_DEFAULT    = 1.5
SIGMA_GAMMA_DEFAULT = 0.08   # THz

# Gamma stepped linewidths (rad/s)
GAMMA_NARROW  = 2.0 * PI * 0.05e12
GAMMA_OPTIMAL = 2.0 * PI * 0.10e12
GAMMA_BROAD   = 2.0 * PI * 0.30e12
GAMMA_STEPS   = [GAMMA_NARROW, GAMMA_OPTIMAL, GAMMA_BROAD]
GAMMA_LABELS  = ["0.05 THz", "0.10 THz", "0.30 THz"]


def _m_jet(Gamma_rad_s: float, A_jet: float = A_JET_DEFAULT) -> float:
    """Jet modulation factor M_jet(Γ) = 1 + A_jet · exp[-(Γ−Γ₀)²/(2σ_Γ²)]."""
    sigma_G = SIGMA_GAMMA_DEFAULT * 2.0 * PI * 1e12
    delta = Gamma_rad_s - GAMMA_PHONON
    return 1.0 + A_jet * math.exp(-delta**2 / (2.0 * sigma_G**2))


def _p_bz(M: float, a_spin: float, B: float) -> float:
    """Blandford-Znajek power: P_BZ = (B²/8π)(r_H/c)²a²c."""
    r_S = 2.0 * G * M / c**2
    r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - a_spin**2, 0.0)))
    return (B**2 / (8.0 * PI)) * (r_H / c)**2 * a_spin**2 * c


# ══════════════════════════════════════════════════════════════════════════════
# §1  BLAZAR ERGOSPHERE PHONON RESONANCE
# ══════════════════════════════════════════════════════════════════════════════

class BlazarErgosphereResonance:
    """
    Phonon-coupled ergosphere superradiance for blazar-class AGN.

    Blazar-specific: relativistic Doppler boosting shifts the observed
    phonon coupling into the jet frame:

        ω_obs = ω_SCm · δ_D
        δ_D = 1 / [Γ_bulk · (1 − β cos θ)]  (Doppler factor)

    Superradiant condition: ω_obs < m · Ω_H where Ω_H = a c / (2 r_H).
    Energy extraction: E_ergo = (a/2) · M c² · Φ · S₂₆ · δ_D.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:          BH mass (kg, default 1e9 M☉)
          a_spin:     dimensionless spin (default 0.95)
          Gamma_bulk: bulk Lorentz factor (default 10)
          theta_obs:  observer angle (rad, default 0.1 ~ 5.7°)
          B:          magnetic field (T, default 1.0)
          omega:      phonon frequency (rad/s, default ω_SCm)
          Gamma:      phonon linewidth (rad/s, default Γ_SCm)
          ssq:        [SSq] (default 0.57)
        """
        M           = dataset.get('M', 1e9 * M_sun)
        a_spin      = dataset.get('a_spin', 0.95)
        Gamma_bulk  = dataset.get('Gamma_bulk', 10.0)
        theta_obs   = dataset.get('theta_obs', 0.1)
        B           = dataset.get('B', 1.0)
        omega       = dataset.get('omega', OMEGA_PHONON)
        Gamma       = dataset.get('Gamma', GAMMA_PHONON)
        ssq         = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)

        # Doppler factor
        beta = math.sqrt(1.0 - 1.0 / max(Gamma_bulk**2, 1.0 + 1e-50))
        delta_D = 1.0 / (Gamma_bulk * (1.0 - beta * math.cos(theta_obs)))

        # Doppler-shifted phonon frequency
        omega_obs = omega * delta_D

        # Horizon parameters
        r_S = 2.0 * G * M / c**2
        r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - a_spin**2, 0.0)))
        Omega_H = a_spin * c / (2.0 * r_H)

        # Superradiant condition (m=1 mode)
        is_superradiant = omega_obs < Omega_H

        # Phonon modulation at observed frequency
        delta_omega = omega_obs - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        # Ergosphere energy extraction
        E_ergo = (a_spin / 2.0) * M * c**2 * Phi_norm * S26 * delta_D

        # Superradiant amplification rate
        Gamma_SR = Omega_H * a_spin * Phi_norm

        return {
            "delta_D": delta_D,
            "omega_obs_rad_s": omega_obs,
            "Omega_H_rad_s": Omega_H,
            "is_superradiant": is_superradiant,
            "Gamma_SR_rad_s": Gamma_SR,
            "E_ergo_J": E_ergo,
            "Phi_norm": Phi_norm,
            "r_H_m": r_H,
            "equation": (
                "Blazar ergosphere phonon resonance:\n"
                f"  δ_D = {delta_D:.4f} (Doppler factor)\n"
                f"  ω_obs = {omega_obs:.6e} rad/s\n"
                f"  Ω_H = {Omega_H:.6e} rad/s\n"
                f"  Superradiant: {is_superradiant}\n"
                f"  E_ergo = {E_ergo:.6e} J"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  BLAZAR JET POWER Γ-STEPPED CURVES
# ══════════════════════════════════════════════════════════════════════════════

class BlazarJetPowerGammaCurve:
    """
    Compute blazar jet power P_jet(Γ) at stepped linewidths with
    relativistic beaming correction:

        P_jet^obs = P_BZ · (1 + M_jet(Γ)) · δ_D^(p+α) / (1+z)

    where p=2 (continuous jet), α=spectral index.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:          BH mass (kg)
          a_spin:     dimensionless spin
          B:          magnetic field (T)
          Gamma_bulk: bulk Lorentz factor
          theta_obs:  observer angle (rad)
          z:          redshift
          alpha:      spectral index (default 0.7)
          A_jet:      modulation amplitude
        """
        M           = dataset.get('M', 1e9 * M_sun)
        a_spin      = dataset.get('a_spin', 0.95)
        B           = dataset.get('B', 1.0)
        Gamma_bulk  = dataset.get('Gamma_bulk', 10.0)
        theta_obs   = dataset.get('theta_obs', 0.1)
        z           = dataset.get('z', 0.1)
        alpha       = dataset.get('alpha', 0.7)
        A_jet       = dataset.get('A_jet', A_JET_DEFAULT)

        beta = math.sqrt(1.0 - 1.0 / max(Gamma_bulk**2, 1.0 + 1e-50))
        delta_D = 1.0 / (Gamma_bulk * (1.0 - beta * math.cos(theta_obs)))
        p = 2  # continuous jet
        beaming = delta_D**(p + alpha) / (1.0 + z)

        P_BZ = _p_bz(M, a_spin, B)

        results = []
        for Gamma, label in zip(GAMMA_STEPS, GAMMA_LABELS):
            M_j = _m_jet(Gamma, A_jet)
            P_jet_intrinsic = P_BZ * (1.0 + M_j)
            P_jet_obs = P_jet_intrinsic * beaming
            enhancement = P_jet_obs / max(P_BZ, 1e-50)

            results.append({
                "Gamma_label": label,
                "Gamma_rad_s": Gamma,
                "M_jet": M_j,
                "P_BZ_W": P_BZ,
                "P_jet_intrinsic_W": P_jet_intrinsic,
                "P_jet_observed_W": P_jet_obs,
                "enhancement": enhancement,
            })

        return {
            "sweep": results,
            "delta_D": delta_D,
            "beaming_factor": beaming,
            "primary_equations": [
                "P_jet^obs = P_BZ · (1 + M_jet(Γ)) · δ_D^(p+α) / (1+z)",
                "M_jet(Γ) = 1 + A_jet · exp[-(Γ−Γ₀)²/(2σ_Γ²)]",
            ],
            "equation": (
                "Blazar jet power stepped curves:\n"
                + "\n".join(
                    f"  Γ={r['Gamma_label']}: M_jet={r['M_jet']:.4f}, "
                    f"P_jet^obs={r['P_jet_observed_W']:.6e} W"
                    for r in results
                )
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  BLAZAR MULTI-MESSENGER PHONON CORRELATION
# ══════════════════════════════════════════════════════════════════════════════

class BlazarMultiMessengerPhononCorrelation:
    """
    Phonon-enhanced multi-messenger correlation for blazars:

    VHE γ-ray luminosity correlation via phonon-boosted pair cascades:
        L_VHE ∝ P_jet · (1 + Φ · S₂₆) · δ_D^4

    Neutrino luminosity from pγ interactions in phonon-enhanced photon field:
        L_ν ∝ L_VHE · f_pγ · (1 + Φ · S₂₆ · [SSq]/N)

    where f_pγ is the photo-pion production efficiency.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:          BH mass (kg)
          a_spin:     spin
          B:          magnetic field (T)
          Gamma_bulk: bulk Lorentz factor
          theta_obs:  observer angle (rad)
          f_pg:       photo-pion efficiency (default 0.05)
          omega:      phonon frequency (rad/s)
          Gamma:      phonon linewidth (rad/s)
          ssq:        [SSq]
        """
        M           = dataset.get('M', 1e9 * M_sun)
        a_spin      = dataset.get('a_spin', 0.95)
        B           = dataset.get('B', 1.0)
        Gamma_bulk  = dataset.get('Gamma_bulk', 10.0)
        theta_obs   = dataset.get('theta_obs', 0.1)
        f_pg        = dataset.get('f_pg', 0.05)
        omega       = dataset.get('omega', OMEGA_PHONON)
        Gamma       = dataset.get('Gamma', GAMMA_PHONON)
        ssq         = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)

        beta = math.sqrt(1.0 - 1.0 / max(Gamma_bulk**2, 1.0 + 1e-50))
        delta_D = 1.0 / (Gamma_bulk * (1.0 - beta * math.cos(theta_obs)))

        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        P_BZ = _p_bz(M, a_spin, B)
        M_j = _m_jet(Gamma)
        P_jet = P_BZ * (1.0 + M_j)

        # VHE γ-ray luminosity (phonon-enhanced pair cascade)
        L_VHE = P_jet * (1.0 + Phi_norm * S26) * delta_D**4
        L_VHE_no_phonon = P_BZ * 2.0 * delta_D**4

        # Neutrino luminosity (pγ interaction)
        neutrino_enhancement = 1.0 + Phi_norm * S26 * ssq / N_LEVELS
        L_nu = L_VHE * f_pg * neutrino_enhancement
        L_nu_no_phonon = L_VHE_no_phonon * f_pg

        # VHE-neutrino correlation coefficient
        r_corr = L_nu / max(L_VHE, 1e-50)

        return {
            "P_jet_W": P_jet,
            "L_VHE_W": L_VHE,
            "L_VHE_no_phonon_W": L_VHE_no_phonon,
            "L_VHE_enhancement": L_VHE / max(L_VHE_no_phonon, 1e-50),
            "L_nu_W": L_nu,
            "L_nu_no_phonon_W": L_nu_no_phonon,
            "L_nu_enhancement": L_nu / max(L_nu_no_phonon, 1e-50),
            "r_VHE_nu": r_corr,
            "delta_D": delta_D,
            "Phi_norm": Phi_norm,
            "primary_equations": [
                "L_VHE ∝ P_jet · (1 + Φ·S₂₆) · δ_D⁴",
                "L_ν ∝ L_VHE · f_pγ · (1 + Φ·S₂₆·[SSq]/N)",
                "r_corr = L_ν / L_VHE",
            ],
            "equation": (
                "Blazar multi-messenger phonon correlation:\n"
                f"  P_jet = {P_jet:.6e} W\n"
                f"  L_VHE = {L_VHE:.6e} W (enhancement: "
                f"{L_VHE / max(L_VHE_no_phonon, 1e-50):.2f}×)\n"
                f"  L_ν   = {L_nu:.6e} W (f_pγ = {f_pg})\n"
                f"  r_corr = {r_corr:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate blazar jet phonon coupling engine."""
    print("=" * 72)
    print("Blazar Jet Phonon Coupling Engine — Session 212")
    print("=" * 72)

    # §1 Ergosphere resonance (BL Lac type)
    print("\n── §1 Blazar Ergosphere Phonon Resonance (BL Lac) ──")
    ergo = BlazarErgosphereResonance()
    result = ergo.compute({'M': 3e8 * M_sun, 'Gamma_bulk': 15.0, 'theta_obs': 0.05})
    print(f"  δ_D = {result['delta_D']:.4f}")
    print(f"  Superradiant: {result['is_superradiant']}")
    print(f"  E_ergo = {result['E_ergo_J']:.6e} J")

    # §2 Jet power Γ-curve (Mrk 421 type)
    print("\n── §2 Blazar Jet Power Γ-Curves (Mrk 421) ──")
    jet = BlazarJetPowerGammaCurve()
    result = jet.compute({'M': 2e8 * M_sun, 'Gamma_bulk': 20.0, 'B': 0.5, 'z': 0.031})
    for r in result["sweep"]:
        print(f"  Γ={r['Gamma_label']}: P_jet^obs = {r['P_jet_observed_W']:.6e} W ({r['enhancement']:.2f}×)")

    # §3 Multi-messenger correlation
    print("\n── §3 Multi-Messenger Phonon Correlation ──")
    mm = BlazarMultiMessengerPhononCorrelation()
    result = mm.compute({'M': 3e8 * M_sun, 'Gamma_bulk': 15.0})
    print(f"  L_VHE = {result['L_VHE_W']:.6e} W (×{result['L_VHE_enhancement']:.2f})")
    print(f"  L_ν   = {result['L_nu_W']:.6e} W (×{result['L_nu_enhancement']:.2f})")
    print(f"  r_corr = {result['r_VHE_nu']:.6e}")

    print(f"\n{'=' * 72}")
    print("BLAZAR JET PHONON COUPLING — Session 211+212+213")
    print(f"{'=' * 72}")


# ═══════════════════════════════════════════════════════════════════════════════
# §6  CENTAURUS A JET PHONON COUPLING (Session 213)
# ═══════════════════════════════════════════════════════════════════════════════

class CentaurusAJetPhononCoupling:
    """Centaurus A (NGC 5128) phonon-modulated jet power.

    M_BH = 5.5×10⁷ M☉, a = 0.70, B = 3000 T.
    FR I/II transition radio galaxy, nearest VHE gamma-ray source.
    P_BZ ≈ 10⁴³ erg/s → enhancements 2.6/2.1/1.4× at Γ = 0.05/0.10/0.30.
    """

    M_BH_MSUN = 5.5e7
    A_SPIN    = 0.70
    B_T       = 3000
    A_JET     = 0.95

    def compute(self, dataset: dict = None) -> dict:
        import math
        d = dataset or {}
        M = float(d.get("M_Msun", self.M_BH_MSUN)) * 1.989e30
        a = float(d.get("a_spin", self.A_SPIN))
        B = float(d.get("B_T", self.B_T))
        A_jet = float(d.get("A_jet", self.A_JET))

        r_S = 2 * 6.674e-11 * M / (3e8)**2
        r_H = r_S / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
        P_BZ = (B**2 / (8 * math.pi)) * (r_H / 3e8)**2 * a**2 * 3e8

        omega_SCm = 2 * math.pi * 1.25e12
        Gamma_0 = 2 * math.pi * 0.1e12
        sigma_G = 0.08 * 2 * math.pi * 1e12

        gammas = [0.05, 0.10, 0.30]
        curves = []
        for gTHz in gammas:
            Gamma = 2 * math.pi * gTHz * 1e12
            delta = Gamma - Gamma_0
            M_j = 1 + A_jet * math.exp(-delta**2 / (2 * sigma_G**2))
            P_jet = P_BZ * (1 + M_j)
            curves.append({"Gamma_THz": gTHz, "enhancement": P_jet / max(P_BZ, 1e-50)})

        return {
            "system": "Centaurus_A",
            "P_BZ_W": P_BZ,
            "curves": curves,
            "primary_equations": [
                f"P_BZ = {P_BZ:.6e} W ({P_BZ * 1e7:.6e} erg/s)",
            ] + [f"Γ={c['Gamma_THz']}: {c['enhancement']:.1f}×" for c in curves],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §7  TXS 0506+056 JET PHONON COUPLING (Session 213)
# ═══════════════════════════════════════════════════════════════════════════════

class TXS0506JetPhononCoupling:
    """TXS 0506+056 phonon-modulated jet power.

    M_BH = 3×10⁸ M☉, a = 0.85, B = 8000 T.
    BL Lac type blazar, IceCube neutrino association (IceCube-170922A).
    P_BZ ≈ 10⁴⁵ erg/s → enhancements 2.9/2.3/1.6× at Γ = 0.05/0.10/0.30.
    """

    M_BH_MSUN = 3e8
    A_SPIN    = 0.85
    B_T       = 8000
    A_JET     = 1.20

    def compute(self, dataset: dict = None) -> dict:
        import math
        d = dataset or {}
        M = float(d.get("M_Msun", self.M_BH_MSUN)) * 1.989e30
        a = float(d.get("a_spin", self.A_SPIN))
        B = float(d.get("B_T", self.B_T))
        A_jet = float(d.get("A_jet", self.A_JET))

        r_S = 2 * 6.674e-11 * M / (3e8)**2
        r_H = r_S / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
        P_BZ = (B**2 / (8 * math.pi)) * (r_H / 3e8)**2 * a**2 * 3e8

        omega_SCm = 2 * math.pi * 1.25e12
        Gamma_0 = 2 * math.pi * 0.1e12
        sigma_G = 0.08 * 2 * math.pi * 1e12

        gammas = [0.05, 0.10, 0.30]
        curves = []
        for gTHz in gammas:
            Gamma = 2 * math.pi * gTHz * 1e12
            delta = Gamma - Gamma_0
            M_j = 1 + A_jet * math.exp(-delta**2 / (2 * sigma_G**2))
            P_jet = P_BZ * (1 + M_j)
            curves.append({"Gamma_THz": gTHz, "enhancement": P_jet / max(P_BZ, 1e-50)})

        return {
            "system": "TXS_0506+056",
            "P_BZ_W": P_BZ,
            "icecube": "IceCube-170922A",
            "curves": curves,
            "primary_equations": [
                f"P_BZ = {P_BZ:.6e} W ({P_BZ * 1e7:.6e} erg/s)",
                "IceCube neutrino association: IceCube-170922A",
            ] + [f"Γ={c['Gamma_THz']}: {c['enhancement']:.1f}×" for c in curves],
        }


if __name__ == "__main__":
    main()
