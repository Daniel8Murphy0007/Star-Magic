#!/usr/bin/env python3
"""
quasar_jet_phonon.py — Quasar Jet Phonon Modulation Engine

Session 211 | Daniel Murphy
PURPOSE: Standalone engine for phonon-modulated jet power in quasar/AGN jets.

         Key physics:
           - M_jet(Γ) = 1 + A_jet · exp[-(Γ−Γ₀)²/(2σ_Γ²)]
             Jet modulation factor as function of phonon linewidth Γ
           - P_jet = P_BZ · (1 + M_jet(Γ))
             Total jet power including phonon enhancement
           - P_BZ = (B²/8π) · (r_H/c)² · a²·c      (Blandford-Znajek)
           - WSTP export: Wolfram Language expressions for M_jet[Γ_]/P_jet[Γ_]

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
# §0  JET-PHONON CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PHI_0_DEFAULT  = 1e20
OMEGA_PHONON   = OMEGA_SCM
GAMMA_PHONON   = GAMMA_DEFAULT
A_JET_DEFAULT  = 1.5         # jet modulation amplitude
SIGMA_GAMMA_DEFAULT = 0.08   # σ_Γ in THz for M_jet(Γ)


# ══════════════════════════════════════════════════════════════════════════════
# §1  JET MODULATION FACTOR M_jet(Γ)
# ══════════════════════════════════════════════════════════════════════════════

class JetModulationFactor:
    """
    Compute the jet modulation factor M_jet(Γ) as function of phonon linewidth:

        M_jet(Γ) = 1 + A_jet · exp[-(Γ − Γ₀)²/(2σ_Γ²)]

    Parameters:
      A_jet:   modulation amplitude (dimensionless, default 1.5)
      Γ₀:      reference linewidth (default Γ_SCm = 2π × 0.1 THz)
      σ_Γ:     Gaussian width of modulation (default 0.08 THz)

    At Γ = Γ₀: M_jet = 1 + A_jet (maximum modulation ~2.5×)
    Far from Γ₀: M_jet → 1 (no modulation)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          Gamma:         phonon linewidth (rad/s, default Γ_SCm)
          A_jet:         amplitude (default 1.5)
          sigma_Gamma_THz: width in THz (default 0.08)
          Gamma_0:       reference linewidth (rad/s, default Γ_SCm)
        """
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        A_jet   = dataset.get('A_jet', A_JET_DEFAULT)
        sigma_G = dataset.get('sigma_Gamma_THz', SIGMA_GAMMA_DEFAULT) * 2.0 * PI * 1e12
        Gamma_0 = dataset.get('Gamma_0', GAMMA_PHONON)

        delta_G = Gamma - Gamma_0
        gauss = math.exp(-delta_G**2 / (2.0 * sigma_G**2))
        M_jet = 1.0 + A_jet * gauss

        return {
            "M_jet": M_jet,
            "A_jet": A_jet,
            "gauss_factor": gauss,
            "Gamma_rad_s": Gamma,
            "Gamma_THz": Gamma / (2.0 * PI * 1e12),
            "Gamma_0_THz": Gamma_0 / (2.0 * PI * 1e12),
            "sigma_Gamma_THz": sigma_G / (2.0 * PI * 1e12),
            "equation": (
                "M_jet(Γ) = 1 + A_jet · exp[-(Γ−Γ₀)²/(2σ_Γ²)]\n"
                f"         = 1 + {A_jet:.2f} × {gauss:.6f}\n"
                f"         = {M_jet:.6f}"
            ),
        }

    def sweep(self, dataset: dict) -> Dict[str, Any]:
        """Sweep Γ from Gamma_min to Gamma_max and compute M_jet at each point."""
        Gamma_min = dataset.get('Gamma_min_THz', 0.02) * 2.0 * PI * 1e12
        Gamma_max = dataset.get('Gamma_max_THz', 0.30) * 2.0 * PI * 1e12
        n_steps   = dataset.get('n_steps', 50)

        results = []
        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            Gamma_i = Gamma_min + frac * (Gamma_max - Gamma_min)

            sub = dict(dataset)
            sub['Gamma'] = Gamma_i
            result = self.compute(sub)

            results.append({
                "Gamma_THz": result["Gamma_THz"],
                "M_jet": result["M_jet"],
            })

        peak = max(results, key=lambda x: x["M_jet"])
        return {
            "sweep": results,
            "peak_M_jet": peak["M_jet"],
            "peak_Gamma_THz": peak["Gamma_THz"],
            "n_steps": n_steps,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  BLANDFORD-ZNAJEK JET POWER
# ══════════════════════════════════════════════════════════════════════════════

class BlandfordZnajekJetPower:
    """
    Blandford-Znajek jet power with phonon modulation:

        P_BZ = (B² / 8π) · (r_H / c)² · a² · c
        P_jet = P_BZ · (1 + M_jet(Γ))

    where B is the magnetic field at the horizon, r_H is the horizon radius,
    and a is the dimensionless spin parameter.

    The phonon enhancement factor (1 + M_jet) gives total modulation
    ranging from ~2× (far off-resonance) to ~3.5× (at Γ = Γ₀).
    """

    def __init__(self):
        self._mjet = JetModulationFactor()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:         BH mass (kg)
          a_spin:    dimensionless spin (0-1, default 0.9)
          B_horizon: magnetic field at horizon (T, default 1e4)
          Gamma:     phonon linewidth (rad/s)
          A_jet, sigma_Gamma_THz: jet modulation parameters
        """
        M       = dataset.get('M', 1e9 * M_sun)
        a_spin  = dataset.get('a_spin', 0.9)
        B       = dataset.get('B_horizon', 1e4)

        # Horizon radius
        r_S = 2.0 * G * M / c**2
        r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - a_spin**2, 0.0)))

        # BZ jet power
        P_BZ = (B**2 / (8.0 * PI)) * (r_H / c)**2 * a_spin**2 * c

        # Jet modulation
        mjet_result = self._mjet.compute(dataset)
        M_jet = mjet_result["M_jet"]

        # Total phonon-enhanced jet power
        P_jet = P_BZ * (1.0 + M_jet)

        # Convert to erg/s for astrophysical convention
        P_BZ_erg = P_BZ * 1e7
        P_jet_erg = P_jet * 1e7

        # Eddington luminosity for comparison
        sigma_T = 6.652e-29  # Thomson cross-section
        m_p = 1.673e-27
        L_Edd = 4.0 * PI * G * M * m_p * c / sigma_T
        L_Edd_erg = L_Edd * 1e7

        return {
            "P_BZ_W": P_BZ,
            "P_BZ_erg_s": P_BZ_erg,
            "P_jet_W": P_jet,
            "P_jet_erg_s": P_jet_erg,
            "M_jet": M_jet,
            "modulation_factor": 1.0 + M_jet,
            "L_Edd_erg_s": L_Edd_erg,
            "P_jet_over_L_Edd": P_jet / max(L_Edd, 1e-50),
            "r_H_m": r_H,
            "B_horizon_T": B,
            "a_spin": a_spin,
            "equation": (
                "P_BZ = (B²/8π) · (r_H/c)² · a² · c\n"
                f"     = {P_BZ_erg:.4e} erg/s\n"
                f"P_jet = P_BZ · (1 + M_jet) = {P_jet_erg:.4e} erg/s\n"
                f"Modulation: ×{1.0 + M_jet:.2f}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  PHONON JET POWER Γ-SWEEP
# ══════════════════════════════════════════════════════════════════════════════

class JetPowerGammaSweep:
    """
    Sweep phonon linewidth Γ and compute jet power P_jet(Γ) to produce
    jet power curves for comparison with observations.
    """

    def __init__(self):
        self._bz = BlandfordZnajekJetPower()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          Gamma_min_THz, Gamma_max_THz: sweep range
          n_steps: number of Γ values
          M, a_spin, B_horizon: BH/jet parameters
          A_jet, sigma_Gamma_THz: modulation parameters
        """
        Gamma_min = dataset.get('Gamma_min_THz', 0.02) * 2.0 * PI * 1e12
        Gamma_max = dataset.get('Gamma_max_THz', 0.30) * 2.0 * PI * 1e12
        n_steps   = dataset.get('n_steps', 50)

        sweep = []
        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            Gamma_i = Gamma_min + frac * (Gamma_max - Gamma_min)

            sub = dict(dataset)
            sub['Gamma'] = Gamma_i

            result = self._bz.compute(sub)
            sweep.append({
                "Gamma_THz": Gamma_i / (2.0 * PI * 1e12),
                "P_jet_erg_s": result["P_jet_erg_s"],
                "P_BZ_erg_s": result["P_BZ_erg_s"],
                "M_jet": result["M_jet"],
                "modulation": result["modulation_factor"],
            })

        peak = max(sweep, key=lambda x: x["P_jet_erg_s"])
        base_P = sweep[-1]["P_BZ_erg_s"]  # off-resonance baseline

        mod_range_min = min(s["modulation"] for s in sweep)
        mod_range_max = max(s["modulation"] for s in sweep)

        return {
            "sweep": sweep,
            "peak_P_jet_erg_s": peak["P_jet_erg_s"],
            "peak_Gamma_THz": peak["Gamma_THz"],
            "P_BZ_baseline_erg_s": base_P,
            "modulation_range": [round(mod_range_min, 2), round(mod_range_max, 2)],
            "n_steps": n_steps,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  WSTP EXPRESSION EXPORT
# ══════════════════════════════════════════════════════════════════════════════

class JetPhononWSTPExporter:
    """
    Generate Wolfram Language expressions for M_jet[Γ_] and P_jet[Γ_]
    suitable for insertion into wstp_kernel_demo_runner.py.
    """

    def export_expressions(self, dataset: dict = None) -> List[Dict[str, str]]:
        """
        Build WL expressions for jet-phonon modulation.

        Returns list of {"label": ..., "code": ...} entries.
        """
        dataset = dataset or {}
        A_jet = dataset.get('A_jet', A_JET_DEFAULT)
        sigma_G_THz = dataset.get('sigma_Gamma_THz', SIGMA_GAMMA_DEFAULT)
        Gamma_0_THz = GAMMA_PHONON / (2.0 * PI * 1e12)

        exprs = []

        # M_jet[Γ] definition
        exprs.append({
            "label": "Jet modulation factor M_jet[Γ]",
            "code": (
                f'Ajet = {A_jet}; Gamma0THz = {Gamma_0_THz:.4f}; '
                f'sigmaGTHz = {sigma_G_THz}; '
                'MjetGamma[GammaTHz_] := 1 + Ajet * '
                'Exp[-(GammaTHz - Gamma0THz)^2 / (2 sigmaGTHz^2)]; '
                'MjetGamma[Gamma0THz] // N'
            ),
        })

        # P_jet[Γ] definition
        exprs.append({
            "label": "Phonon-enhanced jet power P_jet[Γ] (generic, erg/s)",
            "code": (
                'PBZ = 10^44; '  # placeholder BZ power
                'PjetGamma[GammaTHz_] := PBZ * (1 + MjetGamma[GammaTHz]); '
                f'PjetGamma[{Gamma_0_THz:.4f}] // N'
            ),
        })

        # Γ-sweep (table)
        exprs.append({
            "label": "P_jet Γ-sweep table (0.02–0.30 THz)",
            "code": (
                'Table[{g, PjetGamma[g]}, {g, 0.02, 0.30, 0.02}]'
            ),
        })

        return exprs


# ══════════════════════════════════════════════════════════════════════════════
# §5  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate quasar jet phonon modulation."""
    print("=" * 72)
    print("Quasar Jet Phonon Modulation Engine — Session 211")
    print("=" * 72)

    # §1 Jet modulation factor
    print("\n── §1 Jet Modulation Factor M_jet(Γ) ──")
    mjet = JetModulationFactor()
    result = mjet.compute({})
    print(f"  M_jet(Γ₀) = {result['M_jet']:.4f}")
    print(f"  At Γ₀ = {result['Gamma_0_THz']:.4f} THz")

    # §2 BZ jet power
    print("\n── §2 Blandford-Znajek Jet Power with Phonon ──")
    bz = BlandfordZnajekJetPower()
    result = bz.compute({'M': 6.5e9 * M_sun, 'B_horizon': 100, 'a_spin': 0.9})
    print(f"  P_BZ = {result['P_BZ_erg_s']:.4e} erg/s")
    print(f"  P_jet = {result['P_jet_erg_s']:.4e} erg/s")
    print(f"  Modulation: ×{result['modulation_factor']:.2f}")

    # §3 Γ-sweep
    print("\n── §3 Jet Power Γ-Sweep ──")
    sweep = JetPowerGammaSweep()
    result = sweep.compute({'M': 6.5e9 * M_sun, 'B_horizon': 100, 'a_spin': 0.9})
    print(f"  Peak P_jet = {result['peak_P_jet_erg_s']:.4e} erg/s at Γ = {result['peak_Gamma_THz']:.3f} THz")
    print(f"  Modulation range: {result['modulation_range']}")

    # §4 WSTP expressions
    print("\n── §4 WSTP Expression Export ──")
    exporter = JetPhononWSTPExporter()
    exprs = exporter.export_expressions()
    for e in exprs:
        print(f"  [{e['label']}]")
        print(f"    {e['code'][:80]}...")

    print(f"\n{'=' * 72}")
    print("QUASAR JET PHONON MODULATION COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
