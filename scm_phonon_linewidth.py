#!/usr/bin/env python3
"""
scm_phonon_linewidth.py — SCm Phonon Linewidth Exploration Engine

Session 212 | Daniel Murphy
PURPOSE: Standalone engine for systematic linewidth Γ exploration of the
         1.25 THz SCm phonon resonance, deriving Γ→E_net, Γ→F_neutron,
         buoyancy reversal thresholds, and regime classification.

         Key physics:
           - LinewidthEnetEvolution: E_net(Γ) sweep at 0.05/0.1/0.3 THz
           - LinewidthNeutronDrop: F_neutron(ω,Γ) with Kozima coupling
           - LinewidthBuoyancyReversal: sign-flip threshold vs Γ
           - LinewidthRegimeClassifier: narrow/optimal/broad classification

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
# §0  LINEWIDTH-SPECIFIC CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PHI_0_DEFAULT = 1e20           # phonons/m²/s (THz phonon fluence)
OMEGA_PHONON  = OMEGA_SCM      # 2π × 1.25e12 rad/s  (resonance center)
GAMMA_PHONON  = GAMMA_DEFAULT  # 2π × 0.1e12 rad/s   (resonance linewidth)

# Canonical stepped linewidths (THz → rad/s)
GAMMA_NARROW   = 2.0 * PI * 0.05e12   # 0.05 THz (narrow)
GAMMA_OPTIMAL  = 2.0 * PI * 0.10e12   # 0.10 THz (optimal)
GAMMA_BROAD    = 2.0 * PI * 0.30e12   # 0.30 THz (broad)
GAMMA_STEPS    = [GAMMA_NARROW, GAMMA_OPTIMAL, GAMMA_BROAD]
GAMMA_LABELS   = ["0.05 THz (narrow)", "0.10 THz (optimal)", "0.30 THz (broad)"]

# Kozima neutron-drop parameters
N_NEUTRON_DEFAULT = 1e18       # neutron density (m⁻³)
SIGMA_N_0         = 1e-4       # bare neutron cross-section (m²)


# ══════════════════════════════════════════════════════════════════════════════
# §1  LINEWIDTH → E_net EVOLUTION
# ══════════════════════════════════════════════════════════════════════════════

class LinewidthEnetEvolution:
    """
    Sweep linewidth Γ through 0.05 / 0.10 / 0.30 THz and compute
    E_net(Γ) = ρ_SCm · V · (2F_{U,Bi}/F_U − 1) · Φ(ω,Γ) · S₂₆ for each.

    The phonon fluence Φ(ω,Γ) = Φ₀ · exp[−(ω−ω_SCm)²/(2Γ²)] · S₂₆
    narrows with decreasing Γ, concentrating energy near resonance.

    Returns E_net at each stepped Γ with quality factor Q = ω_SCm / (2Γ).
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          V:        region volume (m³, default 1e48)
          F_U_Bi:   buoyancy fraction (default 0.6)
          F_U:      total unified force (default 1.0)
          t:        time (s, default 0)
          omega:    probe frequency (rad/s, default ω_SCm)
          ssq:      [SSq] (default 0.57)
        """
        V       = dataset.get('V', 1e48)
        F_U_Bi  = dataset.get('F_U_Bi', 0.6)
        F_U     = dataset.get('F_U', 1.0)
        t       = dataset.get('t', 0.0)
        omega   = dataset.get('omega', OMEGA_PHONON)
        ssq     = dataset.get('ssq', SSQ)

        net_factor = 2.0 * F_U_Bi / max(F_U, 1e-50) - 1.0
        S26 = S26_accelerated(ssq)
        rho_SCm_t = RHO_VAC_SCM * S26 * math.exp(KAPPA * t + ssq * t / 26.0)

        results = []
        for Gamma, label in zip(GAMMA_STEPS, GAMMA_LABELS):
            delta_omega = omega - OMEGA_PHONON
            gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
            Phi = PHI_0_DEFAULT * gaussian * S26
            E_net = rho_SCm_t * V * net_factor * Phi * S26
            Q_factor = OMEGA_PHONON / (2.0 * Gamma)
            tau_res = 1.0 / (2.0 * PI * Gamma / (2.0 * PI))  # 1/(Γ in Hz)
            results.append({
                "Gamma_label": label,
                "Gamma_rad_s": Gamma,
                "Phi": Phi,
                "E_net_J": E_net,
                "Q_factor": Q_factor,
                "tau_res_s": tau_res,
            })

        return {
            "sweep": results,
            "rho_SCm_t": rho_SCm_t,
            "S26": S26,
            "net_factor": net_factor,
            "primary_equations": [
                "E_net(Γ) = ρ_SCm(t) · V · (2F_{U,Bi}/F_U − 1) · Φ(ω,Γ) · S₂₆",
                "Φ(ω,Γ) = Φ₀ · exp[−(ω−ω_SCm)²/(2Γ²)] · S₂₆",
                "Q = ω_SCm / (2Γ)",
            ],
            "equation": (
                "Linewidth → E_net stepped sweep:\n"
                + "\n".join(
                    f"  Γ={r['Gamma_label']}: E_net={r['E_net_J']:.6e} J, Q={r['Q_factor']:.1f}"
                    for r in results
                )
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  LINEWIDTH → KOZIMA NEUTRON-DROP FORCE
# ══════════════════════════════════════════════════════════════════════════════

class LinewidthNeutronDrop:
    """
    Compute the Kozima neutron-drop force F_neutron(ω,Γ) with explicit
    linewidth dependency:

        σ_n(ω,Γ) = σ₀ · exp[−(ω−ω_SCm)²/(2Γ²)] · (1 + [SSq]·n/N)
        F_neutron(Γ) = N_n · σ_n(ω,Γ) · Φ(ω,Γ) · E_net(t,Γ)

    At narrow Γ (0.05 THz): σ_n peaks sharply, F_neutron maximized
    At broad Γ (0.30 THz):  σ_n broadened, F_neutron reduced
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:    probe frequency (rad/s, default ω_SCm)
          N_n:      neutron density (m⁻³, default 1e18)
          sigma_0:  bare cross-section (m², default 1e-4)
          n_level:  current lattice level (default 13)
          ssq:      [SSq] (default 0.57)
          E_net_bare: pre-computed E_net (J, default 1e20)
        """
        omega      = dataset.get('omega', OMEGA_PHONON)
        N_n        = dataset.get('N_n', N_NEUTRON_DEFAULT)
        sigma_0    = dataset.get('sigma_0', SIGMA_N_0)
        n_level    = dataset.get('n_level', 13)
        ssq        = dataset.get('ssq', SSQ)
        E_net_bare = dataset.get('E_net_bare', 1e20)

        S26 = S26_accelerated(ssq)
        results = []

        for Gamma, label in zip(GAMMA_STEPS, GAMMA_LABELS):
            delta_omega = omega - OMEGA_PHONON
            exponent = -(delta_omega**2) / (2.0 * Gamma**2)
            gaussian = math.exp(min(exponent, 0.0))
            vds_factor = 1.0 + ssq * n_level / N_LEVELS

            sigma_n = sigma_0 * gaussian * vds_factor
            Phi = PHI_0_DEFAULT * gaussian * S26
            F_neutron = N_n * sigma_n * Phi * E_net_bare
            F_ratio = F_neutron / max(abs(E_net_bare), 1e-50)

            results.append({
                "Gamma_label": label,
                "Gamma_rad_s": Gamma,
                "sigma_n": sigma_n,
                "Phi": Phi,
                "F_neutron_N": F_neutron,
                "F_ratio": F_ratio,
            })

        return {
            "sweep": results,
            "N_n": N_n,
            "sigma_0": sigma_0,
            "primary_equations": [
                "σ_n(ω,Γ) = σ₀ · exp[−(ω−ω_SCm)²/(2Γ²)] · (1 + [SSq]·n/N)",
                "F_neutron(Γ) = N_n · σ_n(ω,Γ) · Φ(ω,Γ) · E_net",
            ],
            "equation": (
                "Kozima neutron-drop vs. linewidth:\n"
                + "\n".join(
                    f"  Γ={r['Gamma_label']}: F_neutron={r['F_neutron_N']:.6e} N, "
                    f"σ_n={r['sigma_n']:.6e} m²"
                    for r in results
                )
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  LINEWIDTH → BUOYANCY REVERSAL THRESHOLD
# ══════════════════════════════════════════════════════════════════════════════

class LinewidthBuoyancyReversal:
    """
    Find the buoyancy sign-flip threshold as function of linewidth Γ.

    The net energy E_net changes sign when:
        2 F_{U,Bi}/F_U − 1 = F_neutron / (ρ_SCm · V · Φ · S₂₆)

    At each Γ, compute the critical F_{U,Bi}/F_U ratio where E_net = 0
    and classify whether the system is expanding or eroding.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          V:          region volume (m³, default 1e48)
          F_U:        total unified force (default 1.0)
          omega:      probe frequency (rad/s, default ω_SCm)
          ssq:        [SSq] (default 0.57)
          N_n:        neutron density (m⁻³, default 1e18)
          sigma_0:    bare cross-section (m², default 1e-4)
          n_level:    lattice level (default 13)
        """
        V         = dataset.get('V', 1e48)
        F_U       = dataset.get('F_U', 1.0)
        omega     = dataset.get('omega', OMEGA_PHONON)
        ssq       = dataset.get('ssq', SSQ)
        N_n       = dataset.get('N_n', N_NEUTRON_DEFAULT)
        sigma_0   = dataset.get('sigma_0', SIGMA_N_0)
        n_level   = dataset.get('n_level', 13)

        S26 = S26_accelerated(ssq)
        rho_SCm = RHO_VAC_SCM * S26

        results = []
        for Gamma, label in zip(GAMMA_STEPS, GAMMA_LABELS):
            delta_omega = omega - OMEGA_PHONON
            gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
            Phi = PHI_0_DEFAULT * gaussian * S26

            vds_factor = 1.0 + ssq * n_level / N_LEVELS
            sigma_n = sigma_0 * gaussian * vds_factor

            # Critical ratio where E_net = 0:
            # net_factor = F_neutron / (rho · V · Phi · S26)
            # 2 F_UBi/F_U − 1 = (N_n · sigma_n · Phi) / (rho · V · S26)
            denominator = rho_SCm * V * S26
            if denominator > 0:
                critical_net_factor = N_n * sigma_n * Phi / denominator
                critical_ratio = (critical_net_factor + 1.0) / 2.0
            else:
                critical_ratio = 0.5

            is_physical = 0.0 < critical_ratio < 1.0

            results.append({
                "Gamma_label": label,
                "Gamma_rad_s": Gamma,
                "critical_FUBi_FU": critical_ratio,
                "is_physical": is_physical,
                "classification": (
                    "expansion" if critical_ratio > 0.5 else "erosion"
                ),
            })

        return {
            "sweep": results,
            "rho_SCm": rho_SCm,
            "primary_equations": [
                "Sign flip at: 2F_{U,Bi}/F_U − 1 = F_neutron/(ρ·V·Φ·S₂₆)",
                "Critical ratio: F_{U,Bi}/F_U = (critical_net_factor + 1)/2",
            ],
            "equation": (
                "Buoyancy reversal threshold vs. linewidth:\n"
                + "\n".join(
                    f"  Γ={r['Gamma_label']}: critical F_UBi/F_U = {r['critical_FUBi_FU']:.6f}"
                    for r in results
                )
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  LINEWIDTH REGIME CLASSIFIER
# ══════════════════════════════════════════════════════════════════════════════

class LinewidthRegimeClassifier:
    """
    Classify the phonon linewidth regime:
      - Narrow (Γ < 0.07 THz): Q > 8.9, sharp resonance, max F_neutron
      - Optimal (0.07 ≤ Γ ≤ 0.15 THz): Q ~ 4–9, balanced coupling
      - Broad (Γ > 0.15 THz): Q < 4.2, delocalized phonon field

    Outputs regime label, quality factor, and coupling efficiency
    η = F_neutron(Γ) / F_neutron(Γ_optimal).
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          Gamma:     linewidth (rad/s, default Γ_SCm)
          omega:     probe frequency (rad/s, default ω_SCm)
          ssq:       [SSq] (default 0.57)
        """
        Gamma_in = dataset.get('Gamma', GAMMA_PHONON)
        omega    = dataset.get('omega', OMEGA_PHONON)
        ssq      = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)
        Gamma_THz = Gamma_in / (2.0 * PI * 1e12)
        Q = OMEGA_PHONON / (2.0 * Gamma_in)

        # Classify regime
        if Gamma_THz < 0.07:
            regime = "narrow"
            description = "Sharp resonance, maximum neutron-drop force, high Q"
        elif Gamma_THz <= 0.15:
            regime = "optimal"
            description = "Balanced phonon-gravity coupling, moderate Q"
        else:
            regime = "broad"
            description = "Delocalized phonon field, reduced specificity, low Q"

        # Coupling efficiency relative to optimal Γ
        delta_opt = omega - OMEGA_PHONON
        gaussian_opt = math.exp(-delta_opt**2 / (2.0 * GAMMA_OPTIMAL**2))
        delta_in = omega - OMEGA_PHONON
        gaussian_in = math.exp(-delta_in**2 / (2.0 * Gamma_in**2))
        eta = gaussian_in / max(gaussian_opt, 1e-50)

        # Damping timescale
        tau_damping = 1.0 / (Gamma_in / (2.0 * PI)) if Gamma_in > 0 else float('inf')

        return {
            "Gamma_rad_s": Gamma_in,
            "Gamma_THz": Gamma_THz,
            "Q_factor": Q,
            "regime": regime,
            "description": description,
            "eta_coupling": eta,
            "tau_damping_s": tau_damping,
            "primary_equations": [
                "Q = ω_SCm / (2Γ)",
                "η = Φ(ω,Γ) / Φ(ω,Γ_optimal)",
                "τ_damping = 2π / Γ",
            ],
            "equation": (
                f"Linewidth regime classification:\n"
                f"  Γ = {Gamma_THz:.3f} THz → regime = {regime}\n"
                f"  Q = {Q:.2f}, η = {eta:.6f}\n"
                f"  τ_damping = {tau_damping:.6e} s\n"
                f"  {description}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate SCm phonon linewidth exploration engine."""
    print("=" * 72)
    print("SCm Phonon Linewidth Exploration Engine — Session 212")
    print("=" * 72)

    # §1 E_net evolution
    print("\n── §1 Linewidth → E_net Evolution ──")
    enet = LinewidthEnetEvolution()
    result = enet.compute({})
    for r in result["sweep"]:
        print(f"  Γ={r['Gamma_label']}: E_net={r['E_net_J']:.6e} J, Q={r['Q_factor']:.1f}")

    # §2 Neutron-drop force
    print("\n── §2 Linewidth → Kozima Neutron-Drop Force ──")
    neutron = LinewidthNeutronDrop()
    result = neutron.compute({})
    for r in result["sweep"]:
        print(f"  Γ={r['Gamma_label']}: F_neutron={r['F_neutron_N']:.6e} N")

    # §3 Buoyancy reversal
    print("\n── §3 Linewidth → Buoyancy Reversal ──")
    reversal = LinewidthBuoyancyReversal()
    result = reversal.compute({})
    for r in result["sweep"]:
        print(f"  Γ={r['Gamma_label']}: critical F_UBi/F_U = {r['critical_FUBi_FU']:.6f}")

    # §4 Regime classifier
    print("\n── §4 Linewidth Regime Classifier ──")
    classifier = LinewidthRegimeClassifier()
    for Gamma_test in GAMMA_STEPS:
        result = classifier.compute({"Gamma": Gamma_test})
        print(f"  Γ={result['Gamma_THz']:.3f} THz → {result['regime']} (Q={result['Q_factor']:.1f})")

    print(f"\n{'=' * 72}")
    print("SCm PHONON LINEWIDTH EXPLORATION COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
