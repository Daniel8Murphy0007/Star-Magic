#!/usr/bin/env python3
"""
scm_phonon_resonance.py — SCm Phonon Resonance Exploration Engine

Session 211 | Daniel Murphy
PURPOSE: Standalone exploration engine for the 1.25 THz SCm phonon resonance
         in the superconductive vacuum manifold.

         Derives the resonance acceleration term a_res and explores linewidth
         Γ effects, vacuum density coupling strength, and resonance Q-factor.

         Key physics:
           - a_res = (F_{U,Bi}/F_U) · Φ_{1.25 THz}(ω) · S₂₆([SSq])
           - Linewidth sweep: Γ ∈ [0.05, 0.3] THz → quality factor Q(Γ)
           - Vacuum coupling: C_vac = ρ_SCm / ρ_UA × Φ × S₂₆
           - Phonon pressure: P_phonon = n_phonon · ℏω · Φ / V
           - Resonance damping: τ_res = 1 / (2πΓ)

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

PHI_0_DEFAULT = 1e20          # phonons/m²/s (THz phonon fluence)
OMEGA_PHONON  = OMEGA_SCM     # 2π × 1.25e12 rad/s  (resonance center)
GAMMA_PHONON  = GAMMA_DEFAULT # 2π × 0.1e12 rad/s   (resonance linewidth)

# Derived constants
Q_RESONANCE = OMEGA_PHONON / (2.0 * GAMMA_PHONON)  # Quality factor ~6.25
TAU_RES     = 1.0 / (2.0 * PI * GAMMA_DEFAULT / (2.0 * PI))  # Damping time ~1e-11 s


# ══════════════════════════════════════════════════════════════════════════════
# §1  RESONANCE ACCELERATION TERM  a_res
# ══════════════════════════════════════════════════════════════════════════════

class ResonanceAcceleration:
    """
    Compute the phonon resonance acceleration term:

        a_res = (F_{U,Bi} / F_U) · Φ_{1.25 THz}(ω) · S₂₆([SSq])

    This is the missing link between the phonon modulation factor Φ and
    the gravitational field equations. When F_{U,Bi}/F_U > 0.5, a_res
    drives expansion; when < 0.5, a_res drives contraction.

    The acceleration has units of m/s² when Φ is dimensionless (normalized
    to Φ₀) and F_{U,Bi}/F_U is the buoyancy-to-gravity ratio.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:     angular frequency (rad/s, default ω_SCm)
          Gamma:     linewidth (rad/s, default Γ_SCm)
          Phi_0:     fluence normalization (default 1e20)
          F_U_Bi:    buoyancy force (N, default 1.0)
          F_U:       gravity force (N, default 1.0)
          ssq:       [SSq] parameter (default 0.57)
          normalize_phi: whether to normalize Φ to [0,1] (default True)
        """
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0   = dataset.get('Phi_0', PHI_0_DEFAULT)
        F_U_Bi  = dataset.get('F_U_Bi', 1.0)
        F_U     = dataset.get('F_U', 1.0)
        ssq     = dataset.get('ssq', SSQ)
        norm    = dataset.get('normalize_phi', True)

        # S₂₆([SSq])
        S26 = S26_accelerated(ssq)

        # Φ_{1.25 THz}(ω) = Φ₀ · exp[-(ω−ω_SCm)²/(2Γ²)] · S₂₆
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_raw = Phi_0 * gaussian * S26

        # Optionally normalize
        Phi_peak = Phi_0 * S26  # peak value at ω = ω_SCm
        Phi_norm = gaussian * S26 if norm else Phi_raw

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

        # Quality factor
        Q = omega / (2.0 * Gamma) if Gamma > 0 else float('inf')

        return {
            "a_res": a_res,
            "F_UBi_over_FU": ratio,
            "Phi_125_THz": Phi_norm,
            "Phi_raw": Phi_raw,
            "Phi_peak": Phi_peak,
            "S26": S26,
            "gaussian": gaussian,
            "Q_factor": Q,
            "regime": regime,
            "equation": (
                "a_res = (F_{U,Bi}/F_U) · Φ_{1.25 THz}(ω) · S₂₆([SSq])\n"
                f"      = {ratio:.6f} × {Phi_norm:.6e} × {S26:.6e}\n"
                f"      = {a_res:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  LINEWIDTH Γ SWEEP
# ══════════════════════════════════════════════════════════════════════════════

class LinewidthGammaSweep:
    """
    Sweep the phonon resonance linewidth Γ from 0.05 THz to 0.3 THz
    and compute the quality factor Q(Γ), damping time τ(Γ), and
    resonance acceleration a_res(Γ).

    Narrow linewidths (high Q) → sharp resonance, large a_res at center
    Broad linewidths (low Q) → diffuse resonance, smaller a_res but wider band
    """

    def __init__(self):
        self._accel = ResonanceAcceleration()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          Gamma_min_THz: minimum Γ in THz (default 0.05)
          Gamma_max_THz: maximum Γ in THz (default 0.30)
          n_steps:       number of Γ values (default 26)
          omega:         probe frequency (default ω_SCm)
          F_U_Bi, F_U:  buoyancy/gravity forces
          Phi_0:         fluence normalization
          ssq:           [SSq] parameter
        """
        Gamma_min = dataset.get('Gamma_min_THz', 0.05) * 2.0 * PI * 1e12
        Gamma_max = dataset.get('Gamma_max_THz', 0.30) * 2.0 * PI * 1e12
        n_steps   = dataset.get('n_steps', 26)
        omega     = dataset.get('omega', OMEGA_PHONON)

        sweep = []
        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            Gamma_i = Gamma_min + frac * (Gamma_max - Gamma_min)
            Gamma_THz = Gamma_i / (2.0 * PI * 1e12)

            sub = dict(dataset)
            sub['Gamma'] = Gamma_i
            sub['omega'] = omega

            result = self._accel.compute(sub)

            Q_i = omega / (2.0 * Gamma_i) if Gamma_i > 0 else float('inf')
            tau_i = 1.0 / Gamma_i if Gamma_i > 0 else float('inf')

            sweep.append({
                "Gamma_THz": round(Gamma_THz, 4),
                "Gamma_rad_s": Gamma_i,
                "Q_factor": round(Q_i, 4),
                "tau_damping_s": tau_i,
                "a_res": result["a_res"],
                "Phi_125": result["Phi_125_THz"],
                "regime": result["regime"],
            })

        # Find peak a_res
        peak = max(sweep, key=lambda x: abs(x["a_res"]))

        return {
            "sweep": sweep,
            "n_steps": n_steps,
            "Gamma_range_THz": [
                Gamma_min / (2.0 * PI * 1e12),
                Gamma_max / (2.0 * PI * 1e12),
            ],
            "peak_a_res": peak["a_res"],
            "peak_Gamma_THz": peak["Gamma_THz"],
            "peak_Q": peak["Q_factor"],
            "description": (
                "Linewidth Γ sweep for SCm 1.25 THz phonon resonance.\n"
                f"Γ range: {sweep[0]['Gamma_THz']:.3f}–{sweep[-1]['Gamma_THz']:.3f} THz\n"
                f"Q range: {sweep[-1]['Q_factor']:.2f}–{sweep[0]['Q_factor']:.2f}\n"
                f"Peak |a_res|: {abs(peak['a_res']):.6e} at Γ = {peak['Gamma_THz']:.3f} THz"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  VACUUM DENSITY COUPLING STRENGTH
# ══════════════════════════════════════════════════════════════════════════════

class VacuumDensityCoupling:
    """
    Compute the vacuum density coupling coefficient:

        C_vac = (ρ_SCm / ρ_UA) · Φ_{1.25 THz}(ω) · S₂₆([SSq])

    This coefficient quantifies how strongly the phonon resonance couples
    to the vacuum energy density ratio. At resonance (ω = ω_SCm), C_vac
    is maximized. Off-resonance, C_vac decays as a Gaussian.

    Also computes the phonon pressure:
        P_phonon = n_phonon · ℏω_SCm · Φ / V
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega:       angular frequency (rad/s, default ω_SCm)
          Gamma:       linewidth (rad/s, default Γ_SCm)
          Phi_0:       fluence normalization
          rho_scm:     SCm vacuum density (kg/m³, default 7.09e-37)
          rho_ua:      UA vacuum density (kg/m³, default 7.09e-36)
          V_region:    spatial volume (m³, default 1e48)
          n_phonon:    phonon occupation number (default 1e26)
          ssq:         [SSq] parameter
        """
        omega    = dataset.get('omega', OMEGA_PHONON)
        Gamma    = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0    = dataset.get('Phi_0', PHI_0_DEFAULT)
        rho_scm  = dataset.get('rho_scm', RHO_SCM)
        rho_ua   = dataset.get('rho_ua', RHO_UA)
        V_region = dataset.get('V_region', 1e48)
        n_phonon = dataset.get('n_phonon', 1e26)
        ssq      = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)

        # Gaussian envelope
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))

        # Normalized Φ
        Phi = gaussian * S26

        # Vacuum density coupling
        rho_ratio = rho_scm / max(rho_ua, 1e-50)
        C_vac = rho_ratio * Phi * S26

        # Phonon pressure (thermodynamic)
        E_phonon = hbar * OMEGA_PHONON  # single phonon energy
        P_phonon = n_phonon * E_phonon * Phi / max(V_region, 1e-50)

        # Induced vacuum energy density
        rho_phonon_induced = P_phonon / c**2

        return {
            "C_vac": C_vac,
            "rho_ratio": rho_ratio,
            "Phi_normalized": Phi,
            "S26": S26,
            "gaussian": gaussian,
            "P_phonon_Pa": P_phonon,
            "E_single_phonon_J": E_phonon,
            "rho_phonon_induced_kg_m3": rho_phonon_induced,
            "equation": (
                "C_vac = (ρ_SCm/ρ_UA) · Φ_{1.25THz}(ω) · S₂₆([SSq])\n"
                f"      = {rho_ratio:.4e} × {Phi:.6e} × {S26:.6e}\n"
                f"      = {C_vac:.6e}\n\n"
                f"P_phonon = n · ℏω · Φ / V = {P_phonon:.6e} Pa"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  RESONANCE FREQUENCY SCAN
# ══════════════════════════════════════════════════════════════════════════════

class ResonanceFrequencyScan:
    """
    Scan the probe frequency ω across the 1.25 THz resonance and compute
    the response profile: a_res(ω), Φ(ω), C_vac(ω).

    Produces a Lorentzian-like profile (Gaussian envelope × S₂₆ factor)
    centered on ω_SCm = 2π × 1.25 THz.
    """

    def __init__(self):
        self._accel = ResonanceAcceleration()
        self._coupling = VacuumDensityCoupling()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          omega_min_THz: start frequency in THz (default 0.5)
          omega_max_THz: end frequency in THz (default 2.0)
          n_freq:    number of frequency points (default 100)
          Gamma:     linewidth (rad/s, default Γ_SCm)
        """
        omega_min = dataset.get('omega_min_THz', 0.5) * 2.0 * PI * 1e12
        omega_max = dataset.get('omega_max_THz', 2.0) * 2.0 * PI * 1e12
        n_freq    = dataset.get('n_freq', 100)

        scan = []
        for i in range(n_freq):
            frac = i / max(n_freq - 1, 1)
            omega_i = omega_min + frac * (omega_max - omega_min)
            freq_THz = omega_i / (2.0 * PI * 1e12)

            sub = dict(dataset)
            sub['omega'] = omega_i

            a_result = self._accel.compute(sub)
            c_result = self._coupling.compute(sub)

            scan.append({
                "freq_THz": round(freq_THz, 6),
                "omega_rad_s": omega_i,
                "a_res": a_result["a_res"],
                "Phi": a_result["Phi_125_THz"],
                "C_vac": c_result["C_vac"],
                "P_phonon": c_result["P_phonon_Pa"],
            })

        # FWHM estimation
        peak_a = max(scan, key=lambda x: abs(x["a_res"]))
        half_max = abs(peak_a["a_res"]) / 2.0
        above_half = [s for s in scan if abs(s["a_res"]) >= half_max]
        if len(above_half) >= 2:
            fwhm_THz = above_half[-1]["freq_THz"] - above_half[0]["freq_THz"]
        else:
            fwhm_THz = 0.0

        return {
            "scan": scan,
            "n_freq": n_freq,
            "freq_range_THz": [
                omega_min / (2.0 * PI * 1e12),
                omega_max / (2.0 * PI * 1e12),
            ],
            "peak_freq_THz": peak_a["freq_THz"],
            "peak_a_res": peak_a["a_res"],
            "FWHM_THz": round(fwhm_THz, 6),
            "description": (
                "Resonance frequency scan across 1.25 THz SCm phonon band.\n"
                f"Peak at {peak_a['freq_THz']:.4f} THz, "
                f"|a_res| = {abs(peak_a['a_res']):.6e}\n"
                f"FWHM ≈ {fwhm_THz:.4f} THz"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  PHONON DAMPING TIME-EVOLUTION
# ══════════════════════════════════════════════════════════════════════════════

class PhononDampingEvolution:
    """
    Time-domain evolution of the phonon resonance amplitude at fixed ω = ω_SCm:

        Φ(t) = Φ₀ · exp(-Γ t) · cos(ω_SCm t) · S₂₆([SSq])

    This gives the ring-down of the phonon mode after excitation,
    capturing the damping timescale τ = 1/Γ.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          t_max:     maximum time (s, default 10/Γ)
          n_steps:   time steps (default 200)
          Gamma:     linewidth (rad/s, default Γ_SCm)
          Phi_0:     initial fluence (default 1e20)
          ssq:       [SSq] parameter
        """
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        Phi_0   = dataset.get('Phi_0', PHI_0_DEFAULT)
        ssq     = dataset.get('ssq', SSQ)
        t_max   = dataset.get('t_max', 10.0 / max(Gamma, 1e-20))
        n_steps = dataset.get('n_steps', 200)

        S26 = S26_accelerated(ssq)
        tau = 1.0 / Gamma if Gamma > 0 else float('inf')

        evolution = []
        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            t = frac * t_max

            envelope = Phi_0 * math.exp(-Gamma * t) * S26
            oscillation = math.cos(OMEGA_PHONON * t)
            Phi_t = envelope * oscillation

            evolution.append({
                "t_s": t,
                "t_over_tau": t * Gamma,
                "Phi_t": Phi_t,
                "envelope": envelope,
                "oscillation": oscillation,
            })

        # Energy decay: E(t) ∝ Φ(t)² → E(t) = E₀ exp(-2Γt)
        E_half_life = math.log(2.0) / (2.0 * Gamma) if Gamma > 0 else float('inf')

        return {
            "evolution": evolution,
            "tau_damping_s": tau,
            "t_max_s": t_max,
            "n_steps": n_steps,
            "E_half_life_s": E_half_life,
            "Q_factor": OMEGA_PHONON / (2.0 * Gamma),
            "description": (
                "Phonon ring-down at ω_SCm = 1.25 THz.\n"
                f"Damping time τ = 1/Γ = {tau:.4e} s\n"
                f"Energy half-life = {E_half_life:.4e} s\n"
                f"Q = {OMEGA_PHONON / (2.0 * Gamma):.2f}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  MULTI-LAYER PHONON-GRAVITY COUPLING
# ══════════════════════════════════════════════════════════════════════════════

class MultiLayerPhononGravityCoupling:
    """
    Couple the phonon resonance into the 26-layer compressed gravity framework:

        g_phonon(r,t) = Σ_{i=1}^{26} [Ug1_i + Ug2_i + Ug3_i + Ug4_i] × a_res_i

    Each layer receives a layer-dependent Q_i modulation of the phonon:
        a_res_i = (F_{U,Bi}/F_U) · Φ(ω, Γ) · S₂₆ · Q_i

    where Q_i = 1/(1 + i·0.01) is the layer quantum factor.
    """

    def __init__(self):
        self._accel = ResonanceAcceleration()

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M:        central mass (kg)
          r:        radial distance (m)
          mu_s:     magnetic dipole moment (A·m²)
          omega_s:  spin angular frequency (rad/s)
          t:        time (s)
          omega:    phonon probe frequency (rad/s)
          Gamma:    phonon linewidth (rad/s)
          F_U_Bi, F_U: buoyancy/gravity
          ssq:      [SSq]
        """
        M       = dataset.get('M', 4e6 * M_sun)
        r       = dataset.get('r', 1e12)
        mu_s    = dataset.get('mu_s', 1e30)
        omega_s = dataset.get('omega_s', 1e-3)
        t       = dataset.get('t', 0.0)

        # Get base a_res
        a_result = self._accel.compute(dataset)
        a_res_base = a_result["a_res"]

        g_phonon_total = 0.0
        layers = []

        for layer in range(1, 27):
            Q_i = 1.0 / (1.0 + layer * 0.01)
            theta_i = 2.0 * PI * layer / 26.0
            proj = math.cos(theta_i)

            # Standard 4-channel gravity
            Ug1_i = G * mu_s**2 / max(r, 1e-10)**4 * Q_i * proj
            Ug2_i = G * mu_s * 1e46 / (max(r, 1e-10)**3 * c**2) * Q_i * proj
            Ug3_i = G * M / max(r, 1e-10)**2 * math.sin(omega_s * t) * Q_i * proj
            Ug4_i = G * M / max(r, 1e-10)**2 * SSQ * (RHO_SCM / max(RHO_UA, 1e-50)) * Q_i * proj

            Ug_i = Ug1_i + Ug2_i + Ug3_i + Ug4_i

            # Phonon-modulated gravity for this layer
            a_res_i = a_res_base * Q_i
            g_phonon_i = Ug_i * a_res_i

            g_phonon_total += g_phonon_i

            layers.append({
                "layer": layer,
                "Q_i": Q_i,
                "Ug_i": Ug_i,
                "a_res_i": a_res_i,
                "g_phonon_i": g_phonon_i,
            })

        return {
            "g_phonon_total": g_phonon_total,
            "a_res_base": a_res_base,
            "n_layers": 26,
            "layers": layers,
            "description": (
                "26-layer phonon-gravity coupling.\n"
                f"g_phonon(r,t) = Σ Ug_i · a_res_i = {g_phonon_total:.6e} m/s²\n"
                f"Base a_res = {a_res_base:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §7  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate SCm phonon resonance exploration."""
    print("=" * 72)
    print("SCm Phonon Resonance Exploration Engine — Session 211")
    print("=" * 72)

    # §1 Resonance acceleration
    print("\n── §1 Resonance Acceleration a_res ──")
    accel = ResonanceAcceleration()
    result = accel.compute({'F_U_Bi': 0.6, 'F_U': 1.0})
    print(f"  a_res = {result['a_res']:.6e}")
    print(f"  Regime: {result['regime']}")
    print(f"  Q = {result['Q_factor']:.2f}")

    # §2 Linewidth sweep
    print("\n── §2 Linewidth Γ Sweep ──")
    sweep = LinewidthGammaSweep()
    result = sweep.compute({'F_U_Bi': 0.6, 'F_U': 1.0})
    print(f"  Peak a_res = {result['peak_a_res']:.6e} at Γ = {result['peak_Gamma_THz']:.3f} THz")
    print(f"  Q at peak = {result['peak_Q']:.2f}")

    # §3 Vacuum coupling
    print("\n── §3 Vacuum Density Coupling ──")
    coupling = VacuumDensityCoupling()
    result = coupling.compute({})
    print(f"  C_vac = {result['C_vac']:.6e}")
    print(f"  P_phonon = {result['P_phonon_Pa']:.6e} Pa")

    # §4 Frequency scan
    print("\n── §4 Resonance Frequency Scan ──")
    scan = ResonanceFrequencyScan()
    result = scan.compute({'F_U_Bi': 0.6, 'F_U': 1.0, 'n_freq': 50})
    print(f"  Peak at {result['peak_freq_THz']:.4f} THz")
    print(f"  FWHM = {result['FWHM_THz']:.4f} THz")

    # §5 Damping evolution
    print("\n── §5 Phonon Damping Evolution ──")
    damping = PhononDampingEvolution()
    result = damping.compute({})
    print(f"  τ_damping = {result['tau_damping_s']:.4e} s")
    print(f"  E half-life = {result['E_half_life_s']:.4e} s")

    print(f"\n{'=' * 72}")
    print("SCm PHONON RESONANCE EXPLORATION COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
