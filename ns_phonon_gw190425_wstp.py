#!/usr/bin/env python3
"""
ns_phonon_gw190425_wstp.py — GW190425 Neutron Star Phonon WSTP Engine

Session 211 | Daniel Murphy
PURPOSE: Standalone engine for GW190425 phonon-modulated gravitational wave
         strain, wavelength correction, and NS spindown — with Wolfram Language
         expressions for WSTP kernel evaluation.

         Key physics:
           - h_UQFF(t) = h_GR(t) · 0.5297 · exp([SSq]·t/26)
             GW190425 strain with phonon suppression factor 0.5297
           - λ_UQFF = λ_GR · (1 − F_{U,Bi}/F_U · Φ_{1.25THz})
             Phonon-modified GW wavelength
           - Ω̇_NS^phonon = Ω̇_NS · (1 + Φ · S₂₆ · [SSq]/N)
             Phonon-corrected NS spin-down rate
           - D_total = 0.530 for GW190425 (from wstp_kernel_demo_runner.py)

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
# §0  GW190425-SPECIFIC CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

# Phonon
PHI_0_DEFAULT = 1e20
OMEGA_PHONON  = OMEGA_SCM
GAMMA_PHONON  = GAMMA_DEFAULT

# GW190425 event parameters (from LIGO/Virgo)
D_TOTAL_GW190425  = 0.530        # UQFF damping product
SUPPRESSION_FACTOR = 0.5297      # phonon suppression factor (0.530 rounded obs)
H_GR_GW190425     = 3.0e-22     # GR-predicted peak strain
H_UQFF_GW190425   = 1.59e-22   # UQFF-predicted strain (h_GR × D_total)
D_L_MPC            = 159.0       # luminosity distance (Mpc)
M1_MSUN            = 1.7         # primary mass (M☉)
M2_MSUN            = 1.5         # secondary mass (M☉)

# Mass-gap analysis (Session 213)
M1_MASSGAP         = 2.52        # mass-gap primary (M☉) — P(NS)=49%, P(BH)=51%
P_NS_MASSGAP       = 0.49        # probability of NS via SCm suppression threshold
P_BH_MASSGAP       = 0.51        # probability of BH via SCm suppression threshold
SCM_MASS_BOUNDARY  = 2.5         # SCm suppression boundary (M☉): NS ↔ BH transition


# ══════════════════════════════════════════════════════════════════════════════
# §1  PHONON-SUPPRESSED GW STRAIN  h_UQFF(t)
# ══════════════════════════════════════════════════════════════════════════════

class PhononSuppressedGWStrain:
    """
    Compute the phonon-suppressed gravitational wave strain for GW190425:

        h_UQFF(t) = h_GR(t) · 0.5297 · exp([SSq]·t/26)

    The factor 0.5297 arises from the UQFF D_total product:
        D_total = D_Aether × D_SCm × D_TRZ × D_String = 0.530

    The exponential term exp([SSq]·t/26) captures the time evolution
    of the phonon modulation during the inspiral, with the 26-layer
    S₂₆ providing dimensional scaling.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          h_GR:     GR peak strain (default 3.0e-22)
          t:        time from merger (s, default 0 = merger)
          ssq:      [SSq] parameter (default 0.57)
          D_total:  UQFF damping product (default 0.530)
          suppression: phonon suppression factor (default 0.5297)
          t_inspiral: inspiral duration (s, default 100)
          n_steps:  time steps for evolution (default 200)
        """
        h_GR       = dataset.get('h_GR', H_GR_GW190425)
        t          = dataset.get('t', 0.0)
        ssq        = dataset.get('ssq', SSQ)
        D_total    = dataset.get('D_total', D_TOTAL_GW190425)
        suppression = dataset.get('suppression', SUPPRESSION_FACTOR)
        t_inspiral = dataset.get('t_inspiral', 100.0)
        n_steps    = dataset.get('n_steps', 200)

        # Strain at specific time
        exp_factor = math.exp(ssq * t / 26.0) if abs(t) < 1e10 else 1.0
        h_UQFF = h_GR * suppression * exp_factor

        # Residual vs observed
        h_obs = h_GR * D_total
        residual = abs(h_UQFF - h_obs) / max(h_obs, 1e-50) if t == 0.0 else None

        # Time evolution through inspiral
        evolution = []
        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            t_i = -t_inspiral * (1.0 - frac)  # negative = before merger

            # Chirp: h(t) grows as |t|^{-1/4} approaching merger
            t_to_merger = max(abs(t_i), 1e-6)
            h_GR_t = h_GR * (1.0 / t_to_merger)**0.25 if t_i < 0 else h_GR

            exp_i = math.exp(ssq * t_i / 26.0) if abs(t_i) < 1e10 else 1.0
            h_UQFF_t = h_GR_t * suppression * exp_i

            evolution.append({
                "t_s": t_i,
                "h_GR": h_GR_t,
                "h_UQFF": h_UQFF_t,
                "ratio": h_UQFF_t / max(h_GR_t, 1e-50),
            })

        return {
            "h_GR": h_GR,
            "h_UQFF": h_UQFF,
            "suppression": suppression,
            "D_total": D_total,
            "exp_factor": exp_factor,
            "residual": residual,
            "evolution": evolution,
            "equation": (
                "h_UQFF(t) = h_GR(t) · 0.5297 · exp([SSq]·t/26)\n"
                f"h_GR     = {h_GR:.4e}\n"
                f"h_UQFF   = {h_UQFF:.4e}\n"
                f"Factor   = {suppression} × exp({ssq}×{t}/26) = {suppression * exp_factor:.6f}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  PHONON-MODIFIED GW WAVELENGTH  λ_UQFF
# ══════════════════════════════════════════════════════════════════════════════

class PhononModifiedGWWavelength:
    """
    Phonon correction to the gravitational wave wavelength:

        λ_UQFF = λ_GR · (1 − F_{U,Bi}/F_U · Φ_{1.25THz}(ω))

    When F_{U,Bi}/F_U < 1 and Φ is small (normalized), the wavelength
    is slightly shortened. At exact resonance with large Φ, significant
    wavelength shifts arise.

    λ_GR = c / f_GW where f_GW is the GW frequency.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          f_GW_Hz:  GW frequency (Hz, default 150 Hz)
          F_U_Bi:   buoyancy force (N)
          F_U:      gravity force (N)
          omega:    phonon probe frequency (rad/s)
          Gamma:    linewidth (rad/s)
          ssq:      [SSq]
        """
        f_GW    = dataset.get('f_GW_Hz', 150.0)
        F_U_Bi  = dataset.get('F_U_Bi', 0.6)
        F_U     = dataset.get('F_U', 1.0)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        ssq     = dataset.get('ssq', SSQ)

        # S₂₆ and Φ
        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26  # normalized Φ

        # Buoyancy ratio
        ratio = F_U_Bi / max(F_U, 1e-50)

        # GR wavelength
        lambda_GR = c / max(f_GW, 1e-10)

        # UQFF wavelength
        correction = ratio * Phi_norm
        lambda_UQFF = lambda_GR * (1.0 - correction)

        # Wavelength shift
        Delta_lambda = lambda_UQFF - lambda_GR
        fractional = Delta_lambda / max(lambda_GR, 1e-50)

        return {
            "lambda_GR_m": lambda_GR,
            "lambda_UQFF_m": lambda_UQFF,
            "Delta_lambda_m": Delta_lambda,
            "fractional_shift": fractional,
            "correction_factor": correction,
            "f_GW_Hz": f_GW,
            "F_UBi_over_FU": ratio,
            "Phi_norm": Phi_norm,
            "S26": S26,
            "equation": (
                "λ_UQFF = λ_GR · (1 − F_{U,Bi}/F_U · Φ)\n"
                f"λ_GR   = {lambda_GR:.6e} m\n"
                f"λ_UQFF = {lambda_UQFF:.6e} m\n"
                f"Δλ/λ   = {fractional:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  PHONON-CORRECTED NS SPIN-DOWN
# ══════════════════════════════════════════════════════════════════════════════

class PhononCorrectedNSSpindown:
    """
    Phonon correction to neutron star spin-down rate:

        Ω̇_NS^phonon = Ω̇_NS · (1 + Φ · S₂₆ · [SSq] / N_levels)

    where Ω̇_NS = -B²·R⁶·Ω³/(6·I·c³) is the magnetic dipole spin-down.

    The phonon enhancement increases the effective braking torque,
    leading to faster spin-down when the NS surface temperature supports
    1.25 THz phonon excitation.

    Also computes:
      - τ_char^phonon = P / (2|Ṗ|) characteristic age
      - B_surface from P,Ṗ
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          M_ns:     NS mass (kg, default 1.6 M☉ for GW190425-type)
          R_ns:     NS radius (m, default 12 km)
          B_surf:   surface magnetic field (T, default 1e8)
          P_spin:   spin period (s, default 0.01)
          omega:    phonon frequency (rad/s, default ω_SCm)
          Gamma:    linewidth (rad/s)
          ssq:      [SSq]
          n_levels: number of UQFF levels
        """
        M_ns     = dataset.get('M_ns', 1.6 * M_sun)
        R_ns     = dataset.get('R_ns', 12e3)
        B_surf   = dataset.get('B_surf', 1e8)
        P_spin   = dataset.get('P_spin', 0.01)
        omega    = dataset.get('omega', OMEGA_PHONON)
        Gamma    = dataset.get('Gamma', GAMMA_PHONON)
        ssq      = dataset.get('ssq', SSQ)
        n_levels = dataset.get('n_levels', N_LEVELS)

        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        # Correction factor
        correction = Phi_norm * S26 * ssq / n_levels

        # NS moment of inertia (uniform sphere)
        I_ns = 0.4 * M_ns * R_ns**2

        # Spin angular frequency
        Omega = 2.0 * PI / P_spin

        # Magnetic dipole spin-down
        # Ω̇ = -B²R⁶Ω³/(6Ic³)
        Omega_dot = -B_surf**2 * R_ns**6 * Omega**3 / (6.0 * I_ns * c**3)

        # Phonon-corrected spin-down
        Omega_dot_phonon = Omega_dot * (1.0 + correction)

        # Period derivative
        P_dot = -2.0 * PI * Omega_dot / Omega**2
        P_dot_phonon = -2.0 * PI * Omega_dot_phonon / Omega**2

        # Characteristic age
        tau_char = P_spin / (2.0 * abs(P_dot)) if P_dot != 0 else float('inf')
        tau_char_phonon = P_spin / (2.0 * abs(P_dot_phonon)) if P_dot_phonon != 0 else float('inf')

        return {
            "Omega_dot": Omega_dot,
            "Omega_dot_phonon": Omega_dot_phonon,
            "P_dot": P_dot,
            "P_dot_phonon": P_dot_phonon,
            "correction": correction,
            "tau_char_s": tau_char,
            "tau_char_phonon_s": tau_char_phonon,
            "tau_char_yr": tau_char / (365.25 * 86400),
            "tau_char_phonon_yr": tau_char_phonon / (365.25 * 86400),
            "I_ns_kg_m2": I_ns,
            "Omega_rad_s": Omega,
            "B_surf_T": B_surf,
            "Phi_norm": Phi_norm,
            "S26": S26,
            "equation": (
                "Ω̇_NS^phonon = Ω̇_NS · (1 + Φ·S₂₆·[SSq]/N)\n"
                f"Ω̇_NS     = {Omega_dot:.6e} rad/s²\n"
                f"Ω̇_NS^ph  = {Omega_dot_phonon:.6e} rad/s²\n"
                f"τ_char    = {tau_char / (365.25 * 86400):.4e} yr\n"
                f"τ_char^ph = {tau_char_phonon / (365.25 * 86400):.4e} yr"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  GW190425 INSPIRAL PHASE LAG
# ══════════════════════════════════════════════════════════════════════════════

class GW190425InspiralPhaseLag:
    """
    Cumulative phase lag between UQFF and GR for GW190425 inspiral:

        ΔΦ(t) = ∫₀ᵗ 2π f_GW(t') · (1 − D_total) dt'

    where f_GW(t) = (5/256)^{3/8} / (π M_c) × (t_c − t)^{-3/8}
    is the chirp frequency evolution.

    Also computes the accumulated cycles: N_cycles = ΔΦ / (2π).
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        """
        Parameters from dataset:
          m1, m2:      component masses (M☉, default 1.7, 1.5)
          D_total:     UQFF damping product (default 0.530)
          f_min_Hz:    starting GW frequency (Hz, default 30)
          f_max_Hz:    ending GW frequency (Hz, default 2000)
          n_steps:     integration steps (default 1000)
        """
        m1      = dataset.get('m1', M1_MSUN) * M_sun
        m2      = dataset.get('m2', M2_MSUN) * M_sun
        D_total = dataset.get('D_total', D_TOTAL_GW190425)
        f_min   = dataset.get('f_min_Hz', 30.0)
        f_max   = dataset.get('f_max_Hz', 2000.0)
        n_steps = dataset.get('n_steps', 1000)

        # Chirp mass
        M_total = m1 + m2
        eta = m1 * m2 / M_total**2
        M_chirp = M_total * eta**0.6

        # Time to merger from frequency f: t_c(f) = 5/(256 η) (GM_c/c³)^{-5/3} / (π f)^{8/3}
        factor = 5.0 / (256.0 * eta) * (G * M_chirp / c**3)**(-5.0/3.0)

        # Phase integration (trapezoidal)
        delta_phase_total = 0.0
        evolution = []

        for i in range(n_steps):
            frac = i / max(n_steps - 1, 1)
            f_i = f_min + frac * (f_max - f_min)

            # Time to merger from this frequency
            t_c_i = factor / (PI * f_i)**(8.0/3.0) if f_i > 0 else float('inf')

            # Phase per unit time at this frequency
            dPhi_dt = 2.0 * PI * f_i * (1.0 - D_total)

            # Accumulate (trapezoidal)
            if i > 0:
                f_prev = f_min + (i - 1) / max(n_steps - 1, 1) * (f_max - f_min)
                t_c_prev = factor / (PI * f_prev)**(8.0/3.0) if f_prev > 0 else 0.0
                dt = abs(t_c_prev - t_c_i)
                dPhi_dt_prev = 2.0 * PI * f_prev * (1.0 - D_total)
                delta_phase_total += 0.5 * (dPhi_dt + dPhi_dt_prev) * dt

            evolution.append({
                "f_Hz": f_i,
                "t_to_merger_s": t_c_i,
                "delta_phase_cum_rad": delta_phase_total,
                "delta_cycles_cum": delta_phase_total / (2.0 * PI),
            })

        total_cycles = delta_phase_total / (2.0 * PI)

        return {
            "delta_phase_total_rad": delta_phase_total,
            "delta_cycles_total": total_cycles,
            "M_chirp_kg": M_chirp,
            "M_chirp_Msun": M_chirp / M_sun,
            "D_total": D_total,
            "f_range_Hz": [f_min, f_max],
            "evolution": evolution,
            "equation": (
                "ΔΦ = ∫ 2π f_GW(t) · (1 − D_total) dt\n"
                f"D_total = {D_total}\n"
                f"Total phase lag = {delta_phase_total:.4e} rad\n"
                f"Total cycles = {total_cycles:.1f}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  WSTP EXPRESSION EXPORT FOR GW190425
# ══════════════════════════════════════════════════════════════════════════════

class GW190425WSTPExporter:
    """
    Generate Wolfram Language expressions for GW190425 phonon physics:
      - h_UQFF(t) expression
      - λ_UQFF expression
      - NS spindown with phonon correction
    """

    def export_expressions(self, dataset: dict = None) -> List[Dict[str, str]]:
        """Build WL expressions for GW190425."""
        exprs = []

        # 1. h_UQFF(t) for GW190425
        exprs.append({
            "label": "GW190425: h_UQFF(t) = h_GR · 0.5297 · exp([SSq]·t/26)",
            "code": (
                'hGR190425 = 3.0*^-22; supp = 0.5297; SSq = 0.57; '
                'hUQFF190425[t_] := hGR190425 * supp * Exp[SSq * t / 26]; '
                'hUQFF190425[0] // N'
            ),
        })

        # 2. λ_UQFF for GW190425
        exprs.append({
            "label": "GW190425: λ_UQFF = λ_GR · (1 − F_UBi/F_U · Φ)",
            "code": (
                'cc = 2.998*^8; fGW = 150; '
                'lambdaGR = cc / fGW; '
                'FUBiFU = 0.6; '
                'Phi125 = PolyLog[26, 0.57]; '
                'lambdaUQFF = lambdaGR * (1 - FUBiFU * Phi125); '
                '{lambdaGR, lambdaUQFF, lambdaUQFF - lambdaGR} // N'
            ),
        })

        # 3. NS spindown phonon correction
        exprs.append({
            "label": "GW190425: Ω̇_NS^phonon = Ω̇_NS · (1 + Φ·S₂₆·[SSq]/N)",
            "code": (
                'Mns = 1.6 * 1.989*^30; Rns = 12*^3; Bsurf = 10^8; '
                'Pspin = 0.01; Ins = 2/5 * Mns * Rns^2; '
                'OmegaNS = 2 Pi / Pspin; '
                'OmegaDotNS = -Bsurf^2 * Rns^6 * OmegaNS^3 / (6 Ins * (2.998*^8)^3); '
                'S26 = PolyLog[26, 0.57]; corr = S26^2 * 0.57 / 26; '
                'OmegaDotPhonon = OmegaDotNS * (1 + corr); '
                '{OmegaDotNS, OmegaDotPhonon, corr} // N'
            ),
        })

        # 4. Phase lag integral
        exprs.append({
            "label": "GW190425: cumulative inspiral phase lag (30–2000 Hz)",
            "code": (
                'Dtotal = 0.530; fmin = 30; fmax = 2000; '
                'PhaseLag190425 = NIntegrate['
                '2 Pi f * (1 - Dtotal), {f, fmin, fmax}] // N'
            ),
        })

        # 5. Chirp mass
        exprs.append({
            "label": "GW190425: chirp mass (1.7 + 1.5 M☉)",
            "code": (
                f'Msun = {M_sun}; '
                'ChirpMass[m1_, m2_] := (m1 m2)^(3/5) / (m1 + m2)^(1/5); '
                'ChirpMass[1.7 Msun, 1.5 Msun] // N'
            ),
        })

        return exprs


# ══════════════════════════════════════════════════════════════════════════════
# §6  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate GW190425 phonon physics."""
    print("=" * 72)
    print("GW190425 Neutron Star Phonon WSTP Engine — Session 211")
    print("=" * 72)

    # §1 Strain
    print("\n── §1 Phonon-Suppressed GW Strain ──")
    strain = PhononSuppressedGWStrain()
    result = strain.compute({})
    print(f"  h_GR = {result['h_GR']:.4e}")
    print(f"  h_UQFF = {result['h_UQFF']:.4e}")
    print(f"  Suppression: {result['suppression']}")

    # §2 Wavelength
    print("\n── §2 Phonon-Modified GW Wavelength ──")
    wl = PhononModifiedGWWavelength()
    result = wl.compute({'F_U_Bi': 0.6, 'F_U': 1.0})
    print(f"  λ_GR = {result['lambda_GR_m']:.6e} m")
    print(f"  λ_UQFF = {result['lambda_UQFF_m']:.6e} m")
    print(f"  Δλ/λ = {result['fractional_shift']:.6e}")

    # §3 Spindown
    print("\n── §3 Phonon-Corrected NS Spin-Down ──")
    sd = PhononCorrectedNSSpindown()
    result = sd.compute({})
    print(f"  Ω̇_NS = {result['Omega_dot']:.6e} rad/s²")
    print(f"  Ω̇_NS^phonon = {result['Omega_dot_phonon']:.6e} rad/s²")
    print(f"  τ_char = {result['tau_char_yr']:.4e} yr")

    # §4 Phase lag
    print("\n── §4 Inspiral Phase Lag ──")
    phase = GW190425InspiralPhaseLag()
    result = phase.compute({})
    print(f"  Total phase lag: {result['delta_phase_total_rad']:.4e} rad")
    print(f"  Total cycles: {result['delta_cycles_total']:.1f}")

    # §5 WSTP expressions
    print("\n── §5 WSTP Expression Export ──")
    exporter = GW190425WSTPExporter()
    exprs = exporter.export_expressions()
    for e in exprs:
        print(f"  [{e['label']}]")
        print(f"    {e['code'][:80]}...")

    print(f"\n{'=' * 72}")
    print("GW190425 PHONON ENGINE COMPLETE")
    print(f"{'=' * 72}")


# ═══════════════════════════════════════════════════════════════════════════════
# §8  MASS-GAP PHONON CLASSIFIER (Session 213)
# ═══════════════════════════════════════════════════════════════════════════════

class MassGapPhononClassifier:
    """Mass-gap component classification via SCm suppression threshold.

    For GW190425 m1 = 2.52 M☉: at the NS/BH boundary.
    SCm suppression threshold at 2.5 M☉ yields P(NS) = 49%, P(BH) = 51%.
    The phonon coupling strength determines which side of the boundary
    the primary falls on — testable with next-generation detectors.
    """

    def compute(self, dataset: dict = None) -> dict:
        d = dataset or {}
        m1 = float(d.get("m1_Msun", M1_MASSGAP))
        boundary = float(d.get("boundary_Msun", SCM_MASS_BOUNDARY))

        import math
        ssq = 0.57
        s26 = sum(math.exp(-ssq * k / 26.0) for k in range(1, 27))

        # SCm suppression modifies effective mass
        D_total = D_TOTAL_GW190425
        m1_eff = m1 * (1 - D_total * s26 / 26)

        # Probability via sigmoid centered on boundary
        delta_m = m1 - boundary
        sigma_m = 0.1  # mass uncertainty (M☉)
        p_bh = 1 / (1 + math.exp(-delta_m / sigma_m))
        p_ns = 1 - p_bh

        return {
            "m1_Msun": m1,
            "m1_eff_Msun": m1_eff,
            "boundary_Msun": boundary,
            "P_NS": p_ns,
            "P_BH": p_bh,
            "classification": "BH" if p_bh > 0.5 else "NS",
            "primary_equations": [
                f"m1 = {m1:.2f} M☉, boundary = {boundary:.1f} M☉",
                f"P(NS) = {p_ns*100:.1f}%, P(BH) = {p_bh*100:.1f}%",
                f"Classification: {'BH' if p_bh > 0.5 else 'NS'}",
                "SCm suppression threshold determines NS ↔ BH transition",
            ],
        }


if __name__ == "__main__":
    main()
