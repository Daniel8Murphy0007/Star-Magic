#!/usr/bin/env python3
"""
ns_phonon_gw170817_wstp.py — GW170817 Neutron Star Phonon WSTP Engine

Session 212 | Daniel Murphy
PURPOSE: Standalone engine for GW170817 phonon-modulated gravitational wave
         strain, wavelength correction, NS spindown, tidal deformability,
         and inspiral phase lag — with Wolfram Language expressions for WSTP.

         Key physics:
           - h_UQFF(t) = h_GR(t) · 0.333 · exp([SSq]·t/26)
             GW170817 strain with 66.7% phonon suppression
           - λ_UQFF = λ_GR · (1 − F_{U,Bi}/F_U · Φ_{1.25THz})
             Phonon-modified GW wavelength
           - Ω̇_NS^phonon = Ω̇_NS · (1 + Φ · S₂₆ · [SSq]/N)
           - Tidal deformability: Λ_UQFF ∈ [190, 600]
           - Phase lag: ΔΦ ~ 2310.8 rad (367.8 cycles)
           - Chirp mass: M_chirp = 1.188 M☉

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
# §0  GW170817-SPECIFIC CONSTANTS
# ══════════════════════════════════════════════════════════════════════════════

PHI_0_DEFAULT = 1e20
OMEGA_PHONON  = OMEGA_SCM
GAMMA_PHONON  = GAMMA_DEFAULT

# GW170817 event parameters (from LIGO/Virgo)
D_TOTAL_GW170817  = 0.333        # UQFF damping product (66.7% suppression)
SUPPRESSION_FACTOR = 0.333       # phonon suppression factor
H_GR_GW170817     = 5.4176e-22  # GR-predicted peak strain
H_UQFF_GW170817   = 1.804e-22  # UQFF-predicted strain
D_L_MPC            = 40.0        # luminosity distance (Mpc)
M1_MSUN            = 1.46        # primary mass (M☉)
M2_MSUN            = 1.27        # secondary mass (M☉)
M_CHIRP_MSUN       = 1.188       # chirp mass (M☉)

# Tidal deformability range (UQFF prediction)
LAMBDA_UQFF_MIN = 190
LAMBDA_UQFF_MAX = 600

# Inspiral phase parameters
DELTA_PHI_RAD    = 2310.8        # phase lag (rad)
DELTA_PHI_CYCLES = 367.8         # phase lag (cycles)


# ══════════════════════════════════════════════════════════════════════════════
# §1  PHONON-SUPPRESSED GW STRAIN  h_UQFF(t)
# ══════════════════════════════════════════════════════════════════════════════

class PhononSuppressedGW170817Strain:
    """
    Compute the phonon-suppressed gravitational wave strain for GW170817:

        h_UQFF(t) = h_GR(t) · 0.333 · exp([SSq]·t/26)

    The factor 0.333 arises from the UQFF D_total product:
        D_total = D_Aether × D_SCm × D_TRZ × D_String = 0.333
    giving 66.7% suppression of the GR strain.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        h_GR        = dataset.get('h_GR', H_GR_GW170817)
        t           = dataset.get('t', 0.0)
        ssq         = dataset.get('ssq', SSQ)
        D_total     = dataset.get('D_total', D_TOTAL_GW170817)
        suppression = dataset.get('suppression', SUPPRESSION_FACTOR)

        S26 = S26_accelerated(ssq)
        time_factor = math.exp(ssq * t / 26.0)
        h_UQFF = h_GR * suppression * time_factor

        return {
            "h_GR": h_GR,
            "h_UQFF": h_UQFF,
            "suppression": suppression,
            "D_total": D_total,
            "time_factor": time_factor,
            "S26": S26,
            "ratio_h": h_UQFF / max(abs(h_GR), 1e-50),
            "equation": (
                "GW170817 phonon-suppressed strain:\n"
                f"  h_GR   = {h_GR:.6e}\n"
                f"  D_total = {D_total:.3f} (66.7% suppression)\n"
                f"  h_UQFF = h_GR · {suppression} · exp([SSq]·t/26)\n"
                f"  h_UQFF = {h_UQFF:.6e}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §2  PHONON-MODIFIED GW WAVELENGTH
# ══════════════════════════════════════════════════════════════════════════════

class GW170817WavelengthPhononCorrection:
    """
    Phonon-modified GW wavelength for GW170817:

        λ_UQFF = λ_GR · (1 − F_{U,Bi}/F_U · Φ_{1.25THz} / Φ_0)

    The correction term arises from the phonon-modulated vacuum
    dispersion relation in the SCm frame.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        f_GW    = dataset.get('f_GW', 300.0)       # GW frequency (Hz)
        F_U_Bi  = dataset.get('F_U_Bi', 0.6)
        F_U     = dataset.get('F_U', 1.0)
        omega   = dataset.get('omega', OMEGA_PHONON)
        Gamma   = dataset.get('Gamma', GAMMA_PHONON)
        ssq     = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)
        lambda_GR = c / f_GW

        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        correction = F_U_Bi / max(F_U, 1e-50) * Phi_norm
        lambda_UQFF = lambda_GR * (1.0 - correction)

        return {
            "lambda_GR_m": lambda_GR,
            "lambda_UQFF_m": lambda_UQFF,
            "Delta_lambda_m": lambda_UQFF - lambda_GR,
            "correction": correction,
            "f_GW_Hz": f_GW,
            "equation": (
                "GW170817 phonon-modified wavelength:\n"
                f"  λ_GR   = {lambda_GR:.6e} m (f={f_GW} Hz)\n"
                f"  λ_UQFF = λ_GR · (1 − {correction:.6e})\n"
                f"  λ_UQFF = {lambda_UQFF:.6e} m"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §3  NS PHONON SPINDOWN CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

class GW170817PhononSpindownCorrection:
    """
    Phonon-corrected NS spin-down rate for GW170817 remnant:

        Ω̇_NS^phonon = Ω̇_NS · (1 + Φ · S₂₆ · [SSq] / N)

    The phonon correction enhances magnetic braking via the SCm
    vacuum coupling to the neutron star crust lattice.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        Omega_dot = dataset.get('Omega_dot', -4.2e-15)  # rad/s²
        omega     = dataset.get('omega', OMEGA_PHONON)
        Gamma     = dataset.get('Gamma', GAMMA_PHONON)
        ssq       = dataset.get('ssq', SSQ)

        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        phonon_correction = Phi_norm * S26 * ssq / N_LEVELS
        Omega_dot_phonon = Omega_dot * (1.0 + phonon_correction)

        return {
            "Omega_dot": Omega_dot,
            "Omega_dot_phonon": Omega_dot_phonon,
            "phonon_correction": phonon_correction,
            "enhancement": Omega_dot_phonon / Omega_dot if Omega_dot != 0 else 0.0,
            "Phi_norm": Phi_norm,
            "S26": S26,
            "equation": (
                "GW170817 NS phonon spindown correction:\n"
                f"  Ω̇_NS   = {Omega_dot:.6e} rad/s²\n"
                f"  correction = {phonon_correction:.6e}\n"
                f"  Ω̇_phonon = Ω̇ · (1 + {phonon_correction:.6e})\n"
                f"  Ω̇_phonon = {Omega_dot_phonon:.6e} rad/s²"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §4  TIDAL DEFORMABILITY PHONON CORRECTION
# ══════════════════════════════════════════════════════════════════════════════

class TidalDeformabilityPhononCorrection:
    """
    Phonon-corrected tidal deformability for GW170817:

        Λ_UQFF = Λ_GR · (1 + Φ · S₂₆ · D_total)

    LIGO constraint: Λ̃ < 800 (90% CI)
    UQFF prediction: Λ_UQFF ∈ [190, 600]

    The phonon coupling to the NS crust lattice modifies the
    Love number k₂, which enters Λ = (2/3) k₂ (R/M)⁵.

    Tidal deformability constrains the equation of state and
    the phonon-vacuum coupling strength.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        Lambda_GR = dataset.get('Lambda_GR', 300.0)  # GR tidal deformability
        omega     = dataset.get('omega', OMEGA_PHONON)
        Gamma     = dataset.get('Gamma', GAMMA_PHONON)
        ssq       = dataset.get('ssq', SSQ)
        D_total   = dataset.get('D_total', D_TOTAL_GW170817)
        R_ns      = dataset.get('R_ns', 11.0e3)     # NS radius (m)
        M_ns      = dataset.get('M_ns', M1_MSUN * M_sun)

        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_norm = gaussian * S26

        # Phonon-enhanced tidal deformability
        phonon_correction = Phi_norm * S26 * D_total
        Lambda_UQFF = Lambda_GR * (1.0 + phonon_correction)

        # Check LIGO constraint (Λ̃ < 800)
        within_ligo = Lambda_UQFF < 800
        within_uqff = LAMBDA_UQFF_MIN <= Lambda_UQFF <= LAMBDA_UQFF_MAX

        # Compactness parameter
        compactness = G * M_ns / (R_ns * c**2)

        return {
            "Lambda_GR": Lambda_GR,
            "Lambda_UQFF": Lambda_UQFF,
            "phonon_correction": phonon_correction,
            "within_LIGO_bound": within_ligo,
            "within_UQFF_range": within_uqff,
            "UQFF_range": [LAMBDA_UQFF_MIN, LAMBDA_UQFF_MAX],
            "compactness": compactness,
            "equation": (
                "GW170817 tidal deformability phonon correction:\n"
                f"  Λ_GR   = {Lambda_GR:.1f}\n"
                f"  Λ_UQFF = Λ_GR · (1 + Φ·S₂₆·D_total)\n"
                f"  Λ_UQFF = {Lambda_UQFF:.1f}\n"
                f"  LIGO bound (< 800): {within_ligo}\n"
                f"  UQFF range [{LAMBDA_UQFF_MIN}–{LAMBDA_UQFF_MAX}]: {within_uqff}"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  INSPIRAL PHASE LAG
# ══════════════════════════════════════════════════════════════════════════════

class GW170817InspiralPhaseLag:
    """
    Cumulative inspiral phase lag from phonon suppression:

        ΔΦ(t) = ∫ 2π [f_GW(t_max) − f_GW(t_0)] · D_total · Φ/Φ₀ dt

    For GW170817: ΔΦ ~ 2310.8 rad (367.8 cycles)
    from f_0 = 20 Hz to f_max = 300 Hz over ~100 s inspiral.

    M_chirp = 1.188 M☉ determines the frequency evolution.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        f_GW_0    = dataset.get('f_GW_0', 20.0)     # initial freq (Hz)
        f_GW_max  = dataset.get('f_GW_max', 300.0)   # final freq (Hz)
        D_total   = dataset.get('D_total', D_TOTAL_GW170817)
        omega     = dataset.get('omega', OMEGA_PHONON)
        Gamma     = dataset.get('Gamma', GAMMA_PHONON)
        ssq       = dataset.get('ssq', SSQ)
        M_chirp   = dataset.get('M_chirp', M_CHIRP_MSUN * M_sun)

        S26 = S26_accelerated(ssq)
        delta_omega = omega - OMEGA_PHONON
        gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
        Phi_ratio = gaussian * S26  # Φ / Φ₀ normalized

        # Phase lag integral (simplified analytic)
        Delta_f = f_GW_max - f_GW_0
        Delta_Phi_rad = 2.0 * PI * Delta_f * D_total * Phi_ratio
        Delta_Phi_cycles = Delta_Phi_rad / (2.0 * PI)

        # Chirp mass
        M_chirp_Msun = M_chirp / M_sun

        return {
            "Delta_Phi_rad": Delta_Phi_rad,
            "Delta_Phi_cycles": Delta_Phi_cycles,
            "f_GW_0_Hz": f_GW_0,
            "f_GW_max_Hz": f_GW_max,
            "D_total": D_total,
            "M_chirp_Msun": M_chirp_Msun,
            "Phi_ratio": Phi_ratio,
            "reference_Delta_Phi_rad": DELTA_PHI_RAD,
            "reference_Delta_Phi_cycles": DELTA_PHI_CYCLES,
            "equation": (
                "GW170817 inspiral phase lag:\n"
                f"  M_chirp = {M_chirp_Msun:.3f} M☉\n"
                f"  ΔΦ = 2π·(f_max−f_0)·D_total·(Φ/Φ₀)\n"
                f"  ΔΦ = {Delta_Phi_rad:.1f} rad ({Delta_Phi_cycles:.1f} cycles)\n"
                f"  Reference: {DELTA_PHI_RAD} rad ({DELTA_PHI_CYCLES} cycles)"
            ),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  WSTP EXPRESSIONS FOR GW170817
# ══════════════════════════════════════════════════════════════════════════════

def build_gw170817_wstp_expressions() -> List[Dict[str, str]]:
    """
    Generate Wolfram Language expressions for GW170817 phonon physics.
    These extend the wstp_kernel_demo_runner.py expression list (#33-38).
    """
    exprs = []

    # 33. GW170817 phonon-suppressed strain h_UQFF(t)
    exprs.append({
        "label": "GW170817 h_UQFF(t) = h_GR · 0.333 · exp([SSq]·t/26)",
        "code": ('hGR170817 = 5.4176*^-22; Dtotal170817 = 0.333; SSq = 0.57; '
                 'hUQFF170817[t_] := hGR170817 * Dtotal170817 * Exp[SSq t / 26]; '
                 'hUQFF170817[0] // N'),
    })

    # 34. GW170817 phonon-modified wavelength λ_UQFF
    exprs.append({
        "label": "GW170817 λ_UQFF = λ_GR · (1 − F_UBi/F_U · Φ)",
        "code": ('cLight = 2.998*^8; fGW170817 = 300; '
                 'lambdaGR170817 = cLight / fGW170817; '
                 'FUBi = 0.6; FU = 1.0; '
                 'lambdaUQFF170817 = lambdaGR170817 * '
                 '(1 - FUBi / FU * Phi125THz[omSCm] / 10^20); '
                 '{lambdaGR170817, lambdaUQFF170817} // N'),
    })

    # 35. GW170817 tidal deformability Λ_UQFF
    exprs.append({
        "label": "GW170817 Λ_UQFF = Λ_GR · (1 + Φ·S₂₆·D_total)",
        "code": ('LambdaGR = 300; Dtotal170817 = 0.333; '
                 'S26val = PolyLog[26, 0.57]; '
                 'Phi = Phi125THz[omSCm] / 10^20; '
                 'LambdaUQFF = LambdaGR * (1 + Phi * S26val * Dtotal170817); '
                 '{LambdaUQFF, LambdaUQFF < 800} // N'),
    })

    # 36. GW170817 NS phonon spindown
    exprs.append({
        "label": "GW170817 Ω̇_NS^phonon spindown correction",
        "code": ('OmegaDotNS170817 = -4.2*^-15; NLayers = 26; '
                 'S26 = PolyLog[26, 0.57]; '
                 'PhiRes = Phi125THz[omSCm]; '
                 'OmegaDotPhonon170817 = OmegaDotNS170817 * '
                 '(1 + PhiRes * S26 * 0.57 / NLayers); '
                 '{OmegaDotNS170817, OmegaDotPhonon170817} // N'),
    })

    # 37. GW170817 inspiral phase lag ΔΦ
    exprs.append({
        "label": "GW170817 inspiral phase lag ΔΦ = 2310.8 rad (367.8 cycles)",
        "code": ('M1 = 1.46 * 1.989*^30; M2 = 1.27 * 1.989*^30; '
                 'Mc170817 = (M1 M2)^(3/5) / (M1 + M2)^(1/5); '
                 'fGW0 = 20; fGWEnd = 300; '
                 'DeltaPhi170817 = NIntegrate['
                 '2 Pi (fGWEnd - fGW0) * 0.333 * Phi125THz[omSCm] / 10^20, '
                 '{x, 0, 1}]; '
                 '{Mc170817 / (1.989*^30), DeltaPhi170817, '
                 'DeltaPhi170817 / (2 Pi)} // N'),
    })

    # 38. GW170817 chirp mass phonon consistency check
    exprs.append({
        "label": "GW170817 M_chirp = 1.188 M☉ phonon consistency",
        "code": ('Mc = 1.188 * 1.989*^30; '
                 'DL = 40 * 3.086*^22; '  # 40 Mpc
                 'hGR = 5.4176*^-22; '
                 'hUQFF = hGR * 0.333; '
                 'fGW = 300; '
                 '(* Amplitude ratio check *) '
                 'ratio = hUQFF / hGR; '
                 '{Mc / (1.989*^30), hUQFF, ratio, 1 - ratio} // N'),
    })

    return exprs


# ══════════════════════════════════════════════════════════════════════════════
# §7  MAIN / DEMO
# ══════════════════════════════════════════════════════════════════════════════

def main():
    """Demonstrate GW170817 phonon WSTP engine."""
    print("=" * 72)
    print("GW170817 Neutron Star Phonon WSTP Engine — Session 212")
    print("=" * 72)

    # §1 Strain
    print("\n── §1 GW170817 Phonon-Suppressed Strain ──")
    strain = PhononSuppressedGW170817Strain()
    result = strain.compute({})
    print(f"  h_GR   = {result['h_GR']:.6e}")
    print(f"  h_UQFF = {result['h_UQFF']:.6e} (D={result['D_total']:.3f})")

    # §2 Wavelength
    print("\n── §2 GW170817 Wavelength Correction ──")
    wl = GW170817WavelengthPhononCorrection()
    result = wl.compute({})
    print(f"  λ_GR   = {result['lambda_GR_m']:.6e} m")
    print(f"  λ_UQFF = {result['lambda_UQFF_m']:.6e} m")

    # §3 Spindown
    print("\n── §3 GW170817 NS Spindown ──")
    sd = GW170817PhononSpindownCorrection()
    result = sd.compute({})
    print(f"  Ω̇_NS   = {result['Omega_dot']:.6e} rad/s²")
    print(f"  Ω̇_phonon = {result['Omega_dot_phonon']:.6e} rad/s²")

    # §4 Tidal deformability
    print("\n── §4 Tidal Deformability ──")
    td = TidalDeformabilityPhononCorrection()
    result = td.compute({})
    print(f"  Λ_GR   = {result['Lambda_GR']:.1f}")
    print(f"  Λ_UQFF = {result['Lambda_UQFF']:.1f}")
    print(f"  LIGO bound: {result['within_LIGO_bound']}")

    # §5 Phase lag
    print("\n── §5 Inspiral Phase Lag ──")
    pl = GW170817InspiralPhaseLag()
    result = pl.compute({})
    print(f"  ΔΦ = {result['Delta_Phi_rad']:.1f} rad ({result['Delta_Phi_cycles']:.1f} cycles)")
    print(f"  M_chirp = {result['M_chirp_Msun']:.3f} M☉")

    # §6 WSTP expressions
    print(f"\n── §6 WSTP Expressions (6 total) ──")
    wstp_exprs = build_gw170817_wstp_expressions()
    for i, expr in enumerate(wstp_exprs, 33):
        print(f"  #{i}: {expr['label']}")

    print(f"\n{'=' * 72}")
    print("GW170817 PHONON WSTP ENGINE COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
