"""
Session 264 - PAPER_1171 fudge-factor audit and replacement
============================================================

PAPER_1171 attempts a first-principles KK regulator derivation:

    rho_KK = (3*zeta(5)/(64*pi^6)) * (D_crit/D_BSFG)^4 * (v_UA/c)^4
              * rho_SCm * c^2 * (c/v_UA)^2 * 10^17

It claims to evaluate to 5.951e-10 J/m^3 with 0.15% residual.  The
paper labels the trailing factor "10^17" as a "dimensional bookkeeping
factor" but the same paper's table gives the bookkeeping value as
"1.0e+8" — these two statements are inconsistent.  Computing the
formula numerically:

    4.857e-5 * 31.605 * 1.235e-2 * 5.748e-19  =  1.090e-23 J/m^3
    -> Multiplied by 10^8:  1.090e-15 J/m^3
    -> Multiplied by 10^17: 1.090e-6  J/m^3
    -> Target was       :   5.95e-10  J/m^3

Neither bookkeeping value reproduces the target.  The "first-principles"
derivation is in fact a fit with a tuneable power-of-10 conversion.

The honest factor needed to make the formula match would be 10^13.74,
which is not a structural primitive — confirming this is a fudge.

REPLACEMENT (carried from Session 263)
======================================
Cosmological constant closes directly from 6 anchors:

    rho_Lambda_obs = 2 * D_BSFG * D_phys * rho_A * f_THz
                   = 6.0e-10 J/m^3   (residual 0.67%)

No KK tower required.  PAPER_1171's (13/3)^4 prefactor IS structurally
meaningful (it appears via L_KK* = (D_BSFG/D_crit)(c/v_UA) = 9/13
canonical), but the way PAPER_1171 plugs it back through QFT zeta-
regularisation does not produce a clean SI number without the 10^17
fudge.

CLEAN m_1 INTERPRETATION
========================
The (13/3) factor IS the ratio D_crit/D_BSFG.  Reading m_1 as a
characteristic energy:

    m_1 c^2 = (13/3) * v_UA^2 / c * c^2 = (13/3) * v_UA^2 * c
            = (13/3) * (1e8)^2 * 3e8 = 13e24 J

This is dimensionally inconsistent — v_UA^2 * c has units m^3/s^3,
not energy.  PAPER_1171 has a units error in equation 3.

CORRECT reading: m_1 should be set by h*f_THz, not v_UA^2/c.

    m_1_c2_correct = (D_crit/D_BSFG) * h * f_THz
                   = (13/3) * 6.626e-34 * 1.25e12
                   = 3.59e-21 J
                   = 22.4 meV

This is the characteristic resonant-pair excitation energy at the
cosmological-boundary fixed point — comparable to the 7.49 meV figure
PAPER_1171 quotes but with a clean structural derivation.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

# 6 anchors + SI scale anchors
RHO_UA = 7.09e-36
RHO_SCm = 7.09e-37
RHO_A = 1.0e-23
V_SCm = 1.0e8
f_THz = 1.25e12

c = 2.99792458e8
h = 6.62607015e-34
hbar = 1.054571817e-34

# Structural primitives
D_crit = 26.0
D_BSFG = 6.0
D_phys = 4.0

RHO_LAMBDA_OBS = 5.96e-10


def main():
    print("=" * 76)
    print("SESSION 264  -  PAPER_1171 fudge-factor audit and replacement")
    print("=" * 76)
    print()

    print("PART 1: Compute PAPER_1171 formula exactly as written")
    print("-" * 76)
    zeta_5 = 1.0369277551433699
    prefactor = 3.0 * zeta_5 / (64.0 * math.pi ** 6)
    ratio4 = (D_crit / D_BSFG) ** 4    # = (13/3)^4
    v_ratio4 = (V_SCm / c) ** 4         # (v_UA/c)^4
    ref_scale = RHO_SCm * c ** 2 * (c / V_SCm) ** 2

    bare = prefactor * ratio4 * v_ratio4 * ref_scale
    print(f"  Prefactor 3 zeta(5)/(64 pi^6)     = {prefactor:.4e}")
    print(f"  (D_crit/D_BSFG)^4 = (13/3)^4      = {ratio4:.4f}")
    print(f"  (v_UA/c)^4                        = {v_ratio4:.4e}")
    print(f"  rho_SCm c^2 (c/v_UA)^2            = {ref_scale:.4e}")
    print(f"  Bare product (no 10^N factor)     = {bare:.4e} J/m^3")
    print()
    print(f"  With 10^8 bookkeeping:  {bare*1e8:.4e}")
    print(f"  With 10^17 bookkeeping: {bare*1e17:.4e}")
    print(f"  Target (paper claim):   5.951e-10")
    print()
    needed_log10 = math.log10(5.951e-10 / bare)
    print(f"  Actual factor needed:   10^{needed_log10:.2f}")
    print(f"  Conclusion: neither 10^8 nor 10^17 reproduces the claim.")
    print(f"  The required factor 10^{needed_log10:.2f} is not a structural")
    print(f"  primitive — this is a hand-tuned fudge.")
    print()

    print("PART 2: Units check on PAPER_1171 equation 3 (m_1 definition)")
    print("-" * 76)
    print()
    print("  Paper writes: m_1 = v_UA / L_KK*  with  L_KK* = (D_BSFG/D_crit)*(c/v_UA)")
    print("  m_1 = v_UA / ((D_BSFG/D_crit)(c/v_UA))")
    print("       = (D_crit/D_BSFG) * v_UA^2 / c        (units of m/s)")
    print()
    print("  But m_1 is supposed to be a MASS in QFT KK calculations.  v_UA^2/c")
    print("  has units of m/s (velocity), not kg.  This is dimensionally wrong.")
    print()
    print("  The paper later treats m_1^4 as having units of J^4/m^4 via")
    print("  rho_SCm * (some factors), but the unit-tracking is inconsistent.")
    print()

    print("PART 3: Clean re-interpretation of m_1")
    print("-" * 76)
    print()
    m1_c2_clean = (D_crit / D_BSFG) * h * f_THz
    print(f"  m_1 c^2 = (D_crit/D_BSFG) * h * f_THz")
    print(f"          = (13/3) * 6.626e-34 * 1.25e12")
    print(f"          = {m1_c2_clean:.4e} J")
    print(f"          = {m1_c2_clean/1.602e-22:.2f} meV")
    print()
    print(f"  Compare paper's m_1 c^2 ~ 7.49 meV (from CondensedPhysics4 line 15955).")
    print(f"  Ratio of clean/paper = {m1_c2_clean/(7.49*1.602e-22):.3f}")
    print(f"  These differ by ~3x; both are O(10 meV).  The clean form is fully")
    print(f"  structural (D_crit/D_BSFG * h * f_THz) and dimensionally consistent.")
    print()

    print("PART 4: Carry forward Session 263 closure as canonical")
    print("-" * 76)
    rho_lambda_S263 = 2.0 * D_BSFG * D_phys * RHO_A * f_THz
    residual = abs(rho_lambda_S263 - RHO_LAMBDA_OBS) / RHO_LAMBDA_OBS * 100.0
    print(f"  S263 canonical: rho_Lambda = 2 D_BSFG D_phys rho_A f_THz")
    print(f"                             = {rho_lambda_S263:.4e} J/m^3")
    print(f"                  residual = {residual:.3f}%   [CLOSED]")
    print()
    print(f"  This supersedes both PAPER_1170 KK tower and PAPER_1171's")
    print(f"  (13/3)^4 + 10^17 derivation.")
    print()

    print("PART 5: Recommended PAPER_1171 revision")
    print("-" * 76)
    print()
    print("  REMOVE:  Sections 3-5 KK-tower zeta-regularisation derivation")
    print("           (units inconsistent, requires 10^13.74 fudge).")
    print()
    print("  REPLACE WITH: Session 263 closure")
    print("    rho_Lambda = 2 * D_BSFG * D_phys * rho_A * f_THz")
    print("    = (paired-buoyancy crossing factor 2) x")
    print("      (BSFG-physical-spacetime cell volume D_BSFG * D_phys = 24) x")
    print("      (ambient aether anchor density rho_A) x")
    print("      (THz reactant-pair oscillation frequency f_THz)")
    print()
    print("  RETAIN: PAPER_1171's identification (13/3) = D_crit/D_BSFG as a")
    print("          structural ratio appearing throughout UQFF.  Used at")
    print("          stage-3 Chandrasekhar bound (Session 262) and Riemann")
    print("          critical-line argument (Session 261, Level_13/D_crit=1/2).")
    print()

    out = {
        "session": 264,
        "paper_1171_audit": {
            "issue": "10^17 'bookkeeping factor' inconsistent with paper's own"
                     " table value 10^8; neither value matches numerical claim."
                     "  Actual factor required: 10^{:.2f}".format(needed_log10),
            "units_error": "m_1 = v_UA^2 / c * (D_crit/D_BSFG) has units of"
                           " m/s (velocity), not mass/energy.  Paper's"
                           " subsequent treatment is dimensionally inconsistent.",
            "bare_formula_result": bare,
            "target_value": 5.951e-10,
            "required_log10_factor": needed_log10,
            "verdict": "FUDGE - hand-tuned conversion factor.",
        },
        "replacement": {
            "form": "rho_Lambda = 2 * D_BSFG * D_phys * rho_A * f_THz",
            "value": rho_lambda_S263,
            "observed": RHO_LAMBDA_OBS,
            "residual_pct": residual,
            "verdict": "CLOSED",
            "source_session": 263,
        },
        "structural_clean_m1": {
            "form": "m_1 c^2 = (D_crit/D_BSFG) * h * f_THz",
            "value_J": m1_c2_clean,
            "value_meV": m1_c2_clean / 1.602e-22,
            "interpretation": "Resonant-pair excitation energy at THz scale"
                              " with structural enhancement ratio 13/3.",
        },
    }
    Path("_session264_closures.json").write_text(json.dumps(out, indent=2))
    print()
    print("Wrote _session264_closures.json")


if __name__ == "__main__":
    main()
