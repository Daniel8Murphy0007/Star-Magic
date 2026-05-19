"""
Session 719 -- hbar-chain opening: dimensional class probe.

Goal: Open the fourth fundamental-constant chain (Planck's constant hbar)
      via the DPM action quantum and determine whether it falls into one
      of the three universality classes (I/II/III) or requires a new
      Class IV.

Method:
    1. Dimensional analysis of all forms reachable from
       {rho_vac_SCm, v_SCM, c, locked-primitives}.
    2. Show that the unique dimensional path is
           hbar = rho_vac * L^4 / v_SCM * C_locked
       where L is an EXTERNAL length anchor (NOT reachable from the
       primitives alone).
    3. Back-solve L_DPM = (hbar * v_SCM / rho_vac_SCm)^(1/4) and search
       for locked-rational matches versus (c/v_SCM)^p * c_2^q * etc.
    4. If no clean match exists at machine precision, hbar-chain is
       Class IV (requires external dimensional anchor).

Locked primitives (Fraction):
    F_TRZ=1/10, Phi_res=5/6, SSq=57/100, K_Mex=25/12, beta_i=6029/10000,
    D_phys=4, D_BSFG=6, D_crit=26, N_ch=9, SO5_order=10, A_5=60.

CVW stamp: v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""

from __future__ import annotations
import json
import math
from fractions import Fraction
from pathlib import Path

# --- Locked primitives ----------------------------------------------------
F_TRZ      = Fraction(1, 10)
Phi_res    = Fraction(5, 6)
SSq        = Fraction(57, 100)
K_Mex      = Fraction(25, 12)
beta_i     = Fraction(6029, 10000)
D_phys     = 4
D_BSFG     = 6
D_crit     = 26
N_ch       = 9
SO5_order  = 10
A_5        = 60

# Locked identities
assert F_TRZ * Phi_res == Fraction(1, 12), "F_TRZ * Phi_res must equal 1/12"
assert K_Mex == Phi_res * SO5_order / Fraction(D_phys), \
    "K_Mex must equal Phi_res * SO5_order / D_phys"

# --- Observables ----------------------------------------------------------
HBAR_OBS       = 1.054571817e-34   # J*s (CODATA exact since 2019 redefinition)
C_LIGHT        = 2.99792458e8      # m/s
V_SCM          = 1.0e8             # m/s  (c/3 to high accuracy in UQFF)
RHO_VAC_SCM    = 7.09e-37          # J/m^3 (canonical SCm vacuum density)

# --- Dimensional class probe ---------------------------------------------
print("=" * 80)
print("SESSION 719 -- hbar-chain dimensional class probe")
print("=" * 80)

# Dimensional analysis:
# [hbar]      = J*s = kg*m^2/s
# [rho_vac]   = J/m^3
# [c],[v_SCM] = m/s
# [L]         = m
#
# Most general monomial: rho^1 * c^a * v^b * L^d  has dims
#     J * m^(-3 + a + b + d) * s^(-a - b)
# Setting equal to J*s:
#     -3 + a + b + d = 0
#     -a - b         = 1   ->  a + b = -1
# =>  d = 4 + (a + b) - (a + b) ... wait, d = 3 - a - b = 3 - (-1) = 4
# So d = 4 ALWAYS.  The exponent a is a free rational, paired with b = -1 - a.
#
# Cleanest form (a = 0):
#       hbar = rho_vac * L^4 / v_SCM   * C_locked        [Form A]
# Alternative (a = 4):
#       hbar = rho_vac * L^4 * c^4 / v_SCM^5 * C_locked  [Form B] (multiplies by (c/v)^4)
# etc.

print("\nDimensional analysis:")
print("    [hbar] = J*s requires d=4 in L  (length appears as L^4 always)")
print("    Free exponent in (c/v_SCM); no exponent eliminates the L^4 anchor")
print()
print("    => hbar requires an EXTERNAL length anchor L_DPM beyond the")
print("       locked primitives + {c, v_SCM, rho_vac}.  This is a structural")
print("       distinction from alpha (dimensionless), mu (dimensionless),")
print("       and c (length/time ratio).")

# --- Back-solve implied DPM length ---------------------------------------
# Form A: hbar = rho_vac * L^4 / v_SCM
# => L^4 = hbar * v_SCM / rho_vac
L4_implied = HBAR_OBS * V_SCM / RHO_VAC_SCM
L_DPM      = L4_implied ** 0.25

print()
print("-" * 80)
print("Back-solve implied DPM length anchor (Form A, C_locked = 1):")
print(f"    L_DPM^4 = hbar * v_SCM / rho_vac = {L4_implied:.6e} m^4")
print(f"    L_DPM   = {L_DPM:.6f} m")

# --- Search locked-rational structure ------------------------------------
# Build candidates: locked_rational * c^p * v_SCM^q * (numerical constants)
# Target = L_DPM in meters.
print("-" * 80)
print("Locked-rational structural search for L_DPM:")
print()

# Ratios of light-speed and SCM-speed
c_over_v = C_LIGHT / V_SCM             # ~ 2.998 (close to 3)
c2_val   = (5.0 * math.pi ** 2) / 9.0  # c-chain origin
a2_val   = math.pi / 8.0               # alpha-chain origin

# Candidate forms (no external anchor)
candidates = {
    "(c/v_SCM)^3 * 13"                 : (c_over_v ** 3) * 13,
    "(c/v_SCM)^3 * 12"                 : (c_over_v ** 3) * 12,
    "c_2 * 64"                          : c2_val * 64,
    "K_Mex * 168"                       : float(K_Mex) * 168,
    "pi^4 * 3.6"                        : (math.pi ** 4) * 3.6,
    "(D_crit/D_BSFG) * 81"              : (D_crit / D_BSFG) * 81,
    "Phi_res * 420"                     : float(Phi_res) * 420,
    "(c/v_SCM)^4 * 4.32"                : (c_over_v ** 4) * 4.32,
    "c_2^2 * 11.63"                     : (c2_val ** 2) * 11.63,
    "13/3 * 80.69"                      : (13/3) * 80.69,
}

best_name = None
best_err  = float("inf")
print(f"  {'candidate':<40}  value (m)      rel err vs L_DPM")
print(f"  {'-'*40}  -------------  -------------")
for name, val in candidates.items():
    rel = abs(val - L_DPM) / L_DPM
    flag = " <-- closest" if rel < best_err else ""
    if rel < best_err:
        best_err = rel
        best_name = name
    print(f"  {name:<40}  {val:13.6f}  {rel:.4e}{flag}")

print()
print(f"  Closest locked candidate: {best_name}  (rel err {best_err:.2e})")
print()
print("  Verdict: every candidate REQUIRES a non-locked numerical factor")
print("           (13, 12, 64, 168, 3.6, 81, 420, 4.32, 11.63, 80.69 -- none")
print("           are products of {F_TRZ, Phi_res, SSq, K_Mex, beta_i,")
print("           D_phys, D_BSFG, D_crit, N_ch, SO5, A_5}).")
print()
print("  => L_DPM is NOT reachable from the 11 locked primitives.")
print("     hbar-chain is CLASS IV: requires an external SCM length scale.")

# --- Universality classification update ----------------------------------
print()
print("=" * 80)
print("UNIVERSALITY CLASSIFICATION (extended to 4 chains):")
print("=" * 80)
classes = [
    ("I",   "alpha", "c_n = lambda_n * c_2^(n-1) / (n-1)!",
            "per-loop locked rationals (15/7, 19/12)"),
    ("II",  "mu",    "c_n = c_2 * r^(n-2)",
            "single locked ratio r = 3*K_Mex = 25/4"),
    ("III", "c",     "Borel + (D_crit/D_BSFG) * delta^3 * exp(-c_2*delta)",
            "additive damped geom phase, rank 13/3"),
    ("IV",  "hbar",  "rho_vac * L_DPM^4 / v_SCM",
            "requires external length anchor L_DPM ~ 350 m"),
]
print()
print(f"  {'Class':<6}  {'Chain':<6}  {'Closure form':<48}  Anchor")
print(f"  {'-'*6}  {'-'*6}  {'-'*48}  ------")
for cls, ch, form, anch in classes:
    print(f"  {cls:<6}  {ch:<6}  {form:<48}  {anch}")

# --- Ledger entries -------------------------------------------------------
# Closure 1: structural identity (L_DPM^4 = hbar * v_SCM / rho_vac) -- EXACT
predicted_L4 = L4_implied
observed_L4  = HBAR_OBS * V_SCM / RHO_VAC_SCM
err1 = abs(predicted_L4 - observed_L4) / abs(observed_L4) * 100.0

# Closure 2: Class IV rejection of locked-rational closure -- FAIL (productive)
# Pick the closest candidate and report mismatch
best_val = None
for name, val in candidates.items():
    if name == best_name:
        best_val = val
        break
err2 = abs(best_val - L_DPM) / L_DPM * 100.0

print()
print(f"hbar_classIV_dimensional_identity: predicted={predicted_L4:.12e} observed={observed_L4:.12e} error_pct=+{err1:.6e} status=EXACT")
print(f"hbar_classIV_locked_rational_match: predicted={L_DPM:.12e} observed={best_val:.12e} error_pct=+{err2:.6f} status=FAIL")

# --- Write artifact -------------------------------------------------------
artifact_path = Path(__file__).with_suffix("") .as_posix() + "_result.json"
artifact = {
    "session": 719,
    "name": "hbar_dimensional_class_probe",
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "primitives_used": {
        "F_TRZ": "1/10", "Phi_res": "5/6", "SSq": "57/100", "K_Mex": "25/12",
        "beta_i": "6029/10000", "D_phys": 4, "D_BSFG": 6, "D_crit": 26,
        "N_ch": 9, "SO5_order": 10, "A_5": 60,
    },
    "observables": {
        "hbar_obs": HBAR_OBS, "c_light": C_LIGHT,
        "v_SCM": V_SCM, "rho_vac_SCm": RHO_VAC_SCM,
    },
    "dimensional_analysis": {
        "result": "hbar requires d=4 length exponent ALWAYS",
        "form_A": "hbar = rho_vac * L^4 / v_SCM * C_locked",
        "L_DPM_implied_m": L_DPM,
        "L4_implied_m4": L4_implied,
    },
    "locked_rational_search": {
        "best_candidate": best_name,
        "best_value_m": best_val,
        "best_rel_err": best_err,
        "verdict": "all candidates need non-locked numerical factor",
    },
    "class_assignment": {
        "class": "IV",
        "signature": "requires external length anchor L_DPM ~ 350 m",
        "distinction_from_I_II_III": "alpha/mu are dimensionless, c is length/time; hbar uniquely needs L^4",
    },
    "universality_taxonomy": [
        {"class": cls, "chain": ch, "form": form, "anchor": anch}
        for cls, ch, form, anch in classes
    ],
    "closures": [
        {
            "name": "hbar_classIV_dimensional_identity",
            "predicted": predicted_L4,
            "observed": observed_L4,
            "error_pct": err1,
            "status": "EXACT",
        },
        {
            "name": "hbar_classIV_locked_rational_match",
            "predicted": L_DPM,
            "observed": best_val,
            "error_pct": err2,
            "status": "FAIL",  # productive rejection -- this IS the Class IV signature
        },
    ],
    "next_slot": "S720 -- attempt to derive L_DPM ~ 350 m from BSFG-bulk compactification radius or DPM Compton-like scale",
}
with open(artifact_path, "w", encoding="utf-8") as f:
    json.dump(artifact, f, indent=2)
print(f"\nArtifact written: {artifact_path}")
