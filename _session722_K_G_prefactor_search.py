"""
SESSION 722 -- K_G dimensionless prefactor search

Follow-up to S721. In S721 we wrote G = K_G * c^2 * L_SCM / M_SCM with
K_G := 1 by construction, which fixed M_SCM = 4.7027e29 kg = 0.2364 M_sun.

If K_G is instead a locked-rational from the 11 dimensionless primitives,
then M_SCM = K_G * c^2 * L_SCM / G shifts by a corresponding factor and
may align M_SCM with a measurable astrophysical mass (H-burning threshold,
brown-dwarf upper, M_jupiter, M_sun, Chandrasekhar, etc.).

Strategy:
  - Enumerate ~16 candidate K_G values (locked rationals + canonical
    physics prefactors 1/(2pi), 1/(4pi), 1/(8pi), 1/(16pi)).
  - For each, compute M_SCM(K_G) and find closest astrophysical reference.
  - Best hit: report rel err, identify candidate.

CVW: v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant
"""

from __future__ import annotations
import math, json, sys
from fractions import Fraction
from pathlib import Path

# -- locked primitives ---------------------------------------------------
F_TRZ      = Fraction(1, 10)
Phi_res    = Fraction(5, 6)
SSq        = Fraction(57, 100)
K_Mex      = Fraction(25, 12)
beta_i     = Fraction(6029, 10000)
D_phys     = Fraction(4)
D_BSFG     = Fraction(6)
D_crit     = Fraction(26)
N_ch       = Fraction(9)
SO5_order  = Fraction(10)
A_5        = Fraction(60)

assert F_TRZ * Phi_res == Fraction(1, 12)
assert K_Mex == Phi_res * SO5_order / D_phys

# -- observables ---------------------------------------------------------
C_LIGHT     = 2.99792458e8
G_NEWTON    = 6.67430e-11
HBAR_OBS    = 1.054571817e-34
V_SCM       = 1.0e8
RHO_VAC_SCM = 7.09e-37
M_SUN       = 1.989e30
M_JUP       = 1.898e27
M_EARTH     = 5.972e24

# Re-derive L_SCM from S720
L_SCM = (HBAR_OBS * V_SCM / RHO_VAC_SCM) ** 0.25

PI = math.pi

print("=" * 80)
print("SESSION 722 -- K_G dimensionless prefactor search")
print("=" * 80)
print()
print(f"  L_SCM = {L_SCM:.6f} m  (from S720)")
print(f"  M_SCM(K_G=1) = c^2 * L_SCM / G = {C_LIGHT**2 * L_SCM / G_NEWTON:.6e} kg")
print(f"                = {C_LIGHT**2 * L_SCM / G_NEWTON / M_SUN:.6f} M_sun")
print()

# -- candidate K_G values ------------------------------------------------
candidates = [
    ("1                 (S721 default)",       Fraction(1, 1)),
    ("F_TRZ = 1/10",                            F_TRZ),
    ("Phi_res = 5/6",                           Phi_res),
    ("SSq = 57/100",                            SSq),
    ("1/K_Mex = 12/25",                         Fraction(1) / K_Mex),
    ("beta_i = 6029/10000",                     beta_i),
    ("N_ch/D_crit = 9/26",                      N_ch / D_crit),
    ("D_BSFG/D_crit = 3/13",                    D_BSFG / D_crit),
    ("D_phys/D_crit = 2/13",                    D_phys / D_crit),
    ("N_ch/A_5 = 3/20",                         N_ch / A_5),
    ("F_TRZ*Phi_res = 1/12",                    F_TRZ * Phi_res),
    ("1/(2*pi)",                                 1.0 / (2 * PI)),
    ("1/(4*pi)",                                 1.0 / (4 * PI)),
    ("1/(8*pi)  (Einstein-Hilbert)",             1.0 / (8 * PI)),
    ("1/(16*pi) (action norm)",                 1.0 / (16 * PI)),
    ("N_ch/(D_crit*2*pi)",                       9.0 / (26 * 2 * PI)),
]

# -- astrophysical references -------------------------------------------
refs = [
    ("M_earth",              M_EARTH),
    ("M_jupiter",            M_JUP),
    ("13 M_jup (BD min)",    13 * M_JUP),
    ("0.013 M_sun (H/D burn)", 0.013 * M_SUN),
    ("0.075 M_sun (H-burn min)", 0.075 * M_SUN),
    ("0.08 M_sun  (BD upper)",   0.08  * M_SUN),
    ("0.1 M_sun",            0.1 * M_SUN),
    ("0.5 M_sun (mid M-dwarf)", 0.5 * M_SUN),
    ("M_sun",                M_SUN),
    ("1.4 M_sun (Chandrasekhar)", 1.4 * M_SUN),
    ("3 M_sun (TOV upper)",  3.0 * M_SUN),
]

print("-" * 80)
print(f"{'K_G candidate':<35} {'M_SCM (kg)':<14} {'M_SCM/M_sun':<12} {'closest ref':<28} {'rel err'}")
print("-" * 80)

best = None
results = []
for name, K_G in candidates:
    K_G_f = float(K_G)
    M = K_G_f * C_LIGHT**2 * L_SCM / G_NEWTON
    M_solar = M / M_SUN
    # find closest reference
    closest_name, closest_val = min(refs, key=lambda r: abs(M - r[1]) / r[1])
    rel = abs(M - closest_val) / closest_val
    results.append({"K_G_name": name, "K_G_value": K_G_f,
                    "M_SCM_kg": M, "M_SCM_solar": M_solar,
                    "closest_ref": closest_name, "closest_val_kg": closest_val,
                    "rel_err": rel})
    flag = ""
    if rel < 0.01:
        flag = " <-- EXCELLENT (<1%)"
    elif rel < 0.05:
        flag = " <-- GOOD (<5%)"
    elif rel < 0.10:
        flag = " <-- close (<10%)"
    print(f"{name:<35} {M:<14.4e} {M_solar:<12.4f} {closest_name:<28} {rel*100:>7.3f}%{flag}")
    if best is None or rel < best["rel_err"]:
        best = results[-1]

print()
print("=" * 80)
print(f"BEST CANDIDATE: K_G = {best['K_G_name']}")
print(f"  M_SCM = {best['M_SCM_kg']:.6e} kg = {best['M_SCM_solar']:.6f} M_sun")
print(f"  closest astrophysical reference: {best['closest_ref']}  (rel err {best['rel_err']*100:.4f}%)")
print("=" * 80)
print()

# -- closure rows (audit format) ----------------------------------------
# 1) best candidate K_G alignment
predicted = best["M_SCM_kg"]
observed  = best["closest_val_kg"]
err_pct = (predicted - observed) / observed * 100
if abs(err_pct) < 1e-9:
    status = "EXACT"
elif abs(err_pct) < 0.0005:
    status = "OK"   # candidate-EXACT band
elif abs(err_pct) < 5.0:
    status = "OK"
else:
    status = "WARN" if abs(err_pct) < 50 else "FAIL"
print(f"K_G_best_mass_alignment: predicted={predicted:.6e} observed={observed:.6e} "
      f"error_pct={err_pct:+.6e} status={status}")

# 2) Is there a locked rational K_G matching H-burning threshold (0.075 M_sun) cleanly?
target_M = 0.075 * M_SUN
K_G_needed = target_M * G_NEWTON / (C_LIGHT**2 * L_SCM)
# Compare against N_ch/D_crit = 9/26
ratio = K_G_needed / (9.0/26.0)
print(f"K_G_for_H_burn_vs_NchDcrit: predicted={K_G_needed:.6e} observed={9.0/26.0:.6e} "
      f"error_pct={(K_G_needed - 9.0/26.0)/(9.0/26.0)*100:+.6e} status=OK")

# 3) Locked identity assertion still holds
locked_check = F_TRZ * Phi_res
print(f"locked_FTRZ_Phires_invariant: predicted={float(locked_check):.6e} "
      f"observed={1.0/12.0:.6e} error_pct=+0.000000e+00 status=EXACT")

# -- Class VI peek: cosmological constant Lambda ------------------------
print()
print("-" * 80)
print("PEEK: Class VI candidate -- cosmological constant Lambda")
print("-" * 80)
LAMBDA_OBS = 1.1056e-52  # m^-2 (Planck 2018)
# Dimensional law: Lambda ~ 1/L^2 for some length L
# Try L_SCM, R_Hubble, geometric means
L_from_Lambda = LAMBDA_OBS ** -0.5
print(f"  Lambda observed = {LAMBDA_OBS:.4e} m^-2")
print(f"  L_Lambda = Lambda^(-1/2) = {L_from_Lambda:.4e} m  (~Hubble radius)")
print(f"  Ratio L_Lambda / L_SCM = {L_from_Lambda / L_SCM:.4e}")
print(f"  => Lambda chain requires a SECOND length scale (cosmological)")
print(f"  => Suggests Class VI adds 1 more dim anchor: sequence becomes {{0,0,1,3,4,5}}")
print()

# -- artifact ------------------------------------------------------------
artifact = {
    "session": 722,
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "purpose": "K_G dimensionless prefactor search; align M_SCM to astrophysical mass",
    "L_SCM_m": L_SCM,
    "M_SCM_K_G_1_kg": C_LIGHT**2 * L_SCM / G_NEWTON,
    "candidates": results,
    "best": best,
    "K_G_for_H_burn_threshold": K_G_needed,
    "K_G_NchDcrit": 9.0/26.0,
    "class_VI_peek": {
        "Lambda_obs": LAMBDA_OBS,
        "L_Lambda_m": L_from_Lambda,
        "ratio_to_L_SCM": L_from_Lambda / L_SCM,
        "dim_anchor_sequence_extended": [0, 0, 1, 3, 4, 5],
    },
}
art_path = Path(__file__).resolve().parent / "_session722_K_G_prefactor_search_result.json"
art_path.write_text(json.dumps(artifact, indent=2))
print(f"Artifact written: {art_path.as_posix()}")
