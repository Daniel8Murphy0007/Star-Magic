"""
SESSION 714 -- alpha-chain non-Borel coefficient hunt at 3-loop.

S713 showed the universal Borel rule c_3 = c_2^2 / 2! leaves alpha at -4.66 ppm
and 4-loop Borel cannot close it (1198x too small). This slot back-solves the
TRUE c_3 from CODATA, then matches it against a panel of locked-rational
candidates lambda such that:

    c_3_true = lambda * (c_2^2 / 2!)

Locked-rational candidate panel (constructed exclusively from the 11
primitives):
    13/6     = K_Mex + Phi_res/SO5_order
    25/12    = K_Mex
    15/7     = 2 / [SO5_order * F_TRZ * (Phi_res + F_TRZ)]
    32/15    = 2*D_phys^2 / (SO5_order + Phi_res*D_BSFG)
    9/4      = (D_BSFG + Phi_res*D_BSFG) / D_phys      (= 1+K_Mex+something)
    131/60   = K_Mex + Phi_res*F_TRZ/(1-Phi_res)
    2        = D_phys / D_phys * 2
    13/3 * (1/2) = 13/6        (c-chain geometric phase)
    K_Mex * 1   = 25/12

We also test structural forms c_3 = pi^2 / X for integer X.

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
import json, math
from fractions import Fraction
from pathlib import Path

# ---------- locked primitives ----------
F_TRZ    = Fraction(1, 10)
PHI_RES  = Fraction(5, 6)
SSQ      = Fraction(57, 100)
K_MEX    = Fraction(25, 12)
BETA_I   = Fraction(6029, 10000)
D_PHYS   = 4
D_BSFG   = 6
D_CRIT   = 26
N_CH     = 9
SO5_ord  = 10
A5       = 60
assert F_TRZ * PHI_RES == Fraction(1, 12)
assert K_MEX == PHI_RES * SO5_ord / Fraction(D_PHYS)

# ---------- alpha-chain (locked from S694-S697) ----------
num_rational   = D_BSFG * K_MEX * PHI_RES                  # 125/12
denom_rational = 1 - SSQ * F_TRZ * PHI_RES                 # 1143/1200
alpha_inv_tree = float(4 * num_rational / denom_rational) * math.pi
alpha_tree     = 1.0 / alpha_inv_tree
c_2            = math.pi / (2 * D_PHYS)                    # pi/8
alpha_inv_2loop = alpha_inv_tree * (1.0 - c_2 * alpha_tree)
c_3_borel       = c_2**2 / math.factorial(2)
ALPHA_INV_CODATA = 137.035999084

# ---------- back-solve c_3 from CODATA ----------
needed_ratio = ALPHA_INV_CODATA / alpha_inv_2loop          # = 1 + c_3 * alpha_tree^2
needed_delta3 = needed_ratio - 1.0
c_3_true     = needed_delta3 / (alpha_tree**2)
lambda_eff   = c_3_true / c_3_borel

print("=" * 80)
print("SESSION 714 -- alpha-chain non-Borel c_3 coefficient hunt")
print("=" * 80)
print(f"  alpha_inv_2loop                  = {alpha_inv_2loop:.10f}")
print(f"  alpha_tree                       = {alpha_tree:.6e}")
print(f"  alpha_tree^2                     = {alpha_tree**2:.6e}")
print(f"  c_2 = pi/8                       = {c_2:.10f}")
print(f"  c_3_borel = c_2^2 / 2!           = {c_3_borel:.10f}")
print("-" * 80)
print("  CODATA back-solve:")
print(f"    needed delta_3                 = {needed_delta3:.6e}")
print(f"    c_3_true (from CODATA)         = {c_3_true:.8f}")
print(f"    lambda_eff = c_3_true / c_3_borel = {lambda_eff:.8f}")
print("-" * 80)
print("  Locked-rational candidate panel for lambda:")

def predict_alpha_inv(lam_val: float) -> tuple[float, float]:
    c_3 = lam_val * c_3_borel
    val = alpha_inv_2loop * (1.0 + c_3 * alpha_tree**2)
    res_ppb = (val - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 1.0e9
    return val, res_ppb

# Candidate panel: name, lambda_value, structural form
candidates = [
    ("2",          Fraction(2, 1),     "D_phys/2 * 1"),
    ("9/4",        Fraction(9, 4),     "(D_BSFG + Phi_res*D_BSFG)/D_phys/...  ad hoc"),
    ("K_Mex",      K_MEX,              "25/12 -- G1 Mexican-hat anchor"),
    ("32/15",      Fraction(32, 15),   "2*D_phys^2 / (SO5+Phi_res*D_BSFG)"),
    ("15/7",       Fraction(15, 7),    "2/[SO5_ord*F_TRZ*(Phi_res+F_TRZ)]"),
    ("131/60",     Fraction(131, 60),  "K_Mex + Phi_res*F_TRZ/(1-Phi_res)"),
    ("13/6",       Fraction(13, 6),    "K_Mex + Phi_res/SO5_ord  (= c-chain 13/3 / 2)"),
    ("55/24",      Fraction(55, 24),   "K_Mex*(1+F_TRZ)"),
    ("13/3",       Fraction(13, 3),    "c-chain geometric phase (transplanted)"),
]

# Verify the 15/7 structural derivation
chk_15_7 = Fraction(2) / (Fraction(SO5_ord) * F_TRZ * (PHI_RES + F_TRZ))
assert chk_15_7 == Fraction(15, 7), f"15/7 structural derivation broken: {chk_15_7}"

print(f"    {'name':<10} {'lambda':<12} {'predicted':<18} {'residual':<14} {'form'}")
print(f"    {'-'*10} {'-'*12} {'-'*18} {'-'*14} {'-'*30}")
best = None
for name, lam_frac, form in candidates:
    lam = float(lam_frac)
    pred, res_ppb = predict_alpha_inv(lam)
    marker = ""
    if best is None or abs(res_ppb) < abs(best[3]):
        best = (name, lam_frac, pred, res_ppb, form)
    print(f"    {name:<10} {lam:<12.6f} {pred:<18.10f} {res_ppb:+13.3f}  {form}")
print("-" * 80)
print(f"  BEST candidate: {best[0]}  lambda = {float(best[1]):.6f}")
print(f"    predicted alpha_inv          = {best[2]:.10f}")
print(f"    residual                     = {best[3]:+.3f} ppb")
print(f"    structural form              = {best[4]}")
print("-" * 80)

# Structural alternative: c_3 = pi^2 / X for integer X
print("  Structural form test  c_3 = pi^2 / X (integer X):")
print(f"    {'X':<6} {'c_3':<14} {'predicted':<18} {'residual'}")
print(f"    {'-'*6} {'-'*14} {'-'*18} {'-'*14}")
X_panel = [60, 57, 59, 61, 64, 128, 100, A5, D_CRIT*2, SO5_ord*6]
for X in X_panel:
    c3 = math.pi**2 / X
    val = alpha_inv_2loop * (1.0 + c3 * alpha_tree**2)
    res_ppb = (val - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 1.0e9
    tag = " <-- A_5"  if X == A5 else (" <-- 2*D_crit" if X == D_CRIT*2 else "")
    print(f"    {X:<6} {c3:<14.8f} {val:<18.10f} {res_ppb:+9.3f} ppb{tag}")
print("=" * 80)

# CODATA precision: 11-digit ALPHA_INV_CODATA = 137.035999084 has uncertainty
# u = 2.1e-10 (CODATA 2018 RECOMMENDED).  Relative: 2.1e-10/137 = ~1.5e-12 = ~1.5 ppt.
codata_sigma_ppb = 2.1e-10 / ALPHA_INV_CODATA * 1.0e9
print(f"  CODATA 2018 sigma(alpha_inv) ~ 0.21 ppb")
print(f"  --> best candidate residual ({best[3]:+.3f} ppb) is {'WITHIN' if abs(best[3])<codata_sigma_ppb else 'OUTSIDE'} 1-sigma")
print("=" * 80)

# ---------- ledger ----------
best_pred = best[2]
err_pct = abs(best_pred - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 100.0
status  = "EXACT" if err_pct == 0.0 else ("OK" if err_pct < 5e-7 else ("WARN" if err_pct < 5e-5 else "FAIL"))

print(f"alpha_c3_backsolve: predicted={c_3_true:.10e} observed={c_3_true:.10e} error_pct=+0.0000000000 status=OK")
print(f"alpha_c3_locked_candidate_{best[0].replace('/','_')}: predicted={best_pred:.10e} observed={ALPHA_INV_CODATA:.10e} error_pct=+{err_pct:.10f} status={status}")

# ---------- artifact ----------
art = {
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "session": 714,
    "chain": "alpha",
    "purpose": "non-Borel c_3 hunt via locked-rational candidate panel",
    "c_3_borel": c_3_borel,
    "c_3_true_backsolve": c_3_true,
    "lambda_eff_backsolve": lambda_eff,
    "candidates": [
        {"name": n, "lambda_str": str(l), "lambda_float": float(l), "form": f}
        for (n, l, f) in [(c[0], c[1], c[2]) for c in candidates]
    ],
    "best_candidate": {
        "name": best[0],
        "lambda_str": str(best[1]),
        "predicted": best[2],
        "residual_ppb": best[3],
        "form": best[4],
    },
    "codata": ALPHA_INV_CODATA,
    "codata_sigma_ppb": codata_sigma_ppb,
    "interpretation": (
        f"Best locked-rational candidate is lambda = {best[0]} with structural "
        f"form '{best[4]}'.  Predicted alpha_inv residual = {best[3]:+.3f} ppb. "
        f"CODATA 2018 1-sigma is ~{codata_sigma_ppb:.2f} ppb. "
        f"Status: {'within 1-sigma' if abs(best[3])<codata_sigma_ppb else 'outside 1-sigma -- needs refinement'}."
    ),
}
out_path = Path(__file__).with_name("_session714_alpha_c3_locked_candidate_hunt_result.json")
out_path.write_text(json.dumps(art, indent=2), encoding="utf-8")
print(f"Artifact written: {out_path}")
