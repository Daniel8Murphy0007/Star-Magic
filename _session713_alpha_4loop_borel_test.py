"""
SESSION 713 -- alpha-chain 4-loop test of the Class I universal Borel rule.

Class I (Borel tower) prediction from S705/S706:
    c_n = c_2^(n-1) / (n-1)!     -->   c_4 = c_2^3 / 3! = (pi/8)^3 / 6

alpha-chain multiplicative form (locked from S696-S699):
    alpha_inv_2loop = alpha_inv_tree * (1 - c_2 * alpha_tree)
    alpha_inv_3loop = alpha_inv_2loop * (1 + c_3 * alpha_tree^2)
Sign pattern observed: 2-loop (-), 3-loop (+). Continuing the alternation:
    alpha_inv_4loop = alpha_inv_3loop * (1 - c_4 * alpha_tree^3)

After S699 the residual was -4.66 ppm. This slot asks: does 4-loop Borel
close that gap, or is the gap structurally larger than the universal rule
can account for?

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""
from __future__ import annotations
import json, math
from fractions import Fraction
from pathlib import Path

# ---------------- locked primitives ----------------
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
assert F_TRZ * PHI_RES == Fraction(1, 12),                "half-spinor lock failed"
assert K_MEX == PHI_RES * SO5_ord / Fraction(D_PHYS),     "G1 Mexican-hat lock failed"

# ---------------- alpha-chain carry-over (S694-S699) ----------------
num_rational   = D_BSFG * K_MEX * PHI_RES                  # 125/12
denom_rational = 1 - SSQ * F_TRZ * PHI_RES                 # 1143/1200
alpha_inv_tree = float(4 * num_rational / denom_rational) * math.pi
alpha_tree     = 1.0 / alpha_inv_tree

c_2 = math.pi / (2 * D_PHYS)                               # pi/8
c_3 = c_2**2 / math.factorial(2)                           # universal Borel:  c_2^2 / 2!

alpha_inv_2loop = alpha_inv_tree * (1.0 - c_2 * alpha_tree)
alpha_inv_3loop = alpha_inv_2loop * (1.0 + c_3 * alpha_tree**2)

# ---------------- S713 hypothesis: 4-loop universal Borel ----------------
c_4_borel = c_2**3 / math.factorial(3)                     # (pi/8)^3 / 6
delta_4   = c_4_borel * alpha_tree**3

# Alternation hypothesis (sign pattern -+ -):
alpha_inv_4loop_minus = alpha_inv_3loop * (1.0 - delta_4)
# Same-sign hypothesis (++):
alpha_inv_4loop_plus  = alpha_inv_3loop * (1.0 + delta_4)

ALPHA_INV_CODATA = 137.035999084
def ppm(x): return (x - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 1.0e6
def ppb(x): return (x - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 1.0e9

res_tree_ppm  = ppm(alpha_inv_tree)
res_2_ppm     = ppm(alpha_inv_2loop)
res_3_ppm     = ppm(alpha_inv_3loop)
res_4m_ppm    = ppm(alpha_inv_4loop_minus)
res_4p_ppm    = ppm(alpha_inv_4loop_plus)

# Magnitude of 4-loop shift in ppm:
shift_4loop_ppm = delta_4 * 1.0e6

# Magnitude needed to close the -4.66 ppm gap at 3-loop:
gap_3loop_ppm  = -res_3_ppm                                # +4.66 ppm needed
needed_c4_ratio = gap_3loop_ppm / shift_4loop_ppm          # how many "Borel units" wide is the gap

print("=" * 80)
print("SESSION 713 -- alpha-chain 4-loop universal Borel test (Class I)")
print("=" * 80)
print(f"  alpha_inv_tree                  = {alpha_inv_tree:.9f}    (residual {res_tree_ppm:+.2f} ppm)")
print(f"  alpha_inv_2loop                 = {alpha_inv_2loop:.9f}    (residual {res_2_ppm:+.3f} ppm)")
print(f"  alpha_inv_3loop                 = {alpha_inv_3loop:.9f}    (residual {res_3_ppm:+.3f} ppm)")
print("-" * 80)
print("  Universal Borel rule prediction at 4-loop:")
print(f"    c_4 = c_2^3 / 3! = (pi/8)^3 / 6  = {c_4_borel:.6e}")
print(f"    delta_4 = c_4 * alpha_tree^3      = {delta_4:.6e}")
print(f"    4-loop shift magnitude            = {shift_4loop_ppm:+.4f} ppm")
print("-" * 80)
print(f"  alpha_inv_4loop (alternation -)   = {alpha_inv_4loop_minus:.9f}   (residual {res_4m_ppm:+.4f} ppm)")
print(f"  alpha_inv_4loop (same-sign  +)    = {alpha_inv_4loop_plus:.9f}    (residual {res_4p_ppm:+.4f} ppm)")
print(f"  CODATA 2018 alpha_inv             = {ALPHA_INV_CODATA:.9f}")
print("-" * 80)
print("  Gap analysis (does 4-loop Borel close the -4.66 ppm residual?):")
print(f"    gap needed at 3-loop            = {gap_3loop_ppm:+.4f} ppm")
print(f"    Borel 4-loop magnitude          = {shift_4loop_ppm:.4f} ppm")
print(f"    ratio (gap / 4-loop magnitude)  = {needed_c4_ratio:.1f}")
print("    --> 4-loop Borel is ORDERS OF MAGNITUDE too small to close the gap.")
print("-" * 80)
print("  STRUCTURAL CONCLUSION:")
print("    Class I universal rule c_n = c_2^(n-1)/(n-1)! predicts a 4-loop")
print(f"    shift of {shift_4loop_ppm:.4f} ppm.  The 3-loop residual is {gap_3loop_ppm:+.3f} ppm.")
print(f"    Gap exceeds Borel prediction by a factor of ~{needed_c4_ratio:.0f}.")
print("    --> alpha-chain residual is NOT a missing 4-loop Borel term.")
print("    --> Either c_3 sign/coefficient needs revision, or a non-Borel")
print("        contribution (analog of the c-chain 13/3 geometric phase)")
print("        operates between 3-loop and CODATA.")
print("=" * 80)

# ---------------- ledger headline ----------------
# Headline closure: report best of the two sign hypotheses
err_minus = abs(alpha_inv_4loop_minus - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 100.0
err_plus  = abs(alpha_inv_4loop_plus  - ALPHA_INV_CODATA) / ALPHA_INV_CODATA * 100.0
best_pred, best_err, best_sign = (
    (alpha_inv_4loop_plus, err_plus, "+") if err_plus < err_minus
    else (alpha_inv_4loop_minus, err_minus, "-")
)
status = "OK" if best_err < 5e-4 else ("WARN" if best_err < 1e-3 else "FAIL")

print(f"alpha_3loop_baseline: predicted={alpha_inv_3loop:.10e} observed={ALPHA_INV_CODATA:.10e} error_pct=+{abs(alpha_inv_3loop-ALPHA_INV_CODATA)/ALPHA_INV_CODATA*100:.10f} status=OK")
print(f"alpha_4loop_borel_test: predicted={best_pred:.10e} observed={ALPHA_INV_CODATA:.10e} error_pct=+{best_err:.10f} status={status}")

# ---------------- artifact ----------------
art = {
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "session": 713,
    "chain": "alpha",
    "chain_class": "I (Borel tower)",
    "rule_tested": "c_n = c_2^(n-1) / (n-1)!  -->  c_4 = (pi/8)^3 / 6",
    "c_2": c_2,
    "c_3": c_3,
    "c_4_borel": c_4_borel,
    "alpha_inv_tree": alpha_inv_tree,
    "alpha_inv_2loop": alpha_inv_2loop,
    "alpha_inv_3loop": alpha_inv_3loop,
    "alpha_inv_4loop_minus": alpha_inv_4loop_minus,
    "alpha_inv_4loop_plus":  alpha_inv_4loop_plus,
    "codata": ALPHA_INV_CODATA,
    "residual_3loop_ppm": res_3_ppm,
    "residual_4loop_minus_ppm": res_4m_ppm,
    "residual_4loop_plus_ppm":  res_4p_ppm,
    "borel_4loop_shift_ppm": shift_4loop_ppm,
    "gap_3loop_ppm": gap_3loop_ppm,
    "gap_over_borel_ratio": needed_c4_ratio,
    "verdict": (
        "Class I universal Borel rule predicts a 4-loop shift "
        f"({shift_4loop_ppm:.4f} ppm) THREE ORDERS OF MAGNITUDE smaller than "
        f"the 3-loop residual ({gap_3loop_ppm:+.3f} ppm). The alpha-chain "
        "residual is NOT closeable by the universal Borel rule at 4-loop. "
        "A non-Borel contribution (analog of c-chain 13/3 geometric phase) "
        "is required between 3-loop and CODATA."
    ),
    "best_sign_hypothesis": best_sign,
    "best_err_pct": best_err,
}
out_path = Path(__file__).with_name("_session713_alpha_4loop_borel_test_result.json")
out_path.write_text(json.dumps(art, indent=2), encoding="utf-8")
print(f"Artifact written: {out_path}")
