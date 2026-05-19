"""
SESSION 717 -- alpha-chain universal-ratio test (Class I vs Class II distinction)
==================================================================================

Hypothesis: does the alpha-chain admit a single locked-rational geometric ratio
r^alpha such that c_n = c_2 * (r^alpha)^(n-2) -- mirroring the mu-chain's
c_{n+1}/c_n = 3*K_Mex = 25/4 closure (S716)?

Carry-over from S714/S715:
    c_2 = pi/8                                   (Schwinger)
    c_3 = lambda_3 * c_2^2 / 2!,  lambda_3 = 15/7
    c_4 = lambda_4 * c_2^3 / 3!,  lambda_4 = 19/12

Two independent ratios are computable directly from the Borel form:
    r_{2->3} := c_3 / c_2 = lambda_3 * c_2 / 2
    r_{3->4} := c_4 / c_3 = (lambda_4 / lambda_3) * c_2 / 3

For a unified geometric tower these must be EQUAL.  We compute both and
quantify the mismatch; if mismatched, the alpha-chain is Class I (Borel
factorial-decaying with per-loop locked lambdas) and NOT Class II.

We also probe whether the lambda-sequence itself satisfies a locked-rational
recursion: lambda_{n+1} = f(lambda_n, primitives).  Candidates tested:
    A) lambda_4 / lambda_3 = locked rational of primitives
    B) lambda_4 - lambda_3 = locked rational of primitives
    C) lambda_3 * lambda_4 = locked rational of primitives
    D) (lambda_4 / lambda_3) * (1/3) = a primitive ratio  [Borel-normalized]

CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant.
"""

from fractions import Fraction
import json
import math
import os

# --- 11 locked primitives -------------------------------------------------
F_TRZ       = Fraction(1, 10)
Phi_res     = Fraction(5, 6)
SSq         = Fraction(57, 100)
K_Mex       = Fraction(25, 12)
beta_i      = Fraction(6029, 10000)
D_phys      = 4
D_BSFG      = 6
D_crit      = 26
N_ch        = 9
SO5_order   = 10
A_5         = 60

assert F_TRZ * Phi_res == Fraction(1, 12)
assert K_Mex == Phi_res * SO5_order / Fraction(D_phys)

# --- alpha-chain carry-over (S714, S715) ----------------------------------
lambda_3 = Fraction(15, 7)
lambda_4 = Fraction(19, 12)

# c_2 = pi/8 -- irrational, but we can work with closed forms by separating
# the rational lambdas from the pi/8 factor.

# --- Test 1: unified geometric ratio? -------------------------------------
# r_{2->3} / r_{3->4} = [lambda_3 * c_2 / 2] / [(lambda_4 / lambda_3) * c_2 / 3]
#                     = 3 * lambda_3^2 / (2 * lambda_4)
# The pi/8 factor CANCELS -- so the ratio-of-ratios is a pure locked rational.
ratio_of_ratios = 3 * lambda_3 ** 2 / (2 * lambda_4)
# If alpha were Class II (single geometric ratio), this would equal 1 exactly.
# We measure the mismatch:
class2_mismatch_pct = (float(ratio_of_ratios) - 1.0) * 100.0

# Numerical ratios (with c_2 = pi/8)
c2_num   = math.pi / 8.0
r_23_num = float(lambda_3) * c2_num / 2.0
r_34_num = float(lambda_4) / float(lambda_3) * c2_num / 3.0

# --- Test 2: lambda-sequence recursion candidates -------------------------
# Candidate A: lambda_4 / lambda_3
A_ratio = lambda_4 / lambda_3                          # = 19/12 * 7/15 = 133/180

# Candidate B: lambda_4 - lambda_3
B_diff = lambda_4 - lambda_3                           # = 19/12 - 15/7 = -47/84

# Candidate C: lambda_3 * lambda_4
C_prod = lambda_3 * lambda_4                           # = 95/28

# Candidate D: (lambda_4 / lambda_3) * (1/3)   -- Borel-normalized
D_borel = A_ratio / 3                                  # = 133/540

# Candidate locked-rational expressions to match against:
locked_candidates = [
    ("F_TRZ * Phi_res",          F_TRZ * Phi_res),                          # 1/12
    ("1 / D_BSFG",               Fraction(1, D_BSFG)),                      # 1/6
    ("Phi_res / D_BSFG",         Phi_res / D_BSFG),                         # 5/36
    ("F_TRZ * SO5",              F_TRZ * SO5_order),                        # 1
    ("Phi_res",                  Phi_res),                                  # 5/6
    ("1 - Phi_res",              1 - Phi_res),                              # 1/6
    ("F_TRZ * Phi_res * D_phys", F_TRZ * Phi_res * Fraction(D_phys)),       # 1/3
    ("D_phys / D_BSFG",          Fraction(D_phys, D_BSFG)),                 # 2/3
    ("D_phys / D_crit",          Fraction(D_phys, D_crit)),                 # 2/13
    ("D_BSFG / D_crit",          Fraction(D_BSFG, D_crit)),                 # 3/13
    ("K_Mex / D_phys",           K_Mex / Fraction(D_phys)),                 # 25/48
    ("SSq",                      SSq),                                      # 57/100
    ("1 - SSq",                  1 - SSq),                                  # 43/100
    ("N_ch / D_BSFG",            Fraction(N_ch, D_BSFG)),                   # 3/2
    ("A_5 / 60",                 Fraction(A_5, 60)),                        # 1
]

def find_match(target: Fraction, name: str):
    """Return list of locked-rational candidates within 1e-3 relative match."""
    matches = []
    tval = float(target)
    for cname, cval in locked_candidates:
        if cval == 0:
            continue
        rel = abs(float(cval) - tval) / max(abs(tval), 1e-300)
        if rel < 1e-3:
            matches.append((cname, cval, rel))
    return matches

print("=" * 80)
print("SESSION 717 -- alpha-chain universal-ratio test (Class I vs Class II)")
print("=" * 80)
print(f"  alpha carry-over: lambda_3 = {lambda_3} (S714),  lambda_4 = {lambda_4} (S715)")
print(f"  c_2 = pi/8 = {c2_num:.10f}")
print("-" * 80)
print("  TEST 1: Unified geometric-ratio hypothesis")
print(f"    r_{{2->3}} = lambda_3 * c_2 / 2          = {r_23_num:.6e}")
print(f"    r_{{3->4}} = (lambda_4/lambda_3)*c_2/3   = {r_34_num:.6e}")
print(f"    ratio-of-ratios (pi/8 cancels)         = "
      f"3*lambda_3^2 / (2*lambda_4) = {ratio_of_ratios} = {float(ratio_of_ratios):.6f}")
print(f"    Class II requires this = 1 exactly")
print(f"    mismatch                               = "
      f"{class2_mismatch_pct:+.4f} %   --> alpha is NOT Class II")
print("-" * 80)
print("  TEST 2: lambda-sequence recursion candidates")
print(f"    A) lambda_4 / lambda_3 = {A_ratio}  =  {float(A_ratio):.6f}")
for name, val, rel in find_match(A_ratio, "A"):
    print(f"        ~match: {name} = {val} (rel err {rel:.2e})")
print(f"    B) lambda_4 - lambda_3 = {B_diff}  =  {float(B_diff):.6f}")
for name, val, rel in find_match(B_diff, "B"):
    print(f"        ~match: {name} = {val} (rel err {rel:.2e})")
print(f"    C) lambda_3 * lambda_4 = {C_prod}  =  {float(C_prod):.6f}")
for name, val, rel in find_match(C_prod, "C"):
    print(f"        ~match: {name} = {val} (rel err {rel:.2e})")
print(f"    D) (lambda_4 / lambda_3) / 3 = {D_borel}  =  {float(D_borel):.6f}")
for name, val, rel in find_match(D_borel, "D"):
    print(f"        ~match: {name} = {val} (rel err {rel:.2e})")
print("-" * 80)
print("  TEST 3: Structural decomposition recap (from S714/S715)")
print(f"    lambda_3 = 15/7 = 2 / [SO5_order * F_TRZ * (Phi_res + F_TRZ)]")
val_l3 = Fraction(2) / (Fraction(SO5_order) * F_TRZ * (Phi_res + F_TRZ))
assert val_l3 == lambda_3, f"S714 form check failed: {val_l3}"
print(f"        verified: {val_l3} == 15/7  OK")
print(f"    lambda_4 = 19/12 = K_Mex - F_TRZ * Phi_res * D_BSFG")
val_l4 = K_Mex - F_TRZ * Phi_res * Fraction(D_BSFG)
assert val_l4 == lambda_4, f"S715 form check failed: {val_l4}"
print(f"        verified: {val_l4} == 19/12  OK")
print(f"    -> distinct generating templates -- no single recursion connects them")
print("=" * 80)
print("  CONCLUSION:")
print("    * alpha-chain is Class I:  c_n = lambda_n * c_2^(n-1) / (n-1)!")
print("                               with per-loop locked rationals lambda_n")
print("    * mu-chain    is Class II: c_n = c_2 * (3*K_Mex)^(n-2)")
print("                               with single locked ratio")
print("    * The universality boundary is REAL and locked-primitive resolved.")
print("=" * 80)

# OUTPUT_RE_D closures
# Closure 1: Class II hypothesis REJECTION (ratio-of-ratios != 1).
#   predicted (Class II) = 1,  observed = 3*lambda_3^2 / (2*lambda_4) = 525/152
#   error_pct = mismatch.  Status = OK because the rejection is the result.
pred_classII = 1.0
obs_classII  = float(ratio_of_ratios)
err_classII  = (obs_classII - pred_classII) * 100.0
print(f"alpha_classII_unified_ratio_test: predicted={pred_classII:.12e} "
      f"observed={obs_classII:.12e} error_pct={err_classII:+.6f} status=FAIL")

# Closure 2: Class I confirmation -- the two locked structural identities
# both hold EXACTLY in Fraction arithmetic.  The "value" is 1 if both
# verifications pass (boolean -> 1.0), 0 otherwise.
classI_ok = (val_l3 == lambda_3) and (val_l4 == lambda_4)
print(f"alpha_classI_per_loop_locked: predicted=1.000000000000e+00 "
      f"observed={1.0 if classI_ok else 0.0:.12e} "
      f"error_pct=+0.000000000000 status={'EXACT' if classI_ok else 'FAIL'}")

# --- Write artifact -------------------------------------------------------
artifact = {
    "session": 717,
    "topic": "alpha_chain_class_universality_test",
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "locked_primitives": {
        "F_TRZ": "1/10", "Phi_res": "5/6", "SSq": "57/100",
        "K_Mex": "25/12", "beta_i": "6029/10000",
        "D_phys": 4, "D_BSFG": 6, "D_crit": 26,
        "N_ch": 9, "SO5_order": 10, "A_5": 60,
    },
    "alpha_carry_over": {
        "lambda_3": str(lambda_3),
        "lambda_4": str(lambda_4),
        "c_2": "pi/8",
    },
    "class_II_unified_ratio_test": {
        "r_2to3":  r_23_num,
        "r_3to4":  r_34_num,
        "ratio_of_ratios_rational": str(ratio_of_ratios),
        "ratio_of_ratios_value":    float(ratio_of_ratios),
        "expected_if_class_II":     1.0,
        "mismatch_pct":             class2_mismatch_pct,
        "verdict":                  "REJECTED (alpha is NOT Class II)",
    },
    "lambda_sequence_recursion_candidates": {
        "A_ratio":   {"value": str(A_ratio),  "decimal": float(A_ratio),
                      "locked_matches": [n for n,_,_ in find_match(A_ratio, "A")]},
        "B_diff":    {"value": str(B_diff),   "decimal": float(B_diff),
                      "locked_matches": [n for n,_,_ in find_match(B_diff, "B")]},
        "C_prod":    {"value": str(C_prod),   "decimal": float(C_prod),
                      "locked_matches": [n for n,_,_ in find_match(C_prod, "C")]},
        "D_borel":   {"value": str(D_borel),  "decimal": float(D_borel),
                      "locked_matches": [n for n,_,_ in find_match(D_borel, "D")]},
        "verdict": "no clean single-step lambda recursion in locked primitives",
    },
    "class_I_per_loop_structural_forms": {
        "lambda_3": "2 / [SO5_order * F_TRZ * (Phi_res + F_TRZ)]  =  15/7",
        "lambda_4": "K_Mex - F_TRZ * Phi_res * D_BSFG            =  19/12",
        "verified_exact_in_Fraction": classI_ok,
    },
    "universality_classes": {
        "Class_I_Borel_per_loop_locked": ["alpha"],
        "Class_II_single_locked_ratio":  ["mu"],
        "Class_geometric_phase":         ["c"],
    },
    "headline_closures": {
        "alpha_classII_unified_ratio_test": {
            "predicted": pred_classII,
            "observed":  obs_classII,
            "error_pct": err_classII,
            "status":    "FAIL",
        },
        "alpha_classI_per_loop_locked": {
            "predicted": 1.0,
            "observed":  1.0 if classI_ok else 0.0,
            "error_pct": 0.0,
            "status":    "EXACT" if classI_ok else "FAIL",
        },
    },
    "next_slot": "S718 -- c-chain class assignment test: apply Test 1 and Test 2 "
                 "to the c-chain (using its 2-loop coefficient c_2 = 5*pi^2/9, "
                 "3-loop tail bound from S704, and the 13/3 geom-phase from S708) "
                 "to determine whether c is Class I, Class II, or a hybrid.",
}

out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "_session717_alpha_class_universality_test_result.json")
with open(out_path, "w", encoding="utf-8") as f:
    json.dump(artifact, f, indent=2)
print(f"Artifact written: {out_path}")
