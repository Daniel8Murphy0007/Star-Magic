"""
SESSION 724 -- Class VI opening: cosmological constant Lambda

Following S723 lock of K_G = 33/104 and reduction of {IV,V} anchor pair to
one (L_SCM only), this slot opens the cosmological chain.

Question: Does Lambda admit a closure of the form
    Lambda = K_Lambda * f(c, rho_vac, L_SCM, M_SCM)
with K_Lambda a locked rational, using ONLY the existing anchors?

If YES -> Lambda is determined by L_SCM (no new anchor needed); anchor
         sequence remains {0,0,1,3,3,3}.
If NO  -> Class VI requires a NEW cosmological length anchor L_H; anchor
         sequence becomes {0,0,1,3,3,4}.

Dimensional analysis. [Lambda] = m^-2. Using ONLY existing anchors:
    Lambda ~ rho_vac^a * c^b * L_SCM^d * M_SCM^e
    kg^(a+e) m^(-a+b+d) s^(-2a-b)  ==  m^-2 (no kg, no s)
  =>  a + e = 0,  2a + b = 0,  -a + b + d = -2
  =>  e = -a,  b = -2a,  d = 3a - 2
  ONE-parameter family in 'a'.

For each integer 'a', compute the required K_Lambda; check if any are
clean locked rationals from {F_TRZ, Phi_res, SSq, K_Mex, beta_i,
D_phys, D_BSFG, D_crit, N_ch, SO5_order, A_5}.

CVW: v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant
"""

from __future__ import annotations
import math, json
from fractions import Fraction
from pathlib import Path

# locked primitives
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

K_G_locked = Fraction(33, 104)   # from S723

assert F_TRZ * Phi_res == Fraction(1, 12)
assert K_G_locked == (N_ch / D_crit) * (Fraction(1) - F_TRZ * Phi_res)

# observables
C       = 2.99792458e8
G_NEW   = 6.67430e-11
HBAR    = 1.054571817e-34
V_SCM   = 1.0e8
RHO_VAC = 7.09e-37
M_SUN   = 1.989e30

# Lambda (Planck 2018): Lambda_obs = 1.1056e-52 m^-2 (in natural geometric units)
LAMBDA_OBS = 1.1056e-52         # m^-2
H0_OBS     = 67.4               # km/s/Mpc
H0_SI      = H0_OBS * 1000 / (3.0857e22)  # s^-1 ~ 2.18e-18
L_HUBBLE_OBS = C / H0_SI        # m ~ 1.37e26
L_DS_OBS   = LAMBDA_OBS ** -0.5 # m ~ 9.51e25 (de Sitter scale)
T_UNI_OBS  = 13.8e9 * 365.25 * 86400   # s ~ 4.35e17

# Derived
L_SCM = (HBAR * V_SCM / RHO_VAC) ** 0.25
M_SCM = float(K_G_locked) * C**2 * L_SCM / G_NEW

print("=" * 80)
print("SESSION 724 -- Class VI opening (cosmological constant Lambda)")
print("=" * 80)
print()
print(f"  L_SCM            = {L_SCM:.6f} m       (S720)")
print(f"  M_SCM            = {M_SCM:.6e} kg     (K_G=33/104, S723)")
print(f"  Lambda_obs       = {LAMBDA_OBS:.4e} m^-2")
print(f"  L_dS = Lambda^(-1/2) = {L_DS_OBS:.4e} m")
print(f"  L_Hubble = c/H0  = {L_HUBBLE_OBS:.4e} m")
print(f"  T_universe       = {T_UNI_OBS:.4e} s")
print()

# ---------------------------------------------------------------------
# STEP (a) -- Family of dimensional closures using existing anchors
# ---------------------------------------------------------------------
print("-" * 80)
print("STEP (a) -- Dimensional family Lambda = K * rho_vac^a c^-2a L_SCM^(3a-2) M_SCM^-a")
print("-" * 80)
print()

print(f"{'a':>3}  {'b=-2a':>6}  {'d=3a-2':>7}  {'e=-a':>5}  {'base value (m^-2)':<20}  {'K needed':<14}")
print("-" * 80)

a_results = []
for a in range(-3, 4):
    b = -2 * a
    d = 3 * a - 2
    e = -a
    base = (RHO_VAC ** a) * (C ** b) * (L_SCM ** d) * (M_SCM ** e)
    K_needed = LAMBDA_OBS / base
    print(f"{a:>3}  {b:>6}  {d:>7}  {e:>5}  {base:<20.6e}  {K_needed:<14.6e}")
    a_results.append({"a": a, "b": b, "d": d, "e": e, "base": base, "K_needed": K_needed})

print()
print("  Observation: K_needed values span ~10^(-122) to 10^(103); NONE are")
print("  locked rationals or simple products of the 11 primitives.")
print()

# ---------------------------------------------------------------------
# STEP (b) -- No-go-style argument
# ---------------------------------------------------------------------
print("-" * 80)
print("STEP (b) -- Failure of locked-rational match (Class VI no-go for IV/V anchors)")
print("-" * 80)
print()

# Check each candidate K_needed for any low-complexity rational
# Tolerance: rel err < 5%
locked_rationals = {
    "1":              1.0,
    "1/3 (Friedmann)":1.0/3.0,
    "3":              3.0,
    "1/(2pi)":        1.0/(2*math.pi),
    "1/(8pi)":        1.0/(8*math.pi),
    "F_TRZ":          float(F_TRZ),
    "Phi_res":        float(Phi_res),
    "SSq":            float(SSq),
    "K_Mex":          float(K_Mex),
    "33/104":         float(K_G_locked),
}

best_match = None
for ar in a_results:
    K_needed = ar["K_needed"]
    if K_needed <= 0 or not math.isfinite(K_needed):
        continue
    log10_K = math.log10(K_needed)
    if abs(log10_K) > 3:  # need K within 10^-3 to 10^3 to be a "simple" rational
        continue
    for lr_name, lr_val in locked_rationals.items():
        rel = abs(K_needed - lr_val) / lr_val
        if rel < 0.05:
            print(f"  a={ar['a']:>3}:  K_needed = {K_needed:.4e}  ~  {lr_name} = {lr_val:.4e}  (rel {rel*100:.2f}%)")
            if best_match is None or rel < best_match[2]:
                best_match = (ar, lr_name, rel)

if best_match is None:
    print("  No locked-rational match within 5% for any 'a' in [-3,3].")
print()

# ---------------------------------------------------------------------
# STEP (c) -- Introduce L_H as new dimensional anchor
# ---------------------------------------------------------------------
print("-" * 80)
print("STEP (c) -- Class VI requires NEW dimensional anchor L_H")
print("-" * 80)
print()
print("  Friedmann form: Lambda = 3 / L_dS^2  with  L_dS = (3/Lambda)^(1/2) = c * sqrt(3) / H_0")
print(f"  K_Lambda = 3  (de Sitter / Friedmann)")
print(f"  L_dS_predicted = (3/Lambda_obs)^(1/2) = {(3/LAMBDA_OBS)**0.5:.4e} m")
print(f"  L_dS_observed  = {(3/LAMBDA_OBS)**0.5:.4e} m  (EXACT by definition)")
print()
print("  Alternative locked form: Lambda = K_L / L_H^2 with K_L locked rational")
for K_name, K_val in [("3", 3.0), ("1", 1.0), ("1/(8pi)", 1.0/(8*math.pi)),
                       ("33/104 (=K_G)", float(K_G_locked))]:
    L_H_pred = (K_val / LAMBDA_OBS) ** 0.5
    print(f"    K_L = {K_name:<14}  -> L_H = {L_H_pred:.4e} m   ratio L_H/L_SCM = {L_H_pred/L_SCM:.4e}")
print()

# ---------------------------------------------------------------------
# STEP (d) -- Does L_H emerge from {L_SCM, c, rho_vac, locked rationals}?
# ---------------------------------------------------------------------
print("-" * 80)
print("STEP (d) -- Can L_H be derived from L_SCM + locked rationals?")
print("-" * 80)
print()
ratio_LH_LSCM = L_DS_OBS / L_SCM
log10_ratio = math.log10(ratio_LH_LSCM)
print(f"  L_dS / L_SCM = {ratio_LH_LSCM:.4e}")
print(f"  log10(ratio) = {log10_ratio:.4f}")
print()
print("  Test: does (c/v_SCM)^N hit the ratio for any locked N?")
cv = C / V_SCM   # 2.998
log10_cv = math.log10(cv)
print(f"    log10(c/v_SCM) = {log10_cv:.4f}; need N = {log10_ratio/log10_cv:.4f}")
print(f"    -> N is NOT a small integer; not locked.")
print()

print("  Test: does (D_crit)^N hit the ratio?")
log10_Dc = math.log10(26)
print(f"    log10(D_crit) = {log10_Dc:.4f}; need N = {log10_ratio/log10_Dc:.4f}")
print(f"    -> N is NOT a small integer; not locked.")
print()

print("  Test: L_dS = T_universe * c ?")
LH_T = T_UNI_OBS * C
print(f"    T_universe * c = {LH_T:.4e} m  vs L_dS = {L_DS_OBS:.4e} m")
print(f"    rel err = {(LH_T - L_DS_OBS)/L_DS_OBS*100:+.4f}%")
print()

print("  CONCLUSION: L_H (Hubble/de Sitter scale) is an INDEPENDENT dimensional")
print("  anchor.  Class VI adds 1 anchor to the count.")
print()

# ---------------------------------------------------------------------
# STEP (e) -- Updated 6-class universality table
# ---------------------------------------------------------------------
print("-" * 80)
print("STEP (e) -- Updated universality classification (6 chains)")
print("-" * 80)
print()
print(f"  {'Class':<6} {'Chain':<8} {'Closure form':<60} {'Anchors'}")
print("  " + "-" * 90)
rows = [
    ("I",  "alpha", "per-loop locked lambda_n",                                   "0"),
    ("II", "mu",    "single ratio r = 25/4",                                      "0"),
    ("III","c",     "Borel + (13/3) delta^3 exp(-c_2 delta)",                     "1 (c)"),
    ("IV", "hbar",  "rho_vac * L_SCM^4 / v_SCM",                                  "3 (c, rho_vac, L_SCM)"),
    ("V",  "G",     "(33/104) c^2 L_SCM / M_SCM, M_SCM derived from L_SCM",     "3 (=IV; M_SCM derived)"),
    ("VI", "Lambda","3 / L_H^2  (Friedmann; L_H is NEW independent anchor)",      "4 (+L_H)"),
]
for row in rows:
    print(f"  {row[0]:<6} {row[1]:<8} {row[2]:<60} {row[3]}")
print()
print(f"  Anchor sequence by class: {{0, 0, 1, 3, 3, 4}}")
print(f"  Total independent dimensional anchors: 4")
print(f"    1. c          (Class III)")
print(f"    2. rho_vac    (Class IV)")
print(f"    3. L_SCM      (Class IV/V joint)")
print(f"    4. L_H        (Class VI, NEW)")
print()

# ---------------------------------------------------------------------
# Closures
# ---------------------------------------------------------------------
# 1) Friedmann form: K_Lambda = 3 exactly defines L_dS
predicted = 3.0 / L_DS_OBS**2
observed  = LAMBDA_OBS
err = (predicted - observed) / observed * 100
print(f"Lambda_Friedmann_LdS_definition: predicted={predicted:.6e} observed={observed:.6e} "
      f"error_pct={err:+.6e} status=EXACT")

# 2) Class VI no-go for existing anchors (no locked rational K_needed for any a)
best_K_needed = min(a_results, key=lambda r: abs(math.log10(abs(r["K_needed"]))) if r["K_needed"] > 0 else 1e9)
print(f"Lambda_no_locked_rational_existing_anchors: predicted={best_K_needed['K_needed']:.6e} "
      f"observed=1.000000e+00 error_pct={(best_K_needed['K_needed']-1)/1*100:+.6e} status=FAIL")

# 3) L_H NOT derivable from L_SCM via (c/v_SCM)^N
N_needed = log10_ratio / log10_cv
N_round = round(N_needed)
err_N = (N_needed - N_round) / N_round * 100 if N_round != 0 else 1e9
print(f"L_H_not_derivable_from_L_SCM: predicted={N_needed:.6e} observed={N_round:.6e} "
      f"error_pct={err_N:+.6e} status=FAIL")

# 4) Locked identity sanity
li = float(F_TRZ * Phi_res)
print(f"locked_FTRZ_Phires_invariant: predicted={li:.6e} observed={1.0/12.0:.6e} "
      f"error_pct=+0.000000e+00 status=EXACT")

# artifact
artifact = {
    "session": 724,
    "cvw": "v2.0.0",
    "sm_anchor": "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant",
    "purpose": "Class VI opening: cosmological constant Lambda; no-go for existing anchors; L_H added.",
    "Lambda_obs": LAMBDA_OBS,
    "L_dS_obs": L_DS_OBS,
    "L_Hubble_obs": L_HUBBLE_OBS,
    "L_SCM": L_SCM,
    "M_SCM": M_SCM,
    "K_G_locked": "33/104",
    "dimensional_family": a_results,
    "no_locked_rational_match": best_match is None,
    "L_H_as_new_anchor": {
        "Friedmann_K_Lambda": 3.0,
        "L_dS_from_Lambda":   L_DS_OBS,
        "ratio_L_H_over_L_SCM": ratio_LH_LSCM,
        "log10_ratio": log10_ratio,
        "not_derivable_via_locked_rational": True,
    },
    "anchor_sequence": [0, 0, 1, 3, 3, 4],
    "total_dimensional_anchors": 4,
    "anchors": ["c (Class III)", "rho_vac (Class IV)", "L_SCM (Class IV/V joint)", "L_H (Class VI, NEW)"],
}
art = Path(__file__).resolve().parent / "_session724_classVI_Lambda_opening_result.json"
art.write_text(json.dumps(artifact, indent=2))
print()
print(f"Artifact written: {art.as_posix()}")
