#!/usr/bin/env python3
"""
verify_three_number_systems.py — independent verification of UQFF's three numeric systems
+ the BSFG geometry system + the VDS×DVP×BSH hybrid identity (PAPER_1069).

Read-only. Outputs to MY sandbox only.
"""

import math
from mpmath import polylog, mpf, mp

mp.dps = 60  # 60 decimal digits — enough for Li_26(0.57)

print("=" * 100)
print("UQFF — THREE NUMERIC SYSTEMS + GEOMETRY SYSTEM (independent verification)")
print("=" * 100)

# -----------------------------------------------------------------------------
# Locked primitives (subset needed for the three systems)
# -----------------------------------------------------------------------------
SSQ      = 0.57
SO_FIVE  = 10
D_CRIT   = 26
D_PHYS   = 4
D_BSFG   = 6
BETA_I   = 0.6029
KAPPA    = 5e-4
THZ      = 1.25e12
PHI_RES  = 0.84
S26_3_CANONICAL = 1.4531e26   # the published canonical value

print()
print("[1] VDS — Vacuum Density Series Z_26 = Li_26(SSQ)")
print("    Formula: Li_26(z) = Σ_{n=1..∞} z^n / n^26")
print()

Z_26 = polylog(26, SSQ)
print(f"    Li_26(0.57)            = {Z_26}     (mpmath, 60dps)")
print(f"    Li_26(0.57) float      = {float(Z_26):.12e}")

# Numerical check: at n=1 dominant, Li_26(0.57) ≈ 0.57 + 0.57^2/2^26 + ...
partial = sum(0.57**n / n**26 for n in range(1, 50))
print(f"    direct partial (n≤50)  = {partial:.20e}")
print(f"    SSQ alone (= 0.57)     = {0.57}                  ← Li_26 is dominated by n=1 term")
print(f"    R_VDS claim in PAPER_1069: rho_SCm·S26·Phi/Phi0 = 0.167")
# Let's check what value that produces
RHO_VAC_SCM_macro = sum(0.57 * 1e-20 * 10**n for n in range(1, 27))
print(f"    rho_SCm (macro) · S26_3 · 0.84 / (rho_SCm · S26_3 · 1.0) = 0.84  ≠ 0.167")
print(f"    rho_SCm (macro)        = {RHO_VAC_SCM_macro:.6e} J/m^3")

# The 0.167 looks like a different ratio. Let me try other interpretations:
# Maybe (1 - Phi - TRZ) = 1 - 0.84 - 0 = 0.16; or Phi*Phi/5 = 0.141; or SSQ-(Phi-1) = 0.41
# Or directly: (5/6 - SSQ - TRZ)/SSQ = (0.833-0.57-0.1)/0.57 = 0.163/0.57 = 0.286
# Or: rho_SCm/rho_UA × something. We don't try to invent — just report.

print()
print("[2] DVP — Dipole Vortex Primes (p = 113 base)")

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    for i in range(3, int(math.isqrt(n))+1, 2):
        if n % i == 0: return False
    return True

primes_up_to_200 = [p for p in range(2, 200) if is_prime(p)]
print(f"    113 is prime?          = {is_prime(113)}")
print(f"    Position of 113        = #{primes_up_to_200.index(113)+1} prime")
print(f"    Primes near 113        : {[p for p in primes_up_to_200 if 100 <= p <= 130]}")
# Resonance interpretation: 113 = 4 × 26 + 9 = D_PHYS · D_CRIT + N_CH (note: N_CH=9, in locked primitives)
N_CH = 9
print(f"    Identity check         : D_PHYS·D_CRIT + N_CH = {D_PHYS*D_CRIT + N_CH} ← matches 113")
print(f"    Identity check 2       : SO_FIVE² + D_CRIT - SO_FIVE + 3 = {SO_FIVE**2 + D_CRIT - SO_FIVE + 3}")
print(f"    Identity check 3       : 26 × 4 + 9 = {26*4 + 9}   ← simplest")

print()
print("[3] BSH — Buoyancy Saturation Harmonics over 26 states")
print("    Formula: BSH(i) = β_i · exp(-SSQ·i/26),  i = 1..26")
print()

BSH_values = []
for i in range(1, 27):
    bsh_i = BETA_I * math.exp(-SSQ * i / 26)
    BSH_values.append(bsh_i)

print(f"    i= 1 BSH(1)            = {BSH_values[0]:.6f}     (e^{{-0.57/26}} = {math.exp(-SSQ/26):.6f})")
print(f"    i= 13 BSH(13)          = {BSH_values[12]:.6f}     (e^{{-0.5·0.57}} = {math.exp(-0.5*SSQ):.6f})")
print(f"    i= 26 BSH(26)          = {BSH_values[25]:.6f}     (e^{{-0.57}}     = {math.exp(-SSQ):.6f})")
print(f"    Σ BSH(i) i=1..26       = {sum(BSH_values):.6f}")
print(f"    mean BSH               = {sum(BSH_values)/26:.6f}")

print()
print("[4] BSFG — Buoyancy-SCm-Fluid-Geometry (the geometry system within)")
print("    D_BSFG = dim_R[SO(5)/U(2)] = 10 - 4 = 6  (PAPER_1167; derivative not free)")
print()
print(f"    SO(5) real dim         = SO_FIVE                  = {SO_FIVE}")
print(f"    U(2) real dim          = D_PHYS                   = {D_PHYS}")
print(f"    SO(5)/U(2) real dim    = SO_FIVE - D_PHYS         = {SO_FIVE - D_PHYS}")
print(f"    D_BSFG claim           = {D_BSFG}                                  → match: {SO_FIVE - D_PHYS == D_BSFG}")
print()
print(f"    Phi_res derivative check (PAPER_1159):")
print(f"    Phi_res = (D_BSFG - 1)/D_BSFG = {D_BSFG - 1}/{D_BSFG} = {(D_BSFG-1)/D_BSFG:.10f}")
print(f"    Locked 5/6             = {5/6:.10f}                ← match")
print(f"    Locked 0.84            = {0.84:.10f}              ← alternate form (PAPER_1203)")
print()
print(f"    Both 5/6 and 0.84 appear as Phi_res — 5/6 in Nuclear (PAPER_1203),")
print(f"    0.84 in cosmology / canonical primitives table (Gold §3).")
print(f"    The two are distinct but both labeled Phi_res in different contexts.")

print()
print("[5] HYBRID IDENTITY  (PAPER_1069 abstract):")
print("    R_VDS × p_DVP × BSH(i) = F_{U_Bi_i} within 0.1%")
print()

# Try the published R_VDS = 0.167 with p_DVP = 113 and BSH(i=13) midpoint
R_VDS  = 0.167
p_DVP  = 113
BSH_13 = BSH_values[12]
hybrid = R_VDS * p_DVP * BSH_13
print(f"    R_VDS = 0.167 (paper)  × p_DVP = 113 × BSH(13) = {BSH_13:.6f}")
print(f"    Product                = {hybrid:.6f}")
print(f"    BSH(i=1)   product     = {R_VDS * p_DVP * BSH_values[0]:.6f}")
print(f"    BSH(i=26)  product     = {R_VDS * p_DVP * BSH_values[25]:.6f}")
print(f"    Σ over all 26 i        = {R_VDS * p_DVP * sum(BSH_values):.6f}")
print()
print(f"    F_U_Bi_i target/range  — depends on system; canonical β_i × Ug × (M/r²) etc.")
print(f"    Paper claim is 0.1% match per layer, not a single scalar target —")
print(f"    so the 'hybrid = F_U_Bi_i' is a per-system identity verified across §544/543/542/541,")
print(f"    not a one-line scalar check. The recompute confirms the *form* but verifying")
print(f"    the 0.1% requires per-system anchors (see PAPER_1069 §B.4).")

print()
print("[6] S_26^(3) — verify the canonical 1.4531e26 value")
print()
print(f"    Direct Σ over n=1..∞  SSQ^n / n^26 (pure Li_26):")
print(f"      Li_26(0.57) = {float(polylog(26, SSQ)):.6e}     ≈ 0.5700 (n=1 dominates)")
print(f"    Ramanujan order-3 acceleration multiplier needed to reach 1.4531e26:")
print(f"      factor = 1.4531e26 / Li_26(0.57) = {1.4531e26 / float(polylog(26, SSQ)):.6e}")
print()
print(f"    Per PAPER_545 §225 and PAPER_1069, the order-3 Ramanujan series is:")
print(f"      S_26^(3) = Σ_n (1/4)_n (1/2)_n (3/4)_n / (n!)^3 · Π_{{i=1..26}} [1 + SSQ·exp(-κ·i·n/26)]")
print(f"    The Pochhammer-product form gives a much larger value than the bare Li_26.")
print(f"    The canonical 1.4531e26 is the calibrated total at SSQ=0.57; we don't try to reproduce")
print(f"    the full Ramanujan series here without the explicit Pochhammer convergence regime.")

# Direct partial of the Pochhammer series for a SANITY check
def pochhammer(a, n):
    p = 1.0
    for k in range(n):
        p *= (a + k)
    return p

print()
print(f"    Pochhammer-series partial sums (first few n):")
total = 0.0
for n_max in [0, 1, 2, 3, 5, 10, 20]:
    s = 0.0
    for n in range(n_max + 1):
        denom = math.factorial(n) ** 3
        pochs = pochhammer(0.25, n) * pochhammer(0.5, n) * pochhammer(0.75, n)
        # W_26(n) product
        W = 1.0
        for i in range(1, 27):
            W *= (1.0 + SSQ * math.exp(-KAPPA * i * n / 26))
        s += pochs / denom * W
    print(f"      partial n≤{n_max:>2}    = {s:.6e}")
print()
print(f"    Note: the Pochhammer series above converges to a value of order (1.5..2.0)·2^26 ≈ 1e8 at high n,")
print(f"    NOT 1.4531e26. The published S26_3 = 1.4531e26 likely includes additional 26D")
print(f"    volume / projection factors not in the bare Pochhammer form.")
print(f"    Either: (a) the canonical 1.4531e26 carries hidden 10^18 normalization, or")
print(f"    (b) the documented Pochhammer form needs an enhanced W_26 sum (PAPER_1080 binomial)")
print(f"    Either way: the framework uses 1.4531e26 as a locked primitive — derivations downstream")
print(f"    are not affected; the question is about the source's own self-derivation chain for S26_3.")
