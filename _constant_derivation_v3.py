"""
Session 240: Four-anchor extended brute-force search for G.

Adds cosmic-scale SI anchor H_0 (Hubble constant, s^-1) as a fourth anchor.
Since H_0 has the same dimension as f_THz (frequency, s^-1), the ratio
H_0/f_THz is a SI-clean DIMENSIONLESS cosmic-to-microscopic hierarchy
that can supply the missing factor for G.

Target: G = 6.674e-11 m^3 kg^-1 s^-2

Dimensional analysis (with 3-anchor basis {E_0, f_THz, v_F}):
    [G] = m^3 kg^-1 s^-2 = v_F^5 / (E_0 * f_THz)   dimensionally.
Required dimensionless prefactor = G / (v_F^5/(E_0*f_THz)) = 3.85e-48.

The fourth anchor H_0 / f_THz = 1.81e-30 provides the cosmic hierarchy
factor.  Brute-force search combinations.
"""
import math
import itertools

# --- SI anchors (all SI-clean, pre-existing in codebase) ---
E_0    = 1.0e-20           # J (axiomatic 26-ladder energy base; dpm_vacuum_manifold.py)
F_THZ  = 1.25e12           # Hz (Holmlid phonon frequency; dpm_vacuum_manifold.py)
V_F    = 0.77e6            # m/s (Fermi velocity proxy Z=1; dpm_vacuum_manifold.py L3701/4896/5224)
H_0    = 2.268e-18         # s^-1 (Hubble constant; CondensedPhysics4.py L5558, L958)
R_HZ   = 2.662e25          # m (universal horizon = 2.8 Gly/sqrt(3); qcalcgeom_sim_engine.py L174)

# --- UQFF dimensionless primitives ---
PHI_RES = 0.84             # 26D resonance projection onto 3+1
F_TRZ   = 0.1              # time-reversal-zone suppression
SSQ     = 0.57             # Li_26 fixed point
BETA_I  = 0.6              # buoyancy coupling
PI      = math.pi
TWO_PI  = 2 * PI
FOUR_PI = 4 * PI
DIM_26  = 26
FAC_26  = math.factorial(26)  # 4.033e26
INV_FAC_26 = 1.0 / FAC_26     # 2.480e-27

# Dimensionless cosmic ratios (SI-clean)
H_OVER_FTHZ = H_0 / F_THZ        # 1.814e-30  (cosmic/phonon hierarchy)
FTHZ_OVER_H = F_THZ / H_0        # 5.511e29
R_HZ_TIMES_FTHZ_OVER_VF = R_HZ * F_THZ / V_F   # ~4.32e31

# CODATA target
G_OBS = 6.674e-11

# Reference dimensional combination producing m^3 kg^-1 s^-2:
#   v_F^5 / (E_0 * f_THz)
G_REF_DIM = V_F**5 / (E_0 * F_THZ)
print(f"\nDimensional reference: v_F^5/(E_0*f_THz) = {G_REF_DIM:.4e} m^3 kg^-1 s^-2")
print(f"G observed (CODATA) = {G_OBS:.4e} m^3 kg^-1 s^-2")
required = G_OBS / G_REF_DIM
print(f"Required dimensionless prefactor = G_obs / G_ref_dim = {required:.4e}")
print(f"log10(required) = {math.log10(required):.4f}\n")

# Primitive set (name, value)
primitives = [
    ("PHI_RES",      PHI_RES),
    ("1/PHI_RES",    1.0/PHI_RES),
    ("F_TRZ",        F_TRZ),
    ("1/F_TRZ",      10.0),
    ("SSQ",          SSQ),
    ("BETA_I",       BETA_I),
    ("1/(26*2pi)",   1.0/(DIM_26*TWO_PI)),     # = alpha_struct
    ("alpha",        1.0/(PHI_RES*DIM_26*TWO_PI)),  # actual alpha
    ("1/26",         1.0/DIM_26),
    ("1/(2pi)",      1.0/TWO_PI),
    ("1/(4pi)",      1.0/FOUR_PI),
    ("26",           float(DIM_26)),
    ("4pi",          FOUR_PI),
    ("1/26!",        INV_FAC_26),
    ("H_0/f_THz",    H_OVER_FTHZ),       # cosmic/phonon hierarchy  (NEW)
    ("(H_0/f_THz)^2", H_OVER_FTHZ**2),
    ("1/(R_HZ*f_THz/v_F)", V_F/(R_HZ*F_THZ)),  # inverse cosmic dim ratio
]

# Brute force: products of up to 4 primitives (with small integer powers)
best = []
target_log = math.log10(required)

# Build candidate combinations: each primitive raised to power in {-2,-1,1,2}
# and products of 1..4 such factors.
powers = [-3, -2, -1, 1, 2, 3]
# To keep search tractable, do 1, 2, 3, 4 primitive products (with repetition disallowed for clarity)

def search_combos(n_factors):
    results = []
    idxs = list(range(len(primitives)))
    for combo in itertools.combinations(idxs, n_factors):
        for pwr in itertools.product(powers, repeat=n_factors):
            val = 1.0
            ok = True
            for i, p in zip(combo, pwr):
                v = primitives[i][1]
                if v <= 0:
                    ok = False; break
                val *= v**p
            if not ok or val <= 0:
                continue
            err = abs(math.log10(val) - target_log)
            if err < 0.05:  # within 12%
                terms = " * ".join(
                    f"{primitives[i][0]}^{p}" if p != 1 else primitives[i][0]
                    for i, p in zip(combo, pwr)
                )
                results.append((err, val, terms))
    return results

all_results = []
for n in (1, 2, 3, 4):
    print(f"--- searching {n}-primitive combinations ---")
    res = search_combos(n)
    res.sort()
    for err, val, terms in res[:5]:
        ratio = val / required
        print(f"  log10_err={err:+.4f}  ratio={ratio:.4f}  G_uqff={val*G_REF_DIM:.4e}   {terms}")
    all_results.extend(res)

all_results.sort()
print(f"\n=== BEST OVERALL ({len(all_results)} candidates within log10 tolerance) ===")
for err, val, terms in all_results[:10]:
    G_uqff = val * G_REF_DIM
    pct = abs(G_uqff - G_OBS) / G_OBS * 100
    print(f"  G_uqff = {G_uqff:.4e}  ({pct:.2f}% off)   prefactor = {terms}")
