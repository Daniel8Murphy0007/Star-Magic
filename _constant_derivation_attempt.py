# -*- coding: utf-8 -*-
"""
_constant_derivation_attempt.py

Honest, exhaustive attempt to derive h, alpha, c, G from the canonical UQFF
primitives that already exist in dpm_vacuum_manifold.py and the related
Quantum Chain code. NO curve-fitting knobs are introduced. NO new free
parameters. Only the primitives that the user has already built up are used.

Strategy:
  1. Pull ALL canonical primitives from dpm_vacuum_manifold.py.
  2. Identify which are SI-clean (independent of h/c/G) and which are circular
     (built from h/c/G).
  3. For each target constant, try the formulas already encoded in CP4
     #177-#180 with the canonical primitives substituted in.
  4. Also try a brute-force algebraic recombination: any ratio/product of
     SI-clean primitives that has the right dimensions and lands within
     factor 100 of the observed value.
  5. Report concrete numbers. Whatever they are.
"""

import math
import itertools

# ============================================================================
# 1. CANONICAL UQFF PRIMITIVES (from dpm_vacuum_manifold.py)
# ============================================================================

# --- SI-clean (no h/c/G dependency in their definitions) ---
SSQ        = 0.57                     # [SSq] dimensionless
KAPPA_DAY  = 5.0e-4                   # 1/day
KAPPA_S    = KAPPA_DAY / 86400.0      # 1/s  = 5.787e-9
THZ_PHONON = 1.25e12                  # Hz   (Holmlid measured)
OMEGA_S    = 2.5e-6                   # rad/s (stellar rotation)
BETA_I     = 0.6                      # buoyancy coupling
LAMBDA_I   = 1.0
F_TRZ      = 0.1
PHI_RES    = 0.84

# --- Quantum-chain energy hierarchy (E0 = 1e-20 J is a free axiom) ---
E0_J       = 1.0e-20                  # J  (axiom: base of 26-level magnitude ladder)
N_LEVELS   = 26

# --- Dimensionless 26D Ramanujan factor (calibrated to Holmlid 630 eV in
#     dpm_vacuum_manifold.py via scaling_factor; here we treat the *raw*
#     Li_26(0.57) polylog value, which is honest and circularity-free) ---
from mpmath import polylog
LI26_SSQ_RAW = float(polylog(26, SSQ))  # ~ 0.5700... per dpm_vacuum_manifold

# --- vacuum energy densities from quantum chain (J/m^3) ---
def derive_rho(f_SCm, n_levels=26, V=1.0):
    return sum(f_SCm * (E0_J * 10**n) for n in range(1, n_levels + 1)) / V
RHO_VAC_SCM = derive_rho(SSQ)               # 0.57 * sum(1e-19..1e6) ~ 5.7e5 J/m^3
RHO_VAC_UA  = derive_rho(SSQ * 10)          # 10x higher

# --- Circular primitives (defined via h or c, NOT usable as inputs) ---
# E_phonon = h * THZ_PHONON         <- circular for h
# v_SCm    = 0.99 * c               <- circular for c
# S26_3    = 1.4531e26 (tuned to 630 eV w/ h*THZ; circular)
# These are EXCLUDED from the derivation attempt.

# ============================================================================
# 2. OBSERVED VALUES (CODATA / PDG)
# ============================================================================
H_OBS     = 6.62607015e-34   # J*s
HBAR_OBS  = H_OBS / (2 * math.pi)
ALPHA_OBS = 7.2973525693e-3  # ~1/137.036
C_OBS     = 2.99792458e8     # m/s
G_OBS     = 6.6743e-11       # m^3 kg^-1 s^-2

# ============================================================================
# 3. ATTEMPT EACH DERIVATION USING ENCODED FORMULAS
# ============================================================================

print("=" * 78)
print("UQFF Constant Derivation Attempt -- canonical primitives, no fit knobs")
print("=" * 78)
print(f"\nSI-clean primitives:")
print(f"  SSQ              = {SSQ}")
print(f"  KAPPA_S          = {KAPPA_S:.4e} 1/s")
print(f"  THZ_PHONON       = {THZ_PHONON:.4e} Hz")
print(f"  OMEGA_S          = {OMEGA_S:.4e} rad/s")
print(f"  E0_J             = {E0_J:.4e} J  (axiomatic base)")
print(f"  RHO_VAC_SCM      = {RHO_VAC_SCM:.4e} J/m^3")
print(f"  RHO_VAC_UA       = {RHO_VAC_UA:.4e} J/m^3")
print(f"  Li_26(SSQ) raw   = {LI26_SSQ_RAW:.6f}")
print(f"  BETA_I           = {BETA_I}")
print(f"  F_TRZ            = {F_TRZ}")
print(f"  PHI_RES          = {PHI_RES}")

# ----------------------------------------------------------------------------
# 3a. Encoded formula attempts (CP4 #177-#180 source-canonical forms)
# ----------------------------------------------------------------------------

def report(label, derived, observed, units=""):
    if observed == 0 or derived == 0:
        ratio_str = "(div0)"
        log_off = "?"
    else:
        ratio = derived / observed
        ratio_str = f"{ratio:.3e}"
        log_off = f"{math.log10(abs(ratio)):+.2f}"
    print(f"  {label:<40} = {derived:+.4e} {units}  (obs={observed:.4e}, ratio={ratio_str}, log10_off={log_off})")


print("\n" + "-" * 78)
print("[1] PLANCK h  --  CP4 #177 source-canonical form")
print("    h = 2*pi * kappa * rho * Grind_opp / r^2")
print("-" * 78)
# Use canonical UQFF primitives (NOT Grok's stated junk parameters).
# kappa: KAPPA_S [1/s].  rho: RHO_VAC_SCM [J/m^3].  Grind: dimensionless ratio.
# r: must be a length scale.  Try Bohr (5.29e-11) and Planck (1.616e-35).
for r_label, r_m in [("Bohr (5.29e-11 m)", 5.29e-11),
                     ("Planck (1.616e-35 m)", 1.616e-35),
                     ("Holmlid wavelength c/f_THz", 2.998e8 / THZ_PHONON)]:
    # Grind_opp from canonical: omega_CW=THZ_PHONON, omega_CC=THZ_PHONON,
    # SCm=1, UA'=10, exp(-E/v)=exp(-E0/v_init).  v_init unknown; try OMEGA_S and 1.
    for v_init in [1.0, OMEGA_S, KAPPA_S]:
        grind = THZ_PHONON * 1 - THZ_PHONON * 10 * math.exp(-E0_J / v_init)
        h_try = 2 * math.pi * KAPPA_S * RHO_VAC_SCM * grind / (r_m * r_m)
        report(f"r={r_label}, v_init={v_init:.2e}", h_try, H_OBS, "J*s")

print("\n  ALTERNATE: try h = E0 * (2*pi/THZ_PHONON) * Li26  (h ~ E*T form)")
h_alt1 = E0_J * (2 * math.pi / THZ_PHONON) * LI26_SSQ_RAW
report("  E0*(2pi/f)*Li26", h_alt1, H_OBS, "J*s")
print("\n  ALTERNATE: h = E0 / THZ_PHONON  (Planck relation E=hf inverted)")
h_alt2 = E0_J / THZ_PHONON
report("  E0/f_THz", h_alt2, H_OBS, "J*s")

# ----------------------------------------------------------------------------
# 3b. Fine-structure alpha
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("[2] FINE STRUCTURE alpha  --  CP4 #178 source-canonical form")
print("    alpha = 2*kappa*rho*Grind^2 * r^24 * Partition / (3*sqrt(g*SCm/UA))")
print("-" * 78)
# Use canonical primitives. kappa=KAPPA_S, rho=RHO_VAC_SCM, Grind=ratio,
# r=Bohr (only sane length), Partition=1e5 (Grok), g=ratio, SCm/UA=0.1.
grind_dim = 1e-3   # Grok's Grind from source (dimensionless)
for r_label, r_m in [("Bohr 5.29e-11", 5.29e-11)]:
    for partition in [1.0, 1e5, LI26_SSQ_RAW, 1.0 / SSQ]:
        for g_couple in [SSQ, BETA_I, RHO_VAC_SCM/RHO_VAC_UA, 1.0]:
            scm_ua = 0.1   # canonical UQFF ratio
            num = 2 * KAPPA_S * RHO_VAC_SCM * grind_dim**2 * r_m**24 * partition
            den = 3 * math.sqrt(abs(g_couple * scm_ua))
            a_try = num / den
            report(f"P={partition:.2e}, g={g_couple:.2e}", a_try, ALPHA_OBS)

print("\n  ALTERNATE: alpha as pure dimensionless ratio of UQFF primitives:")
# Brute-force: search products/ratios of {SSQ, BETA_I, F_TRZ, PHI_RES,
# LI26_SSQ_RAW, 1/26, 1/137, KAPPA_DAY, ...}
dimless = {
    'SSQ': SSQ, 'BETA_I': BETA_I, 'F_TRZ': F_TRZ, 'PHI_RES': PHI_RES,
    'Li26': LI26_SSQ_RAW, '1/26': 1/26, '1/26!': 1/math.factorial(26),
    'SSQ^26': SSQ**26, 'log(SSQ)': abs(math.log(SSQ)),
    '2pi': 2*math.pi, 'pi': math.pi,
}
best_alpha = []
for n in [1, 2]:
    for combo in itertools.combinations(dimless.items(), n):
        keys = [k for k, _ in combo]; vals = [v for _, v in combo]
        for signs in itertools.product([1, -1], repeat=n):
            try:
                p = 1.0
                for v, s in zip(vals, signs):
                    p *= v**s
                if 1e-4 < p < 1e-1:
                    log_off = abs(math.log10(p / ALPHA_OBS))
                    if log_off < 1.0:
                        best_alpha.append((log_off, p, keys, signs))
            except (ValueError, ZeroDivisionError, OverflowError):
                pass
best_alpha.sort()
for log_off, p, keys, signs in best_alpha[:6]:
    print(f"  {keys} signs={signs} -> {p:.4e}  (alpha_obs={ALPHA_OBS:.4e}, log_off={log_off:+.3f})")
if not best_alpha:
    print("  No 1- or 2-primitive dimensionless combo lands within 1 order of alpha_obs.")

# ----------------------------------------------------------------------------
# 3c. Speed of light c
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("[3] SPEED OF LIGHT c  --  CP4 #179: c = sqrt(g * SCm/UA)")
print("-" * 78)
print("  c is an UQFF AXIOM (v_init = c). Any 'derivation' is circular by")
print("  construction. Showing this concretely:")
for g_label, g in [("SSQ", SSQ), ("BETA_I", BETA_I), ("RHO_VAC_SCM/UA", RHO_VAC_SCM/RHO_VAC_UA),
                    ("RHO_VAC_SCM [J/m^3]", RHO_VAC_SCM), ("1.0", 1.0)]:
    c_try = math.sqrt(abs(g * 0.1))
    report(f"g={g_label}", c_try, C_OBS, "m/s")

# A more interesting attempt: c = (E0/RHO_VAC_SCM)^something^^^ for length/time
# units. RHO_VAC_SCM is J/m^3, E0 is J, so E0/RHO_VAC_SCM is m^3.
vol_implied = E0_J / RHO_VAC_SCM
length_implied = vol_implied**(1/3)
print(f"\n  E0 / RHO_VAC_SCM      = {vol_implied:.4e} m^3 (volume)")
print(f"  length scale          = {length_implied:.4e} m")
# If we propose c = length / time using KAPPA_S as the time:
c_try = length_implied / KAPPA_S
report("length/KAPPA_S", c_try, C_OBS, "m/s")
c_try2 = length_implied * THZ_PHONON
report("length*THZ_PHONON", c_try2, C_OBS, "m/s")

# ----------------------------------------------------------------------------
# 3d. Newton G
# ----------------------------------------------------------------------------
print("\n" + "-" * 78)
print("[4] NEWTON G  --  CP4 #180: four candidate forms")
print("-" * 78)
g_couple = RHO_VAC_SCM / RHO_VAC_UA   # = 0.1 (canonical UQFF coupling)
rho_void_kg_per_m3 = RHO_VAC_SCM / C_OBS**2  # ~ 6.3e-12 kg/m^3 (mass-equivalent)
print(f"  g_couple (SCm/UA)    = {g_couple}")
print(f"  rho_void (kg/m^3)    = {rho_void_kg_per_m3:.4e}")
G_triad   = g_couple * 0.1
G_buoyant = g_couple / (4 * math.pi * rho_void_kg_per_m3)
G_full    = g_couple * math.exp(-1e-3) / (4 * math.pi * rho_void_kg_per_m3)
G_gauss   = g_couple / (rho_void_kg_per_m3 * 1e16 / 92e9)
report("G_triad",    G_triad,    G_OBS, "m^3/kg/s^2")
report("G_buoyant",  G_buoyant,  G_OBS, "m^3/kg/s^2")
report("G_full",     G_full,     G_OBS, "m^3/kg/s^2")
report("G_gauss",    G_gauss,    G_OBS, "m^3/kg/s^2")

# ============================================================================
# 4. SUMMARY
# ============================================================================
print("\n" + "=" * 78)
print("HONEST RESULT SUMMARY")
print("=" * 78)
print("""
The canonical UQFF primitives in dpm_vacuum_manifold.py are:
  SSQ=0.57, KAPPA=5e-4/day, THZ=1.25e12, OMEGA_S=2.5e-6,
  E0=1e-20 J (axiom), RHO_VAC_SCM (derived from E0 and SSQ).

WITHOUT introducing curve-fit knobs, the encoded CP4 #177-#180 formulas
do NOT reproduce h, alpha, c, or G to any reasonable accuracy. The two
'closest' approaches (E0/THZ for h; certain dimensionless ratios for
alpha) only work because E0 is an axiomatic free parameter chosen to
make the 26-level ladder span the relevant physics range -- it is itself
calibrated FROM observation, not derived.

CONCRETE conclusion: the 8.5M-line codebase does not currently contain a
non-circular derivation of these four constants. The pieces it contains
are:
  - 26-level energy ladder (E0=1e-20 chosen as axiom)
  - 26D Ramanujan amplification S26_3 (tuned to Holmlid 630 eV)
  - rho_vac (derived from the ladder + SSq)
  - kappa, OMEGA_S, BETA_I (calibrated from astrophysical observations)
  - dimensionless ratios (SSq, F_TRZ, BETA_I)

The honest path: at least one independent SI anchor is required. The
existing E0 = 1e-20 J was that anchor for the energy hierarchy. To
derive h, alpha, c, G without circularity, you need to either:
  (a) pin one of {h, c, G} from observation and derive the other three,
  (b) introduce a new independent UQFF observable (e.g. measured
      vacuum impedance, a measured DPM frequency at a calibrated lab
      energy) that fixes the absolute scale.

Both (a) and (b) are honest. Neither gets you all four for free. This is
not a bug to fix -- it is the actual structure of the framework. UQFF
constrains *relationships* between fundamental constants. It does not
appear to derive their absolute SI values from purely internal axioms.
""")
