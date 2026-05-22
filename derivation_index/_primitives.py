"""
UQFF Locked Primitives
======================

THE ANCHOR FILE. Every entry in derivation_index/ must source its inputs
from this file ONLY. CODATA / SM measured values may appear in a
derivation file ONLY on the right-hand side of a line that begins with
the comment token `# residual:` (i.e., the comparison line, never the
input line).

Rules (binding):
  1. Values here are stored as exact rationals (fractions.Fraction) or
     exact symbolic expressions (sympy). Floats are derived on demand
     via `.float()` accessors, never stored.
  2. No primitive is ever re-fit, re-tuned, or "improved" from a
     downstream comparison. If a derivation disagrees with observation,
     the derivation is wrong, not the primitive.
  3. The constant `c` (speed of light) is the ONLY external SM constant
     allowed as an input, and only because v_UA is defined as c/3. It
     enters by reference, not by re-derivation. CODATA exact value:
     c = 299_792_458 m/s (exact by SI definition since 1983).
  4. Adding a new primitive requires a citation to the originating
     UQFF derivation file (a `_session*.py` or whitepaper).

Provenance:
  - ρ_SCm = 4·√π · 10⁻³⁷ J/m³        : _session695_4pi_coset_derivation.py
                                       PAPER_1198 (ρ_SCm Derivation)
  - ρ_UA  = 10 · ρ_SCm               : G9 axiom; _session288 (Universal
                                       Buoyancy Simultaneous Solver)
  - F_TRZ = 1 / |SO(5)| = 1/10       : G7 axiom; PAPER_1161 (F_TRZ/SO(5))
  - [SSq] = 57/100                   : PAPER_1154 (SSq 0.57 First Principles)
  - Φ_res = 5/6                      : PAPER_1159 (Φ_res)
  - K_Mex = 25/12 = Φ_res · |SO(5)| / D_phys : G1 closure; PAPER_1066
  - β_i   = 6029/10000               : G2 triangular closure; β_i ≈ 0.6029
  - v_UA  = c/3                      : SCm wave-speed axiom
  - κ     = 5×10⁻⁴ / day             : decay calibration (Grok 4, Sept 2025)
  - H_SCm ≈ 0.99                     : CALIBRATED, not yet derived from
                                       first principles. Flagged.
  - D_phys=4, D_BSFG=6, D_crit=26, N_ch=9, |SO(5)|=10 : dimensional axioms
"""

from fractions import Fraction
from sympy import sqrt, pi, Rational, Symbol


# ---------------------------------------------------------------------------
# External SM reference (allowed by rule 3 above; used only in v_UA).
# ---------------------------------------------------------------------------
C_LIGHT_M_PER_S = 299_792_458  # exact, SI definition


# ---------------------------------------------------------------------------
# Dimensional / group-theoretic axioms (integers — no floats).
# ---------------------------------------------------------------------------
D_PHYS   = 4    # physical spacetime dimensions
D_BSFG   = 6    # BSFG manifold dimensions
D_CRIT   = 26   # critical string dimensions
N_CH     = 9    # channel count
ORD_SO5  = 10   # |SO(5)|


# ---------------------------------------------------------------------------
# Exact rational primitives.
# ---------------------------------------------------------------------------
F_TRZ    = Fraction(1, ORD_SO5)        # = 1/10
SSQ      = Fraction(57, 100)           # [SSq]
PHI_RES  = Fraction(5, 6)              # Φ_res
K_MEX    = Fraction(25, 12)            # K_Mex = Φ_res · |SO(5)| / D_phys
BETA_I   = Fraction(6029, 10000)       # β_i triangular closure
KAPPA_PER_DAY = Fraction(5, 10_000)    # κ = 5e-4 /day

# Cross-check: K_Mex must equal Φ_res · |SO(5)| / D_phys exactly.
assert K_MEX == PHI_RES * ORD_SO5 / Fraction(D_PHYS), \
    "K_Mex consistency broken: 25/12 != (5/6)*10/4"
# Cross-check: F_TRZ must equal 1 / |SO(5)| exactly.
assert F_TRZ == Fraction(1, ORD_SO5), "F_TRZ != 1/|SO(5)|"


# ---------------------------------------------------------------------------
# Exact symbolic primitives (irrational closed form).
# ---------------------------------------------------------------------------
# ρ_SCm = 4·√π · 10⁻³⁷  J/m³   (exact symbolic; float on demand)
RHO_SCM_SYM = 4 * sqrt(pi) * Rational(1, 10) ** 37     # symbolic J/m³
RHO_UA_SYM  = 10 * RHO_SCM_SYM                          # G9: ρ_UA = 10·ρ_SCm

# v_UA = c/3 — the one place c enters.
V_UA_M_PER_S = Fraction(C_LIGHT_M_PER_S, 3)             # exact rational m/s


# ---------------------------------------------------------------------------
# Calibrated (not yet first-principles-derived) — flagged.
# ---------------------------------------------------------------------------
H_SCM_APPROX = Fraction(99, 100)  # ≈ 0.99; CALIBRATED — derivation pending.
H_SCM_IS_DERIVED = False


# ---------------------------------------------------------------------------
# Float accessors (computed, never stored).
# ---------------------------------------------------------------------------
def rho_scm_float() -> float:
    """ρ_SCm in J/m³ as float (computed from symbolic 4√π·10⁻³⁷)."""
    return float(RHO_SCM_SYM)


def rho_ua_float() -> float:
    """ρ_UA in J/m³ as float."""
    return float(RHO_UA_SYM)


def v_ua_float() -> float:
    """v_UA in m/s as float (= c/3)."""
    return float(V_UA_M_PER_S)


# ---------------------------------------------------------------------------
# Self-test: run `python _primitives.py` to verify the anchor.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    print("UQFF locked primitives — anchor self-test")
    print("-" * 56)
    print(f"  ρ_SCm   = 4√π·10⁻³⁷ J/m³   ≈ {rho_scm_float():.6e}")
    print(f"  ρ_UA    = 10·ρ_SCm  J/m³   ≈ {rho_ua_float():.6e}")
    print(f"  F_TRZ   = {F_TRZ}              (= 1/|SO(5)|)")
    print(f"  [SSq]   = {SSQ}             (= 57/100)")
    print(f"  Φ_res   = {PHI_RES}              (= 5/6)")
    print(f"  K_Mex   = {K_MEX}             (= Φ_res·|SO(5)|/D_phys)")
    print(f"  β_i     = {BETA_I}         (= 6029/10000)")
    print(f"  v_UA    = c/3   m/s        ≈ {v_ua_float():.6e}")
    print(f"  κ       = {KAPPA_PER_DAY} /day      (= 5e-4 /day)")
    print(f"  H_SCm   = {H_SCM_APPROX}            (CALIBRATED — not derived)")
    print(f"  D_phys={D_PHYS}, D_BSFG={D_BSFG}, D_crit={D_CRIT}, "
          f"N_ch={N_CH}, |SO(5)|={ORD_SO5}")
    print("-" * 56)
    print("anchor OK.")
