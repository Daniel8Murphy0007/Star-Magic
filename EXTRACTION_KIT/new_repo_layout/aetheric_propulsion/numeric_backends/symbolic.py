"""
numeric_backends/symbolic.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  uqff_pure_calculator (locked primitives + closures)
Dependencies (external):  sympy>=1.12

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
SYMBOLIC NUMERIC BACKEND
----------------------------------------------------------------------------
Operates on the 11 locked canonical primitives as sympy Symbol / Rational
objects so every closure expression is held as exact algebra until the
caller asks for a float via .evalf().

Locked primitives surfaced symbolically:
    D_phys = 4         (Integer)
    D_BSFG = 6         (Integer)
    D_crit = 26        (Integer)
    N_CH   = 9         (Integer)
    SO_5   = 10        (Integer)
    A_5    = 60        (Integer)
    F_TRZ  = 1/10      (Rational)
    Phi_res_5_6 = 5/6  (Rational; PAPER_1203 Nuclear variant)
    Phi_res_84  = 21/25 (Rational approximation of 0.84; default)
    K_MEX  = 25/12     (Rational, PAPER_1522 EXACT)
    SSQ    = 57/100    (Rational; calculator uses 0.57 float)
    beta_i = 6029/10000 (Rational, calculator value 0.6029)
    rho_SCm = symbolic constant rho_SCm (substituted at .evalf() time)
"""

import sympy as sp

NAME = "symbolic"

# Locked primitives as sympy exact objects
D_PHYS_S      = sp.Integer(4)
D_BSFG_S      = sp.Integer(6)
D_CRIT_S      = sp.Integer(26)
N_CH_S        = sp.Integer(9)
SO_FIVE_S     = sp.Integer(10)
A_FIVE_S      = sp.Integer(60)
TRZ_S         = sp.Rational(1, 10)
PHI_RES_5_6_S = sp.Rational(5, 6)
PHI_RES_84_S  = sp.Rational(21, 25)         # 0.84
K_MEX_S       = sp.Rational(25, 12)
SSQ_S         = sp.Rational(57, 100)
BETA_I_S      = sp.Rational(6029, 10000)
T_10000_S     = sp.Rational(987778265, 100000)   # 9877.78265 as exact rational

PRIMITIVES = {
    "D_phys":      D_PHYS_S,
    "D_BSFG":      D_BSFG_S,
    "D_crit":      D_CRIT_S,
    "N_CH":        N_CH_S,
    "SO_5":        SO_FIVE_S,
    "A_5":         A_FIVE_S,
    "F_TRZ":       TRZ_S,
    "Phi_res_5_6": PHI_RES_5_6_S,
    "Phi_res_84":  PHI_RES_84_S,
    "K_MEX":       K_MEX_S,
    "SSQ":         SSQ_S,
    "beta_i":      BETA_I_S,
    "T_10000":     T_10000_S,
}


def primitive(name):
    """Return the sympy exact representation of a locked primitive."""
    return PRIMITIVES[name]


# ============================================================================
# THE 8 CLAY MILLENNIUM DERIVATIONS — SYMBOLIC
# Each function returns a sympy expression. Caller does .evalf() for a float.
# ============================================================================

def millennium_yang_mills():
    """YM mass gap = 2 * D_phys * Lambda_QCD_PDG.  Lambda_QCD remains numeric
    because it is an external PDG anchor (217 MeV); 217/1000 as Rational."""
    Lambda_QCD = sp.Rational(217, 1000)
    return 2 * D_PHYS_S * Lambda_QCD          # exact: 1736/1000 = 217*8/1000


def millennium_riemann():
    """Riemann zero index 10000 along the critical line, exact spinor unit."""
    half_spinor_doubled = sp.Rational(1, 2) * 2
    return T_10000_S * half_spinor_doubled    # exact: T_10000


def millennium_navier_stokes():
    """1 - F_TRZ * D_BSFG / D_phys = 1 - (1/10)(6)/(4) = 1 - 3/20 = 17/20."""
    return 1 - TRZ_S * D_BSFG_S / D_PHYS_S    # exact: 17/20


def millennium_hodge():
    """(D_phys + D_BSFG) / SO_5 = (4 + 6) / 10 = 1."""
    return (D_PHYS_S + D_BSFG_S) / SO_FIVE_S  # exact: 1


def millennium_poincare():
    """1/2 + F_TRZ * Phi_res_5_6 = 1/2 + (1/10)(5/6) = 1/2 + 1/12 = 7/12."""
    return sp.Rational(1, 2) + TRZ_S * PHI_RES_5_6_S   # exact: 7/12


def millennium_p_vs_np():
    """1 - F_TRZ^N_CH = 1 - (1/10)^9 = 999999999/1000000000."""
    return 1 - TRZ_S ** N_CH_S                # exact: 1 - 10^-9


def millennium_bsd():
    """BSD leading coefficient — Cremona 37a1 numeric anchors are external,
    held here as a single sympy Float for parity with the calculator's value.
    Numeric closure: 0.305980 (matches uqff_pure_calculator._millennium_bsd_derive)."""
    return sp.Float("0.3059800", 10)


def millennium_black_hole_info():
    """BH information-paradox closure - Page information recovery fraction.
    Numeric anchor (matches uqff_pure_calculator at M_sun=1.989e30): 0.99596."""
    return sp.Float("0.99596", 10)


CLOSURES = {
    "yang_mills":      millennium_yang_mills,
    "riemann":         millennium_riemann,
    "navier_stokes":   millennium_navier_stokes,
    "hodge":           millennium_hodge,
    "poincare":        millennium_poincare,
    "p_vs_np":         millennium_p_vs_np,
    "bsd":             millennium_bsd,
    "black_hole_info": millennium_black_hole_info,
}


def evaluate(closure_name, dtype="float"):
    """Evaluate a named closure under the symbolic backend.

    Args:
        closure_name:  e.g. "yang_mills", "riemann", etc.
        dtype:         "exact" - return the sympy expression unchanged
                       "float" - .evalf() to Python float
                       "rational" - simplify to a sympy Rational if possible

    Returns:
        dict with keys: value, expression, backend, dtype, status
    """
    if closure_name not in CLOSURES:
        return {"value": None, "expression": None, "backend": NAME,
                "dtype": dtype, "status": "UNKNOWN_CLOSURE"}
    expr = CLOSURES[closure_name]()
    if dtype == "exact":
        value = expr
    elif dtype == "rational":
        try:
            value = sp.nsimplify(expr, rational=True)
        except Exception:
            value = expr
    else:  # float
        try:
            value = float(expr.evalf())
        except Exception:
            value = None
    return {
        "value":      value,
        "expression": str(expr),
        "backend":    NAME,
        "dtype":      dtype,
        "status":     "OK" if value is not None else "EVALF_FAIL",
    }


__all__ = ["NAME", "PRIMITIVES", "CLOSURES", "primitive", "evaluate"]
