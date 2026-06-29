"""
numeric_backends/numerical.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  uqff_pure_calculator (locked primitives + closures)
Dependencies (external):  (Python stdlib only)

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
NUMERICAL NUMERIC BACKEND
----------------------------------------------------------------------------
Operates on the 11 locked canonical primitives as Python `float` (IEEE-754
double precision) and delegates each closure to the existing calculator
function. This is the canonical "what the calculator already computes" path
and serves as the cross-validation reference for the symbolic and discrete
backends.

This backend imports `uqff_pure_calculator` and calls its existing public
or private helpers. It DOES NOT recompute or redefine any closure formula
- it routes through the canonical implementation so any future calculator
update flows through here automatically.
"""

NAME = "numerical"

# Defer the heavy import until first use so the package can be inspected
# even on systems where uqff_pure_calculator has not yet been compiled /
# strip-byte-cleaned. (No file writes, no side-effects from import.)
_calc = None


def _calc_module():
    """Lazy import of the calculator."""
    global _calc
    if _calc is None:
        import uqff_pure_calculator as _u
        _calc = _u
    return _calc


def primitive(name):
    """Return the float-precision locked primitive from the calculator."""
    calc = _calc_module()
    # Surface the canonical names; map common aliases.
    aliases = {
        "D_phys":      "D_PHYS",
        "D_BSFG":      "D_BSFG",
        "D_crit":      "D_CRIT",
        "N_CH":        "N_CH",
        "SO_5":        "SO_FIVE",
        "A_5":         "A_FIVE",
        "F_TRZ":       "TRZ",
        "Phi_res_5_6": "PHI_RES_5_6",
        "K_MEX":       "K_MEX",
        "SSQ":         "SSQ",
        "beta_i":      "BETA_I",
        "T_10000":     "T_10000",
        "rho_SCm":     "RHO_SCM",
    }
    attr = aliases.get(name, name)
    return float(getattr(calc, attr))


# ============================================================================
# THE 8 CLAY MILLENNIUM DERIVATIONS — NUMERICAL
# Each function delegates to the calculator's canonical implementation.
# ============================================================================

def millennium_yang_mills():
    return _calc_module()._millennium_yang_mills_derive()


def millennium_riemann():
    return _calc_module()._millennium_riemann_derive()


def millennium_navier_stokes():
    return _calc_module()._millennium_navier_stokes_derive()


def millennium_hodge():
    return _calc_module()._millennium_hodge_derive()


def millennium_poincare():
    return _calc_module()._millennium_poincare_derive()


def millennium_p_vs_np():
    return _calc_module()._millennium_p_vs_np_derive()


def millennium_bsd():
    return _calc_module()._millennium_bsd_derive()


def millennium_black_hole_info():
    return _calc_module()._millennium_black_hole_info_derive()


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
    """Evaluate a named closure under the numerical backend.

    The dtype argument is accepted for interface uniformity with the other
    backends; this backend always returns Python float because that is the
    type of the calculator's canonical output.

    Returns:
        dict with keys: value, expression, backend, dtype, status
    """
    if closure_name not in CLOSURES:
        return {"value": None, "expression": None, "backend": NAME,
                "dtype": "float", "status": "UNKNOWN_CLOSURE"}
    try:
        value = float(CLOSURES[closure_name]())
        status = "OK"
    except Exception as e:
        value = None
        status = "ERROR:" + type(e).__name__
    return {
        "value":      value,
        "expression": "calculator._millennium_" + closure_name + "_derive()",
        "backend":    NAME,
        "dtype":      "float",
        "status":     status,
    }


__all__ = ["NAME", "CLOSURES", "primitive", "evaluate"]
