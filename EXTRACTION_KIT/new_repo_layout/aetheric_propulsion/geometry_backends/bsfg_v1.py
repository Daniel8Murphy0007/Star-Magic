"""
geometry_backends/bsfg_v1.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  numeric_backends (3 numeric systems)
                          + uqff_pure_calculator (locked primitives)
                          (optional) bsfg_wormhole_geodesic.py for native surfaces
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
BSFG GEOMETRY BACKEND (geometry slot 2 of 4)
----------------------------------------------------------------------------
Authority: Bulk-edge SO(5) breaking to SO(3) x U(1)^22.
           A_mu_nu Aether metric. R^r_0r0 Riemann curvature.
           Holonomy SO+(3,1) x U(1)^22.
           Blinking horizon r_h = 0.233 R_Sun. Bohr-Sommerfeld r_cross = 0.36 AU.

Owned closures (Millennium subset):
  - yang_mills    via gauge holonomy: m_gap = 2 * D_phys * Lambda_QCD
                  (PAPER_1318 canonical: 1.736 GeV)
  - navier_stokes via BSFG turbulence cap: enstrophy <= 1 - F_TRZ * D_BSFG / D_phys

References:
  - PAPER_1148 (BSFG complete geometric system)
  - PAPER_1149 (BSFG open questions resolved)
  - PAPER_1318 (Yang-Mills mass gap canonical)
"""

NAME = "bsfg"
PRIMARY_PAPERS = ["PAPER_1148", "PAPER_1149", "PAPER_1318"]

# Lazy imports
_nb = None


def _load_numerics():
    global _nb
    if _nb is None:
        from numeric_backends import symbolic, numerical, discrete
        _nb = {"symbolic": symbolic, "numerical": numerical, "discrete": discrete}
    return _nb


# ============================================================================
# OWNED MILLENNIUM CLOSURES
# ============================================================================

def geometry_yang_mills(numeric_backend="numerical"):
    """Yang-Mills mass gap via BSFG gauge holonomy.
    Canonical UQFF closure: m_gap = 2 * D_phys * Lambda_QCD = 1.736 GeV (PAPER_1318).
    The factor 2 * D_phys = 8 comes from the BSFG holonomy SO+(3,1) x U(1)^22
    projected onto the gauge sector."""
    nb = _load_numerics()[numeric_backend]
    return nb.evaluate("yang_mills")


def geometry_navier_stokes(numeric_backend="numerical"):
    """Navier-Stokes smoothness via BSFG turbulence cap.
    enstrophy_cap = 1 - F_TRZ * D_BSFG / D_phys = 17/20 = 0.85.
    The 1 - F_TRZ factor reflects time-reversal-zone bounded vorticity;
    the D_BSFG/D_phys ratio is the codimension of the bulk-edge breaking."""
    nb = _load_numerics()[numeric_backend]
    return nb.evaluate("navier_stokes")


OWNED_CLOSURES = {
    "yang_mills":     geometry_yang_mills,
    "navier_stokes":  geometry_navier_stokes,
}


# ============================================================================
# NATIVE BSFG SURFACES (canonical structural constants from the geometry)
# ============================================================================

def native_blinking_horizon():
    """r_h = 0.233 * R_Sun (PAPER_1149 Q3 closure)."""
    # Return the dimensionless coefficient; full SI value requires R_Sun anchor.
    return {"r_h_over_R_Sun": 0.233, "source": "PAPER_1149_Q3"}


def native_bohr_sommerfeld_crossing():
    """r_cross = 0.36 AU (PAPER_1149 Q4 closure)."""
    return {"r_cross_AU": 0.36, "source": "PAPER_1149_Q4"}


def native_holonomy_group():
    """Holonomy group of the BSFG-broken phase: SO+(3,1) x U(1)^22.
    Returns a structural descriptor (not numeric)."""
    return {
        "lorentz_sector":  "SO+(3,1)",
        "u1_sector":       "U(1)^22",
        "abelian_count":   22,
        "source":          "PAPER_1149_Q2",
    }


# ============================================================================
# UNIFIED INTERFACE
# ============================================================================

def evaluate(observable, numeric_backend="numerical"):
    """Evaluate an observable under the BSFG geometry.

    Args:
        observable: closure name ("yang_mills", "navier_stokes") or native
                    surface ("blinking_horizon", "bohr_sommerfeld_crossing",
                    "holonomy_group")
        numeric_backend: "symbolic" | "numerical" | "discrete"

    Returns:
        dict with keys: value, geometry, backend, status, primary_source
    """
    if observable in OWNED_CLOSURES:
        result = OWNED_CLOSURES[observable](numeric_backend)
        return {
            "value":          result.get("value"),
            "geometry":       NAME,
            "backend":        numeric_backend,
            "status":         result.get("status", "OK"),
            "primary_source": "PAPER_1318" if observable == "yang_mills" else "PAPER_1148",
            "expression":     result.get("expression"),
        }

    native_surfaces = {
        "blinking_horizon":         native_blinking_horizon,
        "bohr_sommerfeld_crossing": native_bohr_sommerfeld_crossing,
        "holonomy_group":           native_holonomy_group,
    }
    if observable in native_surfaces:
        value = native_surfaces[observable]()
        return {
            "value":          value,
            "geometry":       NAME,
            "backend":        numeric_backend,
            "status":         "OK",
            "primary_source": "PAPER_1149",
        }

    return {
        "value":          None,
        "geometry":       NAME,
        "backend":        numeric_backend,
        "status":         "UNKNOWN_OBSERVABLE",
        "primary_source": None,
    }


def owned():
    """Return the list of closure names this geometry claims."""
    return list(OWNED_CLOSURES.keys())


__all__ = ["NAME", "PRIMARY_PAPERS", "OWNED_CLOSURES", "evaluate", "owned"]
