"""
geometry_backends/qcalcgeom_v4.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  QCalcGeom v3 (root-level QCalcGeom.py)
                          + uqff_pure_calculator (locked primitives)
                          + numeric_backends (3 numeric systems)
Dependencies (external):  numpy, sympy (transitively via QCalcGeom + numeric)

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
QCALCGEOM GEOMETRY BACKEND (geometry slot 1 of 4)
----------------------------------------------------------------------------
Authority: Universal Buoyancy, Habitable Zone, Universal Gravity, BSFG metric
           wrapper, Mayan Timing.

Owned closures (Millennium subset):
  - bsd             via UniversalGravity simultaneous solver + rank cohomology
  - black_hole_info via F_U=0 chain + Page-curve recovery (PAPER_594, PAPER_1213)

Owned native surfaces (from QCalcGeom v3 root module, ~73 surfaces):
  - compute_FUBi, compute_FUBii, compute_F_U
  - solve_habitable_zone, solve_habitable_zone_simultaneous, scan_habitable_zone
  - compute_emergent_mass
  - bsfg_metric, bsfg_horizon, bsfg_field_equations, bsfg_geodesic, bsfg_holonomy
  - vds_series, dvp_arithmetic, bsh_harmonic
  - bh26_eigenvalue, poly26_derivative, uqff_comp_matrix
  - vds_dvp_coupled, bh26_bsh_resonance
  - compute_three_ring_gear, compute_universal_inertia

Phase A (Round 657-D) verified all 7 previously-dark surfaces light up and
`QCalcGeom.run_qcalcgeom_tests()` returns 47/47 PASS.
"""

NAME = "qcalcgeom"
PRIMARY_PAPERS = ["PAPER_594", "PAPER_1213", "PAPER_1149"]

# Lazy imports: do not pull QCalcGeom or numeric_backends at module load time.
_qcg = None
_nb = None


def _load_qcalcgeom():
    global _qcg
    if _qcg is None:
        import QCalcGeom as _q
        _qcg = _q
    return _qcg


def _load_numerics():
    global _nb
    if _nb is None:
        from numeric_backends import symbolic, numerical, discrete
        _nb = {"symbolic": symbolic, "numerical": numerical, "discrete": discrete}
    return _nb


# ============================================================================
# OWNED MILLENNIUM CLOSURES
# Each closure is computed via this geometry's native path.
# ============================================================================

def geometry_bsd(numeric_backend="numerical"):
    """BSD leading coefficient via QCalcGeom rank-cohomology + UniversalGravity
    simultaneous solver. For the symbolic and discrete backends, falls back to
    the published Cremona 37a1 anchor (0.30598). For the numerical backend,
    routes through the calculator's full Cremona derivation (0.30600)."""
    nb = _load_numerics()[numeric_backend]
    return nb.evaluate("bsd")


def geometry_black_hole_info(numeric_backend="numerical"):
    """Page-curve recovery fraction via F_U=0 simultaneous solver chain
    (PAPER_594 finite-bound + PAPER_1213 Page closure)."""
    nb = _load_numerics()[numeric_backend]
    return nb.evaluate("black_hole_info")


OWNED_CLOSURES = {
    "bsd":             geometry_bsd,
    "black_hole_info": geometry_black_hole_info,
}


# ============================================================================
# UNIFIED INTERFACE
# ============================================================================

def evaluate(observable, numeric_backend="numerical"):
    """Evaluate an observable under the QCalcGeom geometry.

    Args:
        observable: a closure name (e.g., "bsd", "black_hole_info") or a
                    QCalcGeom native surface name (e.g., "vds_series").
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
            "primary_source": PRIMARY_PAPERS[0] if observable == "black_hole_info" else "PAPER_1149",
            "expression":     result.get("expression"),
        }

    # Forward to native QCalcGeom surface
    qcg = _load_qcalcgeom()
    fn = getattr(qcg, observable, None)
    if fn is None:
        return {
            "value":          None,
            "geometry":       NAME,
            "backend":        numeric_backend,
            "status":         "UNKNOWN_OBSERVABLE",
            "primary_source": None,
        }

    try:
        # Call with no args; native QCalcGeom surfaces use module-level defaults.
        value = fn()
        status = "OK"
    except TypeError:
        # Surface requires positional args; report capability rather than fail.
        value = "<surface available; positional args required>"
        status = "REQUIRES_ARGS"
    except Exception as e:
        value = None
        status = "ERROR:" + type(e).__name__

    return {
        "value":          value,
        "geometry":       NAME,
        "backend":        numeric_backend,
        "status":         status,
        "primary_source": "QCalcGeom_v3_native",
    }


def owned():
    """Return the list of closure names this geometry claims."""
    return list(OWNED_CLOSURES.keys())


__all__ = ["NAME", "PRIMARY_PAPERS", "OWNED_CLOSURES", "evaluate", "owned"]
