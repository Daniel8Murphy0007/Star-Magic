"""
qcalcgeom_solver.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  geometry_backends (4 geometries)
                          + numeric_backends  (3 numerics)
                          + provenance_recorder
                          (optional) uqff_pure_calculator for target lookup
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
QCALCGEOM v4.0 SOLVER BUS
----------------------------------------------------------------------------
Single public entry point: `solve(observable, geometry, numeric, ...)`.
Implements the 4-geometry x 3-numeric dispatch matrix specified in
EXPANSION_PLAN.md Section 7.

This module is SEPARATE from the legacy `QCalcGeom.py` at repo root. The
legacy module provides the 73 v3 native surfaces and the 47/47 self-test
suite; this module provides the unified `solve()` entry point that wraps
the geometry_backends + numeric_backends architecture built in Phases B+C.

Phase F will integrate `solve()` into `calculate_analytic_closures(dataset)`
so external callers reach it through the existing public API surface.

----------------------------------------------------------------------------
RETURN SHAPE (per EXPANSION_PLAN.md Section 7.2)
----------------------------------------------------------------------------
{
    "observable":          str,
    "value":               float | None,
    "target":              float | None,
    "residual_pct":        float | None,
    "geometry_used":       str,
    "numeric_system":      str,
    "alternate_paths":     dict | None,
    "overdetermination_N": int,
    "provenance_chain":    list[str],
    "primary_source":      str | None,
    "assimilation_status": str,
    "warnings":            list[str],
}
"""

from geometry_backends import qcalcgeom_v4, bsfg_v1, dpm_v1, d26_compactification
from numeric_backends    import symbolic, numerical, discrete
import provenance_recorder

# Phase E dispatch table - lazy-loaded so the solver works even if the
# dispatch file is absent (e.g. during Phase D regression).
try:
    import assimilation_dispatch as _adispatch
except ImportError:
    _adispatch = None

GEOMETRIES = {
    "qcalcgeom": qcalcgeom_v4,
    "bsfg":      bsfg_v1,
    "dpm":       dpm_v1,
    "d26":       d26_compactification,
}

NUMERIC_BACKENDS = {
    "symbolic":  symbolic,
    "numerical": numerical,
    "discrete":  discrete,
}

# Known Millennium target values (the "target" column of master_closures.csv
# for these closures). External targets resolve through this table; future
# Phase E work will load the full master_closures.csv into a richer lookup.
KNOWN_TARGETS = {
    "yang_mills":      1.736,
    "riemann":         9877.78265,
    "navier_stokes":   0.85,
    "hodge":           1.0,
    "poincare":        7.0 / 12.0,
    "p_vs_np":         1.0 - 1e-9,
    "bsd":             0.30598,
    "black_hole_info": 0.99596,
}

# Per-closure tolerance for the assimilation_status classifier.
# Matches the Phase B/C cross-validation harness tolerances.
KNOWN_TOLERANCES_ABS = {
    "yang_mills":      1e-9,
    "riemann":         1e-6,
    "navier_stokes":   1e-12,
    "hodge":           1e-15,
    "poincare":        1e-15,
    "p_vs_np":         1e-15,
    "bsd":             1e-4,
    "black_hole_info": 1e-5,
}


# ============================================================================
# Internal helpers
# ============================================================================

def _find_owner(observable):
    """Return the geometry name that claims ownership, or None."""
    for gname, gmod in GEOMETRIES.items():
        if observable in getattr(gmod, "OWNED_CLOSURES", {}):
            return gname
    return None


def _resolve_geometries(geometry, observable):
    """Resolve the geometry argument into a list of geometry names."""
    if geometry == "all":
        return list(GEOMETRIES.keys())
    if geometry == "auto":
        owner = _find_owner(observable)
        return [owner] if owner else list(GEOMETRIES.keys())
    if geometry in GEOMETRIES:
        return [geometry]
    return []


def _resolve_numerics(numeric):
    """Resolve the numeric argument into a list of backend names."""
    if numeric == "all":
        return list(NUMERIC_BACKENDS.keys())
    if numeric in NUMERIC_BACKENDS:
        return [numeric]
    return []


def _classify_status(residual_pct, tolerance_pct):
    """Classify the assimilation status given a residual."""
    if residual_pct is None:
        return "UNKNOWN_TARGET"
    if residual_pct == 0.0:
        return "EXACT"
    if residual_pct <= tolerance_pct:
        return "OK"
    return "TENSION"


def _check_disagreements(alternate_paths, tol_abs=1e-6):
    """Scan the alternate_paths matrix; produce a warning string per cell pair
    that disagrees beyond the absolute tolerance."""
    warnings = []
    cells = []
    for g, by_n in alternate_paths.items():
        for n, cell in by_n.items():
            v = cell.get("value")
            if v is not None:
                try:
                    cells.append((g, n, float(v)))
                except (TypeError, ValueError):
                    pass
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            g1, n1, v1 = cells[i]
            g2, n2, v2 = cells[j]
            if abs(v1 - v2) > tol_abs and abs(v1) > 0:
                rel = abs(v1 - v2) / abs(v1)
                if rel > 1e-3:
                    warnings.append(
                        f"{g1}/{n1}={v1:.6g} vs {g2}/{n2}={v2:.6g} "
                        f"(abs diff {abs(v1 - v2):.3e})"
                    )
    return warnings


# ============================================================================
# PUBLIC ENTRY: solve()
# ============================================================================

def solve(observable, geometry="auto", numeric="all",
          record_provenance=True, tolerance_pct=1.0, decompose=True):
    """Single public entry point for UQFF Assimilation Geometry.

    Args:
        observable:        closure name (e.g. "alpha", "yang_mills", "omega_b_h2")
        geometry:          "auto" | "all" | "qcalcgeom" | "bsfg" | "dpm" | "d26"
        numeric:           "all" | "symbolic" | "numerical" | "discrete"
        record_provenance: if True, attach a provenance_chain to the result
        tolerance_pct:     classifier threshold for the assimilation_status field
        decompose:         if True, attach the full alternate_paths matrix

    Returns:
        Assimilation result dict per EXPANSION_PLAN.md Section 7.2.
    """
    # Phase E: if the observable is registered in assimilation_dispatch and
    # the geometry backends do not own it natively, route through the dispatch
    # so the canonical UQFF formula + target + paper citation flow through.
    dispatch_rec = _adispatch.lookup(observable) if _adispatch is not None else None
    native_owner = _find_owner(observable)

    if dispatch_rec is not None and native_owner is None:
        return _solve_via_dispatch(observable, dispatch_rec, geometry, numeric,
                                   record_provenance, tolerance_pct, decompose)

    geoms = _resolve_geometries(geometry, observable)
    nums  = _resolve_numerics(numeric)

    # Build the alternate_paths 4x3 (or smaller) matrix
    alternate_paths = {}
    for g in geoms:
        gmod = GEOMETRIES[g]
        alternate_paths[g] = {}
        for n in nums:
            result = gmod.evaluate(observable, numeric_backend=n)
            alternate_paths[g][n] = {
                "value":          result.get("value"),
                "status":         result.get("status"),
                "primary_source": result.get("primary_source"),
                "expression":     result.get("expression"),
            }

    # Pick primary: prefer the canonical owner's numerical backend
    owner = _find_owner(observable)
    primary_geom = owner if owner in alternate_paths else (geoms[0] if geoms else None)
    primary_num  = "numerical" if "numerical" in nums else (nums[0] if nums else None)

    primary_value = None
    primary_source = None
    primary_expression = None
    if primary_geom and primary_num:
        cell = alternate_paths.get(primary_geom, {}).get(primary_num, {})
        primary_value      = cell.get("value")
        primary_source     = cell.get("primary_source")
        primary_expression = cell.get("expression")

    # Target lookup + residual
    target = KNOWN_TARGETS.get(observable)
    residual_pct = None
    if primary_value is not None and target is not None:
        try:
            pv = float(primary_value)
            tv = float(target)
            if tv != 0:
                residual_pct = abs(pv - tv) / abs(tv) * 100.0
            else:
                residual_pct = abs(pv - tv) * 100.0
        except (TypeError, ValueError):
            pass

    # overdetermination_N = number of cells with valid values
    N = sum(
        1
        for g in alternate_paths
        for n in alternate_paths[g]
        if alternate_paths[g][n].get("value") is not None
    )

    # Provenance chain
    chain = []
    if record_provenance:
        chain = provenance_recorder.build_chain(
            observable,
            primary={
                "value":          primary_value,
                "geometry":       primary_geom,
                "numeric":        primary_num,
                "primary_source": primary_source,
            },
            alternate_paths=alternate_paths,
        )

    # Status + warnings
    closure_tol_pct = (
        KNOWN_TOLERANCES_ABS.get(observable, 1e-6) * 100
        if observable in KNOWN_TARGETS
        else tolerance_pct
    )
    status = _classify_status(residual_pct, closure_tol_pct)
    warnings = _check_disagreements(alternate_paths)

    return {
        "observable":          observable,
        "value":               primary_value,
        "target":              target,
        "residual_pct":        residual_pct,
        "geometry_used":       primary_geom,
        "numeric_system":      primary_num,
        "alternate_paths":     alternate_paths if decompose else None,
        "overdetermination_N": N,
        "provenance_chain":    chain,
        "primary_source":      primary_source,
        "assimilation_status": status,
        "warnings":            warnings,
        "_primary_expression": primary_expression,
    }




# ============================================================================
# Phase E dispatch fallback - for observables wired via assimilation_dispatch
# but not yet present in any geometry backend's OWNED_CLOSURES.
# ============================================================================

def _solve_via_dispatch(observable, rec, geometry, numeric,
                       record_provenance, tolerance_pct, decompose):
    """Solve an observable whose canonical formula lives in
    assimilation_dispatch.DISPATCH (Phase E wired closures).
    Builds a synthetic alternate_paths entry under the dispatch's
    owner_geometry slot for the requested numeric backend(s)."""
    owner_geom = rec["owner_geometry"]
    target = rec["target"]
    uqff_value = rec["uqff_value"]
    formula = rec["uqff_formula"]
    paper = rec["primary_source"]

    nums = _resolve_numerics(numeric)

    alt_paths = {owner_geom: {}}
    for n in nums:
        alt_paths[owner_geom][n] = {
            "value":          uqff_value,
            "status":         "OK",
            "primary_source": paper,
            "expression":     formula,
        }

    primary_num = "numerical" if "numerical" in nums else (nums[0] if nums else "numerical")
    primary_value = uqff_value

    residual_pct = None
    if target is not None:
        try:
            tv = float(target)
            if tv != 0:
                residual_pct = abs(uqff_value - tv) / abs(tv) * 100.0
            else:
                residual_pct = abs(uqff_value - tv) * 100.0
        except (TypeError, ValueError):
            pass

    N = sum(1 for g in alt_paths for n in alt_paths[g]
            if alt_paths[g][n].get("value") is not None)

    chain = []
    if record_provenance:
        chain = [
            f"closure   : {observable} (domain: {rec['domain']})",
            f"geometry  : {owner_geom} (canonical owner; primary source: {paper})",
            f"numeric   : {primary_num} backend",
            f"formula   : {formula}",
            f"target    : {target}",
            f"value     : {uqff_value}",
            f"residual  : {residual_pct:.4f}%" if residual_pct is not None else "residual  : (no target)",
            f"session   : {rec.get('session_script') or '(no session script)'}",
        ]
        if rec.get("notes"):
            chain.append(f"notes     : {rec['notes']}")

    closure_tol_pct = 1.0  # default 1% tolerance for E1 closures
    status = _classify_status(residual_pct, closure_tol_pct)

    return {
        "observable":          observable,
        "value":               primary_value,
        "target":              target,
        "residual_pct":        residual_pct,
        "geometry_used":       owner_geom,
        "numeric_system":      primary_num,
        "alternate_paths":     alt_paths if decompose else None,
        "overdetermination_N": N,
        "provenance_chain":    chain,
        "primary_source":      paper,
        "assimilation_status": status,
        "warnings":            [],
        "_primary_expression": formula,
    }


__all__ = ["solve", "GEOMETRIES", "NUMERIC_BACKENDS",
           "KNOWN_TARGETS", "KNOWN_TOLERANCES_ABS",
           "_solve_via_dispatch"]
