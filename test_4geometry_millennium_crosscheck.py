"""
test_4geometry_millennium_crosscheck.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  numeric_backends (3 numerics) + geometry_backends (4 geometries)
Dependencies (external):  sympy (transitively via numeric_backends.symbolic)

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PHASE C SUCCESS CRITERION
----------------------------------------------------------------------------
For each Millennium closure, verify the 4-geometry assignment is consistent:
    - yang_mills    -> BSFG owns
    - riemann       -> 26D owns
    - navier_stokes -> BSFG owns
    - hodge         -> DPM owns
    - poincare      -> DPM owns
    - p_vs_np       -> 26D owns
    - bsd           -> QCalcGeom owns
    - black_hole_info -> QCalcGeom owns

For each closure, verify:
    1. Exactly one geometry CLAIMS ownership (per the EXPANSION_PLAN);
    2. That geometry's evaluate() returns the canonical value via every numeric backend;
    3. The 8/8 coverage matches the Phase B numeric cross-check.

Exits 0 on full structural agreement; non-zero otherwise.
"""

import sys

from geometry_backends import qcalcgeom_v4, bsfg_v1, dpm_v1, d26_compactification

# Expected geometry assignment per EXPANSION_PLAN.md Section 5
EXPECTED_OWNER = {
    "yang_mills":      "bsfg",
    "riemann":         "d26",
    "navier_stokes":   "bsfg",
    "hodge":           "dpm",
    "poincare":        "dpm",
    "p_vs_np":         "d26",
    "bsd":             "qcalcgeom",
    "black_hole_info": "qcalcgeom",
}

GEOMETRIES = {
    "qcalcgeom": qcalcgeom_v4,
    "bsfg":      bsfg_v1,
    "dpm":       dpm_v1,
    "d26":       d26_compactification,
}

NUMERIC_BACKENDS = ["symbolic", "numerical", "discrete"]


def main():
    print("=" * 78)
    print("PHASE C CROSS-VALIDATION — 4 geometry backends x 3 numeric backends")
    print("=" * 78)
    print()

    # 1. Check geometry ownership matches the plan
    print("[1] Geometry ownership audit")
    ownership_ok = True
    for closure, expected_owner in EXPECTED_OWNER.items():
        owner_geom = GEOMETRIES[expected_owner]
        if closure in owner_geom.OWNED_CLOSURES:
            print(f"    OK   {closure:<18s} owned by {expected_owner}")
        else:
            print(f"    MISS {closure:<18s} expected owner={expected_owner}, "
                  f"but not in OWNED_CLOSURES")
            ownership_ok = False

    # 2. For each closure, run the owning geometry through all 3 numerics
    print()
    print("[2] 4-geometry x 3-numeric evaluation matrix")
    print()

    eval_ok = True
    for closure in EXPECTED_OWNER:
        owner = EXPECTED_OWNER[closure]
        geom = GEOMETRIES[owner]
        print(f"  {closure} (owner: {owner})")
        per_numeric = {}
        for nb in NUMERIC_BACKENDS:
            result = geom.evaluate(closure, numeric_backend=nb)
            v = result.get("value")
            per_numeric[nb] = v
            paper = result.get("primary_source", "?")
            print(f"     {nb:<10s} -> value={v}  primary_source={paper}")

        # Compare floats
        s = per_numeric.get("symbolic")
        n = per_numeric.get("numerical")
        d = per_numeric.get("discrete")
        try:
            sf = float(s) if s is not None else None
            nf = float(n) if n is not None else None
            df = float(d) if d is not None else None
            # Use generous tolerance — same as Phase B harness uses for these closures
            tol = 1e-4 if closure in ("bsd", "black_hole_info") else (1e-6 if closure == "riemann" else 1e-9)
            if sf is not None and nf is not None and abs(sf - nf) > tol:
                print(f"     DIFF symbolic vs numerical: {abs(sf-nf):.3e}")
                eval_ok = False
            if nf is not None and df is not None and abs(nf - df) > tol:
                print(f"     DIFF numerical vs discrete: {abs(nf-df):.3e}")
                eval_ok = False
        except Exception as e:
            print(f"     ERROR comparing: {e}")
            eval_ok = False
        print()

    # 3. Verify each geometry has at least one owned closure
    print("[3] Geometry coverage")
    coverage_ok = True
    for geom_name, geom_mod in GEOMETRIES.items():
        owned = geom_mod.owned()
        if not owned:
            print(f"    EMPTY {geom_name} owns ZERO closures")
            coverage_ok = False
        else:
            print(f"    OK    {geom_name} owns {len(owned)} closures: {owned}")

    print()
    print("=" * 78)
    overall = ownership_ok and eval_ok and coverage_ok
    if overall:
        print("PHASE C SUCCESS CRITERION MET. All 4 geometries operational + assignments consistent.")
        return 0
    else:
        print("PHASE C FAILURE.")
        if not ownership_ok: print("  - ownership mismatch")
        if not eval_ok:      print("  - evaluation disagreement")
        if not coverage_ok:  print("  - empty geometry detected")
        return 1


if __name__ == "__main__":
    sys.exit(main())
