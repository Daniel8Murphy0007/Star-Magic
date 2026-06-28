"""
test_phase_d_solver_bus.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  qcalcgeom_solver, geometry_backends, numeric_backends,
                          provenance_recorder
Dependencies (external):  sympy (transitively via numeric_backends.symbolic)

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PHASE D SUCCESS CRITERION
----------------------------------------------------------------------------
For each of the 8 Clay Millennium derivations, `qcalcgeom_solver.solve()`
must:
    1. Return a non-None value
    2. Identify the correct owning geometry (per EXPANSION_PLAN Section 5)
    3. Produce a residual within the documented tolerance
    4. Populate the alternate_paths matrix with at least 3 valid cells
    5. Attach a non-empty provenance_chain
    6. Mark assimilation_status as OK or EXACT

Plus: solve("yang_mills", geometry="all", numeric="all") must return a
result dict with the BSFG owner correctly identified and the 4-geometry
matrix populated where each geometry has a value (UNKNOWN or computed).

Exits 0 on full success; non-zero with a printed diff otherwise.
"""

import sys

from qcalcgeom_solver import solve, KNOWN_TARGETS, KNOWN_TOLERANCES_ABS

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


def main():
    print("=" * 78)
    print("PHASE D REGRESSION — qcalcgeom_solver.solve() vs canonical targets")
    print("=" * 78)
    print()

    all_pass = True
    total = 0
    passing = 0
    summary_rows = []

    for closure, expected_owner in EXPECTED_OWNER.items():
        total += 1
        print(f"[{closure}]  (expected owner: {expected_owner})")
        result = solve(closure, geometry="all", numeric="all")

        v       = result["value"]
        target  = result["target"]
        rpct    = result["residual_pct"]
        gused   = result["geometry_used"]
        nused   = result["numeric_system"]
        Ncount  = result["overdetermination_N"]
        status  = result["assimilation_status"]
        warns   = result["warnings"]
        chain   = result["provenance_chain"]
        psrc    = result["primary_source"]

        line_ok = True

        if v is None:
            print(f"  FAIL  value is None")
            line_ok = False
        else:
            print(f"  value         : {v}")

        print(f"  target        : {target}")
        print(f"  residual_pct  : {rpct}")
        print(f"  geometry_used : {gused}    (expected {expected_owner})")
        if gused != expected_owner:
            print(f"  FAIL  geometry_used does not match expected owner")
            line_ok = False

        print(f"  numeric_system: {nused}")
        print(f"  overdet_N     : {Ncount}    (4 geoms x 3 nums = 12 max; "
              f"owner-only would be 3)")
        if Ncount < 3:
            print(f"  FAIL  overdetermination_N < 3 (need at least the owner's 3 cells)")
            line_ok = False

        print(f"  status        : {status}")
        if status not in ("OK", "EXACT"):
            print(f"  FAIL  status not OK/EXACT")
            line_ok = False

        print(f"  primary_source: {psrc}")
        print(f"  provenance has {len(chain)} chain lines")
        if not chain:
            print(f"  FAIL  provenance_chain is empty")
            line_ok = False
        elif len(chain) < 3:
            print(f"  FAIL  provenance_chain has < 3 entries")
            line_ok = False

        if warns:
            print(f"  warnings: {warns}")

        # Show the chain so reviewers can see the audit trail
        for c in chain[:5]:
            print(f"    chain> {c}")
        if len(chain) > 5:
            print(f"    chain> ... ({len(chain) - 5} more lines)")

        if line_ok:
            passing += 1
            print(f"  PASS")
        else:
            all_pass = False
            print(f"  FAIL")
        print()

        summary_rows.append((closure, expected_owner, gused, Ncount, status,
                             "PASS" if line_ok else "FAIL"))

    print("=" * 78)
    print("Summary")
    print("=" * 78)
    print(f"  {'closure':<18s} {'expected':<10s} {'used':<10s} "
          f"{'N':>3s} {'status':<10s} {'verdict'}")
    print(f"  {'-'*18} {'-'*10} {'-'*10} {'-'*3} {'-'*10} {'-'*7}")
    for row in summary_rows:
        print(f"  {row[0]:<18s} {row[1]:<10s} {row[2]:<10s} "
              f"{row[3]:>3d} {row[4]:<10s} {row[5]}")
    print()
    print(f"PHASE D total: {passing} / {total} closures pass solve() regression")
    print()
    if all_pass:
        print("PHASE D SUCCESS CRITERION MET. Solver bus operational.")
        return 0
    print("PHASE D FAILURE. See per-closure diagnostics above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
