"""
test_phase_e2_sm_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction.

Dependencies (internal):  qcalcgeom_solver + assimilation_dispatch
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PHASE E2 SUCCESS CRITERION
----------------------------------------------------------------------------
For every SM-domain observable in `assimilation_dispatch.DISPATCH` added by
Round 662 (TOTAL_E2 entries):
    1. solve(observable) returns a non-None value
    2. residual_pct is within the documented residual tolerance + 0.5% slack
    3. geometry_used matches the dispatch's owner_geometry
    4. primary_source citation flows through
    5. provenance_chain has at least 4 lines

Exit 0 on full success; non-zero with diff otherwise.
"""

import sys

import assimilation_dispatch as ad
from qcalcgeom_solver import solve


def main():
    print("=" * 78)
    print("PHASE E2 REGRESSION — SM free parameters assimilation via solve()")
    print("=" * 78)
    print()
    print(f"Total observables in DISPATCH: {len(ad.DISPATCH)}")
    print(f"  - TOTAL_E1 (Round 661 SI fundamentals etc.): {ad.TOTAL_E1}")
    print(f"  - TOTAL_E2 (Round 662 SM additions):         {ad.TOTAL_E2}")
    print()

    # Filter to SM observables (E2 additions)
    sm_observables = ad.observables_by_domain("SM")
    print(f"SM-domain observables to test: {len(sm_observables)}")
    print()

    all_pass = True
    rows = []

    for name in sm_observables:
        rec = ad.DISPATCH[name]
        result = solve(name, geometry="auto", numeric="numerical")

        v       = result["value"]
        target  = result["target"]
        rpct    = result["residual_pct"]
        gused   = result["geometry_used"]
        status  = result["assimilation_status"]
        chain   = result["provenance_chain"]
        psrc    = result["primary_source"]

        ok = True
        reasons = []

        if v is None:
            ok = False; reasons.append("value is None")
        if rpct is not None and rpct > rec["residual_pct"] + 0.5:
            ok = False; reasons.append(f"residual {rpct:.4f}% > documented {rec['residual_pct']}% + 0.5% slack")
        if gused != rec["owner_geometry"]:
            ok = False; reasons.append(f"geometry_used={gused} != owner_geometry={rec['owner_geometry']}")
        if not chain or len(chain) < 4:
            ok = False; reasons.append(f"provenance_chain length {len(chain)} < 4")
        if psrc != rec["primary_source"]:
            ok = False; reasons.append(f"primary_source={psrc} != dispatch.{rec['primary_source']}")

        rows.append((name, gused, status, rpct, "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}:")
            for r in reasons:
                print(f"     - {r}")
            print()

    print()
    print("=" * 78)
    print(f"{'observable':<32s} {'owner':<10s} {'status':<10s} {'residual%':>10s} {'verdict'}")
    print("=" * 78)
    for n, g, s, r, v in rows:
        rstr = f"{r:.4f}" if r is not None else "n/a"
        print(f"{n:<32s} {g:<10s} {s:<10s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[4] == "PASS")
    total = len(rows)
    print(f"PHASE E2 total: {passing} / {total} SM observables pass assimilation regression")
    print()
    if all_pass:
        print("PHASE E2 SUCCESS CRITERION MET. SM free parameter dispatch operational.")
        return 0
    print("PHASE E2 FAILURE. See per-observable diagnostics above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
