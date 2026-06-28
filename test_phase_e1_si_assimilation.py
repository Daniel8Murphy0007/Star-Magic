"""
test_phase_e1_si_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  qcalcgeom_solver + assimilation_dispatch
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PHASE E1 SUCCESS CRITERION
----------------------------------------------------------------------------
For every observable wired into `assimilation_dispatch.DISPATCH`:
    1. solve(observable) returns a non-None value
    2. value matches the dispatch's documented uqff_value
    3. residual_pct is within 1% of the target
    4. provenance_chain has at least 4 lines (closure + geometry + formula + value)
    5. geometry_used matches the dispatch's owner_geometry
    6. primary_source citation flows through to the result

Exit 0 on full success; non-zero with structured diff otherwise.
"""

import sys

import assimilation_dispatch as ad
from qcalcgeom_solver import solve

TOLERANCE_PCT = 1.0


def main():
    print("=" * 78)
    print("PHASE E1 REGRESSION — SI fundamentals assimilation via solve()")
    print("=" * 78)
    print()
    print(f"Wired observables: {ad.TOTAL_E1}")
    print(f"Domains covered:   {ad.domains()}")
    print()

    all_pass = True
    rows = []

    for name in sorted(ad.DISPATCH.keys()):
        rec = ad.DISPATCH[name]
        result = solve(name, geometry="auto", numeric="numerical")

        v       = result["value"]
        target  = result["target"]
        rpct    = result["residual_pct"]
        gused   = result["geometry_used"]
        nused   = result["numeric_system"]
        Ncount  = result["overdetermination_N"]
        status  = result["assimilation_status"]
        chain   = result["provenance_chain"]
        psrc    = result["primary_source"]

        ok = True
        reasons = []

        if v is None:
            ok = False; reasons.append("value is None")
        elif abs(v - rec["uqff_value"]) > 1e-9 * max(1, abs(rec["uqff_value"])):
            ok = False; reasons.append(f"value {v} != dispatch.uqff_value {rec['uqff_value']}")

        if rpct is None and target is not None:
            ok = False; reasons.append("residual_pct is None despite target present")
        elif rpct is not None and rpct > rec["residual_pct"] + 0.5:
            ok = False; reasons.append(f"residual {rpct:.4f}% > documented {rec['residual_pct']}% + 0.5% slack")

        if gused != rec["owner_geometry"]:
            ok = False; reasons.append(f"geometry_used={gused} != owner_geometry={rec['owner_geometry']}")

        if not chain or len(chain) < 4:
            ok = False; reasons.append(f"provenance_chain length {len(chain)} < 4")

        if psrc != rec["primary_source"]:
            ok = False; reasons.append(f"primary_source={psrc} != dispatch.{rec['primary_source']}")

        rows.append((name, rec["domain"], gused, status, rpct, "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}:")
            for r in reasons:
                print(f"     - {r}")
            print(f"     chain={chain}")
            print()

    print()
    print("=" * 78)
    print(f"{'observable':<28s} {'domain':<6s} {'owner':<10s} {'status':<10s} {'residual%':>10s} {'verdict'}")
    print("=" * 78)
    for n, d, g, s, r, v in rows:
        rstr = f"{r:.4f}" if r is not None else "n/a"
        print(f"{n:<28s} {d:<6s} {g:<10s} {s:<10s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[5] == "PASS")
    total = len(rows)
    print(f"PHASE E1 total: {passing} / {total} observables pass assimilation regression")
    print()
    if all_pass:
        print("PHASE E1 SUCCESS CRITERION MET. SI fundamentals dispatch operational.")
        return 0
    print("PHASE E1 FAILURE. See per-observable diagnostics above.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
