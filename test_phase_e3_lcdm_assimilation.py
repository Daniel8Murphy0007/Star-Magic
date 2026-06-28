"""
test_phase_e3_lcdm_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction.

Dependencies (internal):  qcalcgeom_solver + assimilation_dispatch
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

PHASE E3 SUCCESS CRITERION: for every LCDM-domain observable in DISPATCH,
solve(observable) must return value, residual within doc tolerance + 0.5%
slack, owner geometry match, primary_source flow-through, provenance >= 4 lines.
"""
import sys
import assimilation_dispatch as ad
from qcalcgeom_solver import solve


def main():
    print("=" * 78)
    print("PHASE E3 REGRESSION — LambdaCDM cosmology via solve()")
    print("=" * 78)
    print()
    print(f"Total observables in DISPATCH: {len(ad.DISPATCH)}")
    print(f"  - TOTAL_E1 (Round 661): {ad.TOTAL_E1}")
    print(f"  - TOTAL_E2 (Round 662): {ad.TOTAL_E2}")
    print(f"  - TOTAL_E3 (Round 663): {ad.TOTAL_E3}")
    print()

    lcdm = ad.observables_by_domain("LCDM")
    print(f"LCDM observables to test: {len(lcdm)}")
    print()

    all_pass = True
    rows = []
    for name in lcdm:
        rec = ad.DISPATCH[name]
        r = solve(name, geometry="auto", numeric="numerical")
        ok = True
        reasons = []

        if r["value"] is None:
            ok = False; reasons.append("value None")
        if r["residual_pct"] is not None and r["residual_pct"] > rec["residual_pct"] + 0.5:
            ok = False; reasons.append(f"residual {r['residual_pct']:.4f}% > doc {rec['residual_pct']}% + slack")
        if r["geometry_used"] != rec["owner_geometry"]:
            ok = False; reasons.append(f"geom {r['geometry_used']} != {rec['owner_geometry']}")
        if not r["provenance_chain"] or len(r["provenance_chain"]) < 4:
            ok = False; reasons.append(f"chain {len(r['provenance_chain'])} < 4")
        if r["primary_source"] != rec["primary_source"]:
            ok = False; reasons.append(f"source mismatch")

        rows.append((name, r["geometry_used"], r["assimilation_status"], r["residual_pct"], "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}: {reasons}")
            print()

    print()
    print("=" * 78)
    print(f"{'observable':<34s} {'owner':<8s} {'status':<10s} {'residual%':>10s} {'verdict'}")
    print("=" * 78)
    for n, g, s, rp, v in rows:
        rstr = f"{rp:.4f}" if rp is not None else "n/a"
        print(f"{n:<34s} {g:<8s} {s:<10s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[4] == "PASS")
    print(f"PHASE E3 total: {passing} / {len(rows)} LCDM observables pass assimilation regression")
    if all_pass:
        print("PHASE E3 SUCCESS CRITERION MET. LambdaCDM cosmology dispatch operational.")
        return 0
    print("PHASE E3 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
