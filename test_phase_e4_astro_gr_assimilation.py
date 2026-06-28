"""
test_phase_e4_astro_gr_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction.

Dependencies (internal):  qcalcgeom_solver + assimilation_dispatch
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

PHASE E4 SUCCESS CRITERION: for every astro and GR observable in DISPATCH,
solve() returns value, residual within doc tolerance + 0.5%, owner geometry
matches, primary_source flows through, provenance >= 4 lines.
"""
import sys
import assimilation_dispatch as ad
from qcalcgeom_solver import solve


def main():
    print("=" * 78)
    print("PHASE E4 REGRESSION — astro + GR observables via solve()")
    print("=" * 78)
    print()
    print(f"DISPATCH cumulative: {len(ad.DISPATCH)} observables")
    print(f"  TOTAL_E1={ad.TOTAL_E1}  TOTAL_E2={ad.TOTAL_E2}  TOTAL_E3={ad.TOTAL_E3}  TOTAL_E4={ad.TOTAL_E4}")
    print()

    targets = ad.observables_by_domain("astro") + ad.observables_by_domain("GR")
    print(f"Observables to test: {len(targets)} (astro + GR)")
    print()

    all_pass = True
    rows = []
    for name in sorted(targets):
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
            ok = False; reasons.append("source mismatch")

        rows.append((name, rec["domain"], r["geometry_used"], r["assimilation_status"], r["residual_pct"], "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}: {reasons}")
            print()

    print()
    print("=" * 78)
    print(f"{'observable':<36s} {'dom':<6s} {'owner':<10s} {'status':<8s} {'residual%':>10s} verdict")
    print("=" * 78)
    for n, d, g, s, rp, v in rows:
        rstr = f"{rp:.4f}" if rp is not None else "n/a"
        print(f"{n:<36s} {d:<6s} {g:<10s} {s:<8s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[5] == "PASS")
    print(f"PHASE E4 total: {passing} / {len(rows)} astro+GR observables pass assimilation regression")
    if all_pass:
        print("PHASE E4 SUCCESS CRITERION MET. Astrophysical + GR dispatch operational.")
        return 0
    print("PHASE E4 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
