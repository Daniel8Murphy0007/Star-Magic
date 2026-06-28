"""
test_phase_e5_cm_bio_geo_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

PHASE E5 SUCCESS CRITERION: every CM, bio, geo observable in DISPATCH passes
the standard regression contract (value, residual within doc + slack, owner
geometry match, primary_source flow-through, chain >= 4 lines).
"""
import sys
import assimilation_dispatch as ad
from qcalcgeom_solver import solve


def main():
    print("=" * 78)
    print("PHASE E5 REGRESSION — CM + bio + geo observables via solve()")
    print("=" * 78)
    print(f"DISPATCH cumulative: {len(ad.DISPATCH)} observables")
    print(f"  E1={ad.TOTAL_E1} E2={ad.TOTAL_E2} E3={ad.TOTAL_E3} E4={ad.TOTAL_E4} E5={ad.TOTAL_E5}")
    print()
    targets = ad.observables_by_domain("CM") + ad.observables_by_domain("bio") + ad.observables_by_domain("geo")
    print(f"Observables to test: {len(targets)} (CM + bio + geo)")
    print()
    all_pass = True
    rows = []
    for name in sorted(targets):
        rec = ad.DISPATCH[name]
        r = solve(name, geometry="auto", numeric="numerical")
        ok = True; reasons = []
        if r["value"] is None: ok=False; reasons.append("value None")
        if r["residual_pct"] is not None and r["residual_pct"] > rec["residual_pct"] + 0.5:
            ok=False; reasons.append(f"residual {r['residual_pct']:.4f}% > doc {rec['residual_pct']}% + slack")
        if r["geometry_used"] != rec["owner_geometry"]:
            ok=False; reasons.append(f"geom {r['geometry_used']} != {rec['owner_geometry']}")
        if not r["provenance_chain"] or len(r["provenance_chain"]) < 4:
            ok=False; reasons.append(f"chain {len(r['provenance_chain'])} < 4")
        if r["primary_source"] != rec["primary_source"]:
            ok=False; reasons.append("source mismatch")
        rows.append((name, rec["domain"], r["geometry_used"], r["assimilation_status"], r["residual_pct"], "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}: {reasons}")

    print()
    print("=" * 78)
    print(f"{'observable':<38s} {'dom':<4s} {'owner':<10s} {'status':<8s} {'residual%':>10s} verdict")
    print("=" * 78)
    for n,d,g,s,rp,v in rows:
        rstr = f"{rp:.4f}" if rp is not None else "n/a"
        print(f"{n:<38s} {d:<4s} {g:<10s} {s:<8s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[5] == "PASS")
    print(f"PHASE E5 total: {passing} / {len(rows)} CM+bio+geo observables pass")
    if all_pass:
        print("PHASE E5 SUCCESS CRITERION MET.")
        return 0
    print("PHASE E5 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
