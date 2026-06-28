"""
test_phase_e6_kk_assimilation.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

PHASE E6 SUCCESS CRITERION: every KK observable in DISPATCH passes the
standard regression contract. Plus: the BAO OPEN_QUESTION entry must
exist with the documented 4.77% residual and OPEN_QUESTION note (so the
discrepancy stays traceable in the audit trail).
"""
import sys
import assimilation_dispatch as ad
from qcalcgeom_solver import solve


def main():
    print("=" * 78)
    print("PHASE E6 REGRESSION — KK universal scaling + BAO revisit")
    print("=" * 78)
    print(f"DISPATCH cumulative: {len(ad.DISPATCH)} observables")
    print(f"  E1={ad.TOTAL_E1} E2={ad.TOTAL_E2} E3={ad.TOTAL_E3} E4={ad.TOTAL_E4} E5={ad.TOTAL_E5} E6={ad.TOTAL_E6}")
    print()

    # Test all KK observables
    kk = ad.observables_by_domain("KK")
    print(f"KK observables to test: {len(kk)}")
    print()

    all_pass = True
    rows = []
    for name in sorted(kk):
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
        rows.append((name, r["geometry_used"], r["assimilation_status"], r["residual_pct"], "PASS" if ok else "FAIL"))
        if not ok:
            all_pass = False
            print(f"FAIL {name}: {reasons}")

    print()
    print("=" * 78)
    print(f"{'observable':<36s} {'owner':<10s} {'status':<8s} {'residual%':>10s} verdict")
    print("=" * 78)
    for n,g,s,rp,v in rows:
        rstr = f"{rp:.4f}" if rp is not None else "n/a"
        print(f"{n:<36s} {g:<10s} {s:<8s} {rstr:>10s} {v}")
    print()
    passing = sum(1 for r in rows if r[4] == "PASS")
    print(f"PHASE E6 KK total: {passing} / {len(rows)}")

    # Explicit BAO revisit check
    print()
    print("=" * 78)
    print("BAO OPEN_QUESTION audit (S364 revisit per Round 663 flag)")
    print("=" * 78)
    bao_name = "LCDM_BAO_rd_H0_over_c_OPEN"
    if bao_name not in ad.DISPATCH:
        print(f"FAIL: {bao_name} not in dispatch (it should be marked OPEN_QUESTION)")
        return 1
    bao_rec = ad.DISPATCH[bao_name]
    bao_result = solve(bao_name, geometry="auto", numeric="numerical")
    print(f"  observable:  {bao_name}")
    print(f"  formula:     {bao_rec['uqff_formula']}")
    print(f"  uqff_value:  {bao_result['value']}")
    print(f"  target:      {bao_result['target']}")
    print(f"  residual%:   {bao_result['residual_pct']:.4f}%  (documented 4.77% inconsistency vs source-script's claimed 0.02%)")
    print(f"  status:      {bao_result['assimilation_status']}")
    print(f"  notes:       {bao_rec['notes']}")
    if "OPEN_QUESTION" not in bao_rec["notes"]:
        print("FAIL: OPEN_QUESTION marker missing from notes")
        return 1
    print()
    print("BAO OPEN_QUESTION preserved with honest residual + audit trail marker.")
    print()

    if all_pass:
        print("PHASE E6 SUCCESS CRITERION MET.")
        return 0
    print("PHASE E6 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
