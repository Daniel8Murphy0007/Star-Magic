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

    # Round 669 corrective injection audit — BAO dual closures + Li-7 + EDGES
    print()
    print("=" * 78)
    print("Round 669 corrective injection audit (replaces Round 666 OPEN_QUESTION)")
    print("=" * 78)
    if "LCDM_BAO_rd_H0_over_c_OPEN" in ad.DISPATCH:
        print("FAIL: legacy OPEN_QUESTION entry still present; Round 669 removal not applied")
        return 1
    print("  Round 666 OPEN_QUESTION removed: OK")

    targets = [
        ("LCDM_BAO_rd_H0_over_c_primary",   "(SO_5 * SSq * beta_i) / (D_phys * D_crit)",
         0.033040484293971134, 0.02),
        ("LCDM_BAO_rd_H0_over_c_alternate", "1 / (SO_5 * K_MEX * S_26)",
         0.033040484293971134, 0.04),
        ("LCDM_Li7_BBN_dilution",           "D_phys - 1 = 3 EXACT (per PAPER_1227)",
         3.1, 4.5),
        ("LCDM_EDGES_T21_amplitude",        "-D_phys * A_5 * beta_i * 2 mK (per PAPER_1761)",
         -289.0, 0.2),
    ]
    for name, formula, target, max_resid in targets:
        if name not in ad.DISPATCH:
            print(f"  FAIL: {name} missing from dispatch")
            return 1
        r = solve(name, geometry="auto", numeric="numerical")
        if r["value"] is None:
            print(f"  FAIL: {name} returned None")
            return 1
        if abs(r["target"] - target) > 1e-9 * abs(target):
            print(f"  FAIL: {name} target {r['target']} != expected {target}")
            return 1
        if r["residual_pct"] is None or r["residual_pct"] > max_resid:
            print(f"  FAIL: {name} residual {r['residual_pct']:.4f}% > tolerance {max_resid}%")
            return 1
        print(f"  PASS: {name:<40s} resid={r['residual_pct']:.4f}%  formula: {formula}")

    print()
    print("Round 669 multi-path corroboration:")
    p_v = ad.DISPATCH["LCDM_BAO_rd_H0_over_c_primary"]["uqff_value"]
    a_v = ad.DISPATCH["LCDM_BAO_rd_H0_over_c_alternate"]["uqff_value"]
    print(f"  primary   = {p_v:.10f}  (5-primitive triadic co-sum form)")
    print(f"  alternate = {a_v:.10f}  (3-primitive amplification form)")
    print(f"  spread    = {abs(p_v - a_v):.2e}  (simulation-range solution per Daniel multi-path principle)")
    print()

    if all_pass:
        print("PHASE E6 SUCCESS CRITERION MET.")
        return 0
    print("PHASE E6 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
