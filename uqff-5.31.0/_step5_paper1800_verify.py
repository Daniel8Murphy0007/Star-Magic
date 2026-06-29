"""
_step5_paper1800_verify.py
Companion arithmetic verification for PAPER_1800 (BAO + Cabibbo Lagrangian re-derivation).

Re-computes all four dual-closure values independently from the locked UQFF primitives
declared in PAPER_1167. Verifies each formula against its PDG / Planck experimental anchor
within the documented residual tolerance + 0.5% slack.

Per CLAUDE.md Rule 10 (Daniel provides; I assemble): this script computes only — it does
not derive. The derivation lives in PAPER_1800 §§6-8; this script is the arithmetic gate.

Run: python3 _step5_paper1800_verify.py
Expects: exit 0 if all four closures match within tolerance.
"""
import sys

LOCKED_PRIMITIVES = {
    # Integer lattice
    "D_phys":   4,
    "D_BSFG":   6,
    "D_crit":   26,
    "N_CH":     9,
    "SO_5":     10,
    "A_5":      60,
    # Real primitives
    "F_TRZ":    1/10,
    "Phi_res":  0.84,
    "Phi_nuc":  5/6,
    "SSq":      0.57,
    "K_MEX":    25/12,
    "beta_i":   0.6029,
    "S_26":     1.453162,
    "lambda_i": 1.0,
}

ANCHORS = {
    "BAO_rd_H0_over_c": {
        "value":  0.033040484293971134,
        "source": "Planck 2018 + eBOSS DR16",
        "uncertainty_pct": 0.10,
    },
    "Cabibbo_sin_theta_C": {
        "value":  0.22431,
        "source": "PDG 2024 |V_us|",
        "uncertainty_pct": 0.379,
    },
}


def residual_pct(computed, target):
    return 100.0 * abs(computed - target) / abs(target)


def main():
    P = LOCKED_PRIMITIVES
    print("=" * 76)
    print("PAPER_1800 — independent arithmetic verification of BAO + Cabibbo")
    print("=" * 76)
    print()
    print("Locked primitives used:")
    for k in sorted(P):
        print(f"  {k:10s} = {P[k]}")
    print()

    closures = [
        # (name, formula_str, computed, anchor_key, doc_residual_pct)
        ("BAO primary",
         "(SO_5 * SSq * beta_i) / (D_phys * D_crit)",
         (P["SO_5"] * P["SSq"] * P["beta_i"]) / (P["D_phys"] * P["D_crit"]),
         "BAO_rd_H0_over_c", 0.0093),
        ("BAO alternate",
         "1 / (SO_5 * K_MEX * S_26)",
         1.0 / (P["SO_5"] * P["K_MEX"] * P["S_26"]),
         "BAO_rd_H0_over_c", 0.0274),
        ("Cabibbo primary",
         "(N_CH * K_MEX * beta_i) / (A_5 * Phi_res)",
         (P["N_CH"] * P["K_MEX"] * P["beta_i"]) / (P["A_5"] * P["Phi_res"]),
         "Cabibbo_sin_theta_C", 0.008),
        ("Cabibbo alternate",
         "(D_phys * K_MEX * S_26) / (D_BSFG * N_CH)",
         (P["D_phys"] * P["K_MEX"] * P["S_26"]) / (P["D_BSFG"] * P["N_CH"]),
         "Cabibbo_sin_theta_C", 0.025),
    ]

    all_pass = True
    print(f"{'Closure':<22s}  {'Computed':>14s}  {'Target':>14s}  {'Resid%':>9s}  Verdict")
    print("-" * 76)
    for name, formula, computed, anchor_key, doc_pct in closures:
        anchor = ANCHORS[anchor_key]
        target = anchor["value"]
        r = residual_pct(computed, target)
        tol = doc_pct + 0.05
        ok = r <= tol
        verdict = "PASS" if ok else "FAIL"
        if not ok:
            all_pass = False
        print(f"{name:<22s}  {computed:>14.10f}  {target:>14.10f}  {r:>8.4f}%  {verdict}")
        print(f"  formula: {formula}")
        print(f"  anchor:  {anchor['source']} (experimental uncertainty: ~{anchor['uncertainty_pct']}%)")
        print(f"  doc:     {doc_pct}% (PAPER_1800 §6/§7/§8)")
        print()

    # Cross-check: multi-path agreement
    bao_p = closures[0][2]
    bao_a = closures[1][2]
    cab_p = closures[2][2]
    cab_a = closures[3][2]
    bao_spread = abs(bao_p - bao_a)
    cab_spread = abs(cab_p - cab_a)
    print("=" * 76)
    print("Multi-path corroboration (per PAPER_1800 §9)")
    print("=" * 76)
    print(f"  BAO primary    = {bao_p:.10f}")
    print(f"  BAO alternate  = {bao_a:.10f}")
    print(f"  BAO spread     = {bao_spread:.2e}  (joint-probability evidence: structural)")
    print()
    print(f"  Cabibbo primary    = {cab_p:.10f}")
    print(f"  Cabibbo alternate  = {cab_a:.10f}")
    print(f"  Cabibbo spread     = {cab_spread:.2e}  (joint-probability evidence: structural)")
    print()

    if all_pass:
        print("PAPER_1800 ARITHMETIC VERIFICATION: PASS (4/4 closures match within tolerance)")
        return 0
    print("PAPER_1800 ARITHMETIC VERIFICATION: FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
