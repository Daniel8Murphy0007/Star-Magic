"""Smoke test for the Aetheric-Propulsion package.

Runs after pip install to verify the bundle is complete and the dispatch routes
through the solver bus correctly. Mirrors Star-Magic's ci_smoke.py contract.
"""
import sys


def main():
    print("=" * 60)
    print(f"Python: {sys.version}")
    print("=" * 60)

    import aetheric_propulsion as ap
    print(f"aetheric_propulsion {ap.__version__}")

    obs_total = len(ap.DISPATCH)
    print(f"DISPATCH observables: {obs_total}")
    if obs_total < 100:
        print(f"ERROR: dispatch regression: only {obs_total} observables (< 100)")
        return 1

    n_domains = len(ap.domains())
    print(f"DISPATCH domains: {n_domains}")
    if n_domains != 10:
        print(f"ERROR: expected 10 domains, got {n_domains}")
        return 1

    # alpha_inverse should return 137.0 via solver bus
    r = ap.calculate_analytic_closures(
        {"qcalcgeom_solve": {"observable": "alpha_inverse"}})
    if r.get("value") != 137.0:
        print(f"ERROR: alpha_inverse expected 137.0, got {r}")
        return 1
    print(f"alpha_inverse: {r['value']}")

    # BAO dual closure
    for obs in ("LCDM_BAO_rd_H0_over_c_primary",
                "LCDM_BAO_rd_H0_over_c_alternate"):
        r = ap.calculate_analytic_closures(
            {"qcalcgeom_solve": {"observable": obs, "decompose": True}})
        v = r.get("value", {})
        resid = v.get("residual_pct", 999)
        if resid > 0.05:
            print(f"ERROR: {obs} residual {resid}% exceeds 0.05% pin")
            return 1
        print(f"{obs}: residual={resid:.4f}%")

    print("\nSMOKE: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
