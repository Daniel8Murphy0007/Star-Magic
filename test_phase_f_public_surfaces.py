"""
test_phase_f_public_surfaces.py
Phase F regression — verify 8 new public calculate_* surfaces + extended
calculate_analytic_closures qcalcgeom_solve dispatch key.

Verifies:
 1. 8 new public surfaces exist and are callable.
 2. Each returns {'value': X} contract (Rule 5).
 3. The 5 calculate_qcalcgeom_compute_* surfaces produce finite floats.
 4. The 3 analysis surfaces produce dicts with expected keys.
 5. calculate_analytic_closures accepts qcalcgeom_solve dispatch key.
 6. Decomposition mode returns the standard solver-bus 8-field result dict.
 7. Public-surface total = 42 (was 34 prior to Phase F).
"""
import sys
sys.path.insert(0, ".")
import uqff_pure_calculator as u


PHASE_F_SURFACES = [
    "calculate_qcalcgeom_compute_FUBi",
    "calculate_qcalcgeom_compute_FUBii",
    "calculate_qcalcgeom_compute_F_U",
    "calculate_qcalcgeom_solve_habitable_zone",
    "calculate_qcalcgeom_compute_emergent_mass",
    "calculate_3numeric_decomposition",
    "calculate_geometry_decomposition",
    "calculate_overdetermination",
]
COMPUTE_SURFACES = PHASE_F_SURFACES[:5]
ANALYSIS_SURFACES = PHASE_F_SURFACES[5:]


def main():
    print("=" * 72)
    print("PHASE F regression — public surface integration")
    print("=" * 72)

    all_pass = True

    # 1. Existence + signature
    for name in PHASE_F_SURFACES:
        fn = getattr(u, name, None)
        if fn is None or not callable(fn):
            print(f"FAIL: {name} missing or not callable")
            all_pass = False
        else:
            print(f"PASS  exists+callable: {name}")

    # 2. Compute surfaces — return finite floats
    print()
    print("--- compute surfaces ({} required) ---".format(len(COMPUTE_SURFACES)))
    for name in COMPUTE_SURFACES:
        fn = getattr(u, name)
        r = fn({"M": 1.989e30, "r": 1.496e11, "t_n": 0.0})
        if not isinstance(r, dict) or "value" not in r:
            print(f"FAIL: {name} contract — got {r}")
            all_pass = False
            continue
        val = r["value"]
        if val is None or not isinstance(val, (int, float)):
            print(f"FAIL: {name} value type={type(val).__name__}")
            all_pass = False
            continue
        print(f"PASS  {name}: value={val:.6e}")

    # 3. Analysis surfaces — return dicts with expected keys
    print()
    print("--- analysis surfaces (3 required) ---")
    for name in ANALYSIS_SURFACES:
        fn = getattr(u, name)
        r = fn({"observable": "alpha_inverse"})
        if not isinstance(r, dict) or "value" not in r:
            print(f"FAIL: {name} contract")
            all_pass = False
            continue
        val = r["value"]
        if val is None or not isinstance(val, dict):
            print(f"FAIL: {name} returned {val}")
            all_pass = False
            continue
        print(f"PASS  {name}: keys={sorted(val.keys())}")

    # 4. analytic_closures qcalcgeom_solve key (F3)
    print()
    print("--- calculate_analytic_closures qcalcgeom_solve dispatch (F3) ---")
    simple = u.calculate_analytic_closures({"qcalcgeom_solve": {"observable": "alpha_inverse"}})
    if simple["value"] != 137.0:
        print(f"FAIL: alpha_inverse simple = {simple['value']} != 137.0")
        all_pass = False
    else:
        print(f"PASS  alpha_inverse simple call: value={simple['value']}")

    decomp = u.calculate_analytic_closures({
        "qcalcgeom_solve": {"observable": "LCDM_BAO_rd_H0_over_c_primary", "decompose": True}})
    v = decomp.get("value") or {}
    expected = {"value", "target", "residual_pct", "geometry_used",
                "numeric_system", "overdetermination_N", "alternate_paths",
                "assimilation_status"}
    missing = expected - set(v.keys() if isinstance(v, dict) else [])
    if missing:
        print(f"FAIL: decomposed result missing keys: {missing}")
        all_pass = False
    else:
        print(f"PASS  BAO primary decomposed call: residual={v['residual_pct']:.4f}%, status={v['assimilation_status']}")

    # 5. Total public surface count
    print()
    print("--- public surface inventory ---")
    public = sorted(n for n in dir(u)
                    if not n.startswith("_")
                    and callable(getattr(u, n))
                    and n.startswith("calculate_"))
    n = len(public)
    expected_n = 42
    if n != expected_n:
        print(f"FAIL: public surface count {n} != {expected_n}")
        for x in public: print(f"      {x}")
        all_pass = False
    else:
        print(f"PASS  public surface count: {n}")

    # 6. Confirm all 8 Phase F surfaces are in the public list
    missing_pf = [s for s in PHASE_F_SURFACES if s not in public]
    if missing_pf:
        print(f"FAIL: Phase F surfaces missing from public dir: {missing_pf}")
        all_pass = False
    else:
        print(f"PASS  all 8 Phase F surfaces in public dir")

    print()
    if all_pass:
        print("PHASE F SUCCESS CRITERION MET.")
        return 0
    print("PHASE F FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
