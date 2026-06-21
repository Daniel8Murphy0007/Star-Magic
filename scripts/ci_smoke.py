#!/usr/bin/env python3
"""CI smoke test for UQFF Star-Magic.

Runs after the fidelity gate to verify the calculator's public surface
matches a soft contract. Exits non-zero on any HARD regression (primitive
count drift, missing imports). Exits zero with warnings on soft drift
(closure-count growth, surface count growth).

Usage:  python scripts/ci_smoke.py
"""
import os
import sys
import traceback

# Make the calculator importable regardless of where this script lives.
# CI runs this from repo root; calculator file is in repo root.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
for _p in (os.getcwd(), _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def main() -> int:
    print("=" * 60)
    print(f"Python: {sys.version}")
    print(f"Platform: {sys.platform}")
    print("=" * 60)

    try:
        import uqff_pure_calculator as u
    except SystemExit as e:
        print(f"ERROR: import raised SystemExit({e.code})")
        return 1
    except Exception:
        print("ERROR: import of uqff_pure_calculator failed:")
        traceback.print_exc()
        return 1

    publics = sorted(n for n in dir(u) if n.startswith("calculate_"))
    print(f"\npublic calculate_* surfaces: {len(publics)}")
    for p in publics:
        print(f"  {p}")

    # SOFT lower bound on surface count
    if len(publics) < 30:
        print(f"\nERROR: public surface regression: only {len(publics)} surfaces (< 30)")
        return 1

    try:
        report = u.calculate_status_report({})
    except Exception:
        print("\nERROR: calculate_status_report({}) raised:")
        traceback.print_exc()
        return 1

    s = report.get("value", {}).get("summary", {})
    if not s:
        print(f"\nERROR: status_report returned no summary: {report!r}")
        return 1

    print(f"\nclosures: {s.get('total_closures')}")
    print(f"truly_independent_primitives: {s.get('truly_independent_primitives')}")
    print(f"derivative_primitives: {s.get('derivative_primitives')}")

    # SOFT lower bound on closures
    if s.get("total_closures", 0) < 700:
        print(f"\nERROR: closure regression: only {s.get('total_closures')} (< 700)")
        return 1

    # HARD invariant: the LANDMARK reduction is locked
    if s.get("truly_independent_primitives") != 9:
        print(f"\nERROR: truly_independent_primitives != 9 (got {s.get('truly_independent_primitives')})")
        return 1

    if s.get("derivative_primitives", 0) < 2:
        print(f"\nERROR: derivative_primitives < 2 (got {s.get('derivative_primitives')})")
        return 1

    # Holmlid KER — sanity bounded
    try:
        ker = u.calculate_lenr({})["value"]["ker_chain"]["E_phonon_eV"]
    except Exception:
        print("\nERROR: calculate_lenr({}) ker_chain access failed:")
        traceback.print_exc()
        return 1
    print(f"\nHolmlid E_phonon: {ker} eV")
    if not (0.001 < ker < 0.020):
        print(f"ERROR: E_phonon outside sane range: {ker}")
        return 1

    # U_i for the Sun at t=0 should be 2.75e-7 per PAPER_646
    try:
        ui = u.calculate_universal_inertial_operator({})["value"]["U_i_dimensionless"]
    except Exception:
        print("\nERROR: calculate_universal_inertial_operator({}) failed:")
        traceback.print_exc()
        return 1
    print(f"U_i: {ui}")
    if abs(ui - 2.75e-7) > 1e-8:
        print(f"ERROR: U_i drift: {ui} (expected 2.75e-7)")
        return 1

    print("\nSMOKE: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
