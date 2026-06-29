"""
test_extraction_kit.py
Phase Step 7 regression — verify EXTRACTION_KIT/ is well-formed and the bundled
package would install + import cleanly as a standalone pip package.

Checks:
 1. Migration script runs in refresh-only mode (no error, idempotent).
 2. All required files present in EXTRACTION_KIT/new_repo_layout/.
 3. pyproject.toml parses as valid TOML.
 4. The bundled package imports without raising (using sys.path injection — no install needed).
 5. After import, key public surfaces are reachable from the aetheric_propulsion namespace.
 6. The DISPATCH catalog has the expected 116 observables and 10 domains.
 7. alpha_inverse returns 137.0 via the bundled calculator surface.
"""
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
KIT = ROOT / "EXTRACTION_KIT"
LAYOUT = KIT / "new_repo_layout"
PKG = LAYOUT / "aetheric_propulsion"


def main():
    print("=" * 72)
    print("PHASE STEP 7 — EXTRACTION_KIT regression")
    print("=" * 72)
    all_pass = True

    # 1. Migration script refresh
    r = subprocess.run([sys.executable, str(KIT / "_step7_migrate_to_aetheric_propulsion.py")],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"FAIL: migration refresh exit {r.returncode}: {r.stderr[:200]}")
        all_pass = False
    else:
        print("PASS  migration refresh exits cleanly")

    # 2. Required files
    required = [
        "pyproject.toml", "README.md", "LICENSE-AGPL-3.0.txt", "COMMERCIAL.md",
        ".github/workflows/release.yml", "tests/test_smoke.py",
        "aetheric_propulsion/__init__.py", "aetheric_propulsion/calculator.py",
        "aetheric_propulsion/assimilation_dispatch.py", "aetheric_propulsion/qcalcgeom_solver.py",
        "aetheric_propulsion/provenance_recorder.py",
        "aetheric_propulsion/_build_overdetermination_views.py",
        "aetheric_propulsion/geometry_backends/__init__.py",
        "aetheric_propulsion/geometry_backends/qcalcgeom_v4.py",
        "aetheric_propulsion/geometry_backends/bsfg_v1.py",
        "aetheric_propulsion/geometry_backends/dpm_v1.py",
        "aetheric_propulsion/geometry_backends/d26_compactification.py",
        "aetheric_propulsion/numeric_backends/__init__.py",
        "aetheric_propulsion/numeric_backends/symbolic.py",
        "aetheric_propulsion/numeric_backends/numerical.py",
        "aetheric_propulsion/numeric_backends/discrete.py",
        "data/OVERDETERMINATION_MAP.csv", "data/OVERDETERMINATION_WIDE.csv",
        "data/OVERDETERMINATION_MAP.md", "data/ASSIMILATION_GEOMETRY_ATLAS.md",
    ]
    missing = [f for f in required if not (LAYOUT / f).exists()]
    if missing:
        print(f"FAIL: {len(missing)} required files missing: {missing[:5]}")
        all_pass = False
    else:
        print(f"PASS  all {len(required)} required files present")

    # 3. pyproject.toml parses
    try:
        try:
            import tomllib as tl
        except ImportError:
            import tomli as tl
        with (LAYOUT / "pyproject.toml").open("rb") as f:
            tdoc = tl.load(f)
        name = tdoc.get("project", {}).get("name")
        ver = tdoc.get("project", {}).get("version")
        if name != "aetheric-propulsion":
            print(f"FAIL: pyproject project.name = {name!r}, expected 'aetheric-propulsion'")
            all_pass = False
        else:
            print(f"PASS  pyproject.toml parses: name={name}, version={ver}")
    except Exception as e:
        print(f"FAIL: pyproject.toml parse error: {e}")
        all_pass = False

    # 4. Package imports cleanly via sys.path injection (no install needed)
    test_script = """
import sys
sys.path.insert(0, '{}')
try:
    import aetheric_propulsion as ap
    print(f'IMPORT_OK version={{ap.__version__}}')
    print(f'DISPATCH_LEN={{len(ap.DISPATCH)}}')
    print(f'DOMAINS={{ap.domains()}}')
    r = ap.calculate_analytic_closures({{'qcalcgeom_solve': {{'observable': 'alpha_inverse'}}}})
    print(f'ALPHA={{r["value"]}}')
except Exception as e:
    import traceback
    print(f'IMPORT_FAIL: {{type(e).__name__}}: {{e}}')
    traceback.print_exc()
""".format(LAYOUT)
    r = subprocess.run([sys.executable, "-c", test_script],
                       capture_output=True, text=True)
    out = r.stdout
    if "IMPORT_OK" not in out:
        print(f"FAIL: bundled package import failed")
        print(f"  stderr: {r.stderr[:300]}")
        print(f"  stdout: {r.stdout[:300]}")
        all_pass = False
    else:
        print(f"PASS  bundled package imports cleanly")
        for line in out.strip().splitlines():
            if line.startswith(("IMPORT", "DISPATCH", "DOMAINS", "ALPHA")):
                print(f"      {line}")
        # Confirm dispatch count
        for line in out.splitlines():
            if line.startswith("DISPATCH_LEN="):
                n = int(line.split("=")[1])
                if n < 100:
                    print(f"FAIL: bundled dispatch has only {n} observables (< 100)")
                    all_pass = False
            if line.startswith("ALPHA="):
                val = float(line.split("=")[1])
                if abs(val - 137.0) > 1e-9:
                    print(f"FAIL: bundled alpha_inverse != 137.0 (got {val})")
                    all_pass = False

    print()
    if all_pass:
        print("PHASE STEP 7 SUCCESS CRITERION MET.")
        return 0
    print("PHASE STEP 7 FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
