"""
gen_muge_all.py — Orchestrator: generates all 17 MUGE .h files + UQFFSource10.h
v4.84 UQFF Integration Hub

This file imports all individual generator modules and runs them in sequence,
producing every .h file. Run from the Star-Magic repository root.

Run:  python gen_muge_all.py
Output: 17 .h header files
"""
import os
import sys
import time

# Ensure the script directory is on the path
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

# ---- Import all generators ----
from gen_muge_sgr0501    import get_content as gc_sgr0501,    OUTPUT_FILE as of_sgr0501
from gen_muge_sgr1745    import get_content as gc_sgr1745,    OUTPUT_FILE as of_sgr1745
from gen_muge_sgra       import get_content as gc_sgra,       OUTPUT_FILE as of_sgra
from gen_muge_tapestry   import get_content as gc_tapestry,   OUTPUT_FILE as of_tapestry
from gen_muge_wd2        import get_content as gc_wd2,        OUTPUT_FILE as of_wd2
from gen_muge_pillars    import get_content as gc_pillars,    OUTPUT_FILE as of_pillars
from gen_muge_rings      import get_content as gc_rings,      OUTPUT_FILE as of_rings
from gen_muge_assessment import get_content as gc_assessment, OUTPUT_FILE as of_assessment
from gen_muge_ngc2525    import get_content as gc_ngc2525,    OUTPUT_FILE as of_ngc2525
from gen_muge_ngc3603    import get_content as gc_ngc3603,    OUTPUT_FILE as of_ngc3603
from gen_muge_bubble     import get_content as gc_bubble,     OUTPUT_FILE as of_bubble
from gen_muge_antennae   import get_content as gc_antennae,   OUTPUT_FILE as of_antennae
from gen_muge_horsehead  import get_content as gc_horsehead,  OUTPUT_FILE as of_horsehead
from gen_muge_ngc1275    import get_content as gc_ngc1275,    OUTPUT_FILE as of_ngc1275
from gen_muge_hudf       import get_content as gc_hudf,       OUTPUT_FILE as of_hudf
from gen_muge_ngc1792    import get_content as gc_ngc1792,    OUTPUT_FILE as of_ngc1792
from gen_source10        import get_content as gc_source10,   OUTPUT_FILE as of_source10


GENERATORS = [
    # (label, get_content_fn, output_path)
    ("Module 01 — MagnetarSGR0501_4516",   gc_sgr0501,    of_sgr0501),
    ("Module 02 — MagnetarSGR1745_2900",   gc_sgr1745,    of_sgr1745),
    ("Module 03 — SMBHSgrAStar",           gc_sgra,       of_sgra),
    ("Module 04 — StarbirthTapestry",      gc_tapestry,   of_tapestry),
    ("Module 05 — Westerlund2",            gc_wd2,        of_wd2),
    ("Module 06 — PillarsOfCreation",      gc_pillars,    of_pillars),
    ("Module 07 — RingsOfRelativity",      gc_rings,      of_rings),
    ("Module 08 — UQFFLearningAssessment", gc_assessment, of_assessment),
    ("Module 09 — GalaxyNGC2525",          gc_ngc2525,    of_ngc2525),
    ("Module 10 — NGC3603",                gc_ngc3603,    of_ngc3603),
    ("Module 11 — BubbleNebula",           gc_bubble,     of_bubble),
    ("Module 12 — AntennaeGalaxies",       gc_antennae,   of_antennae),
    ("Module 13 — HorseheadNebula",        gc_horsehead,  of_horsehead),
    ("Module 14 — NGC1275",                gc_ngc1275,    of_ngc1275),
    ("Module 15 — HUDFGalaxies",           gc_hudf,       of_hudf),
    ("Module 16 — GalaxyNGC1792",          gc_ngc1792,    of_ngc1792),
    ("Integration Hub — UQFFSource10",     gc_source10,   of_source10),
]


def run_all(verbose: bool = True) -> dict:
    """
    Generate all .h files. Returns summary dict.
    """
    results = {}
    total_lines = 0
    t_start = time.time()

    print("=" * 66)
    print("  gen_muge_all.py — UQFF v4.84 Module Generator")
    print("=" * 66)
    print(f"  Output directory: {SCRIPT_DIR}\n")

    for label, get_content_fn, output_path in GENERATORS:
        try:
            content = get_content_fn()
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(content)
            n_lines  = len(content.splitlines())
            n_bytes  = len(content.encode("utf-8"))
            total_lines += n_lines
            status = "OK"
            if verbose:
                fname = os.path.basename(output_path)
                print(f"  [OK]  {fname:<38}  {n_lines:>5} lines  {n_bytes:>8} bytes"
                      f"  ({label})")
            results[output_path] = {"lines": n_lines, "bytes": n_bytes, "status": status}
        except Exception as exc:
            status = f"FAIL: {exc}"
            print(f"  [FAIL] {label} — {exc}")
            results[output_path] = {"lines": 0, "bytes": 0, "status": status}

    elapsed = time.time() - t_start
    ok_count = sum(1 for v in results.values() if v["status"] == "OK")

    print()
    print("=" * 66)
    print(f"  Completed: {ok_count}/{len(GENERATORS)} files generated")
    print(f"  Total lines: {total_lines:,}")
    print(f"  Elapsed: {elapsed:.2f} s")
    print("=" * 66)

    if ok_count < len(GENERATORS):
        print("\n  [WARNING] Some generators failed. Check output above.")
        sys.exit(1)
    else:
        print("\n  All modules generated successfully — ready for v4.84 commit.")

    return results


def print_module_table():
    """Print a compact summary table of all 17 modules."""
    header = f"{'#':<4} {'Output File':<38} {'Description'}"
    print(header)
    print("-" * 80)
    for i, (label, _, output_path) in enumerate(GENERATORS, 1):
        fname = os.path.basename(output_path)
        print(f"{i:<4} {fname:<38} {label}")


if __name__ == "__main__":
    if "--list" in sys.argv:
        print_module_table()
    else:
        run_all(verbose=("--quiet" not in sys.argv))
