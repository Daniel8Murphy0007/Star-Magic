"""
gen_fubiicalc.py  —  Orchestrator for F_U_Bi_i_QCalc.cpp generation
Run:  python gen_fubiicalc.py
Output: F_U_Bi_i_QCalc.cpp  (complete Universal Buoyancy equation catalogue)

Helper modules:
    gen_fubiicalc_header.py   — C++ header, includes, struct, printSection()
    gen_fubiicalc_secA.py     — Section A: 29 per-system g_UQFF equations
    gen_fubiicalc_secB.py     — Section B: Triadic Master equations
    gen_fubiicalc_secC.py     — Section C: Sub-equations
    gen_fubiicalc_secD.py     — Section D: F_U_Bi_i component forces
    gen_fubiicalc_secE.py     — Section E: 79 F_UBii variant types
    gen_fubiicalc_secF.py     — Section F: 68 Um variant types
    gen_fubiicalc_secG.py     — Section G: Numerical constants
    gen_fubiicalc_secH.py     — Section H: Lambda-CDM/MOND notes + main()
"""
import sys
import os

# ---------------------------------------------------------------------------
# Import all section generators
# ---------------------------------------------------------------------------
try:
    from gen_fubiicalc_header import get_header
    from gen_fubiicalc_secA   import get_section_A
    from gen_fubiicalc_secB   import get_section_B
    from gen_fubiicalc_secC   import get_section_C
    from gen_fubiicalc_secD   import get_section_D
    from gen_fubiicalc_secE   import get_section_E
    from gen_fubiicalc_secF   import get_section_F
    from gen_fubiicalc_secG   import get_section_G
    from gen_fubiicalc_secH   import get_section_H
except ImportError as e:
    print(f"[ERROR] Could not import helper module: {e}")
    print("  Ensure all gen_fubiicalc_sec*.py files are in the same directory.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Output paths
# ---------------------------------------------------------------------------
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
OUTPUT_FILE = os.path.join(SCRIPT_DIR, "F_U_Bi_i_QCalc.cpp")


def main():
    print("gen_fubiicalc.py  —  Generating F_U_Bi_i_QCalc.cpp ...")

    # Build all sections
    sections = [
        ("Header",    get_header()),
        ("Section A", get_section_A()),
        ("Section B", get_section_B()),
        ("Section C", get_section_C()),
        ("Section D", get_section_D()),
        ("Section E", get_section_E()),
        ("Section F", get_section_F()),
        ("Section G", get_section_G()),
        ("Section H", get_section_H()),
    ]

    # Report per-section sizes
    for name, content in sections:
        lines = content.count("\n")
        print(f"  {name:<12} {lines:5d} lines  ({len(content):,} bytes)")

    # Concatenate
    cpp_source = "".join(content for _, content in sections)

    # Write output
    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        f.write(cpp_source)

    total_lines = cpp_source.count("\n")
    total_bytes = len(cpp_source.encode("utf-8"))
    print(f"\n  [OK] Written -> {OUTPUT_FILE}")
    print(f"       {total_lines:,} lines   {total_bytes:,} bytes")

    # Count catalogued equation entries in Section E and F data
    try:
        from gen_fubiicalc_secE import _E_DATA
        from gen_fubiicalc_secF import _F_DATA
        print(f"\n  Catalogue statistics:")
        print(f"    Section E (F_UBii variants): {len(_E_DATA):3d} entries")
        print(f"    Section F (Um variants):     {len(_F_DATA):3d} entries")
    except Exception:
        pass

    print("\n  Done. Compile and run with:")
    print("    cl /EHsc /std:c++20 F_U_Bi_i_QCalc.cpp /Fe:F_U_Bi_i_QCalc.exe")
    print("    F_U_Bi_i_QCalc.exe")


if __name__ == "__main__":
    main()
