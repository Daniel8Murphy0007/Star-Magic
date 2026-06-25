#!/usr/bin/env python3
"""Python <-> C++ cross-verification for UQFF.

Auto-generates a C++ harness that calls EVERY function in
uqff_exact_closures.cpp, compiles and runs it, captures every name=value
result, and compares against the Python closure of the same name.

Reports four categories:
  MATCH       Python and C++ agree (within tolerance)
  DRIFT       Both return scalars but they differ
  MISSING     C++ has a function with no matching Python closure key
  UNCALLABLE  Python closure exists but doesn't return a comparable scalar

Usage:
    python scripts/cpp_python_crosscheck.py [--tol 1e-6] [--out FILE]
"""
import argparse
import os
import re
import subprocess
import sys
import tempfile
from typing import Any

# Make calculator importable regardless of cwd
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
for _p in (os.getcwd(), _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import uqff_pure_calculator as u


def extract_cpp_functions(cpp_path: str) -> list[str]:
    """Extract names of zero-arg numeric functions from the C++ file.

    Matches `<type> <name>()` where type is double/int/float/long.
    Skips lines starting with `//` (commented-out duplicates).
    """
    pattern = re.compile(
        r"^\s*(?:double|int|float|long)\s+([a-zA-Z_]\w*)\s*\(\s*\)\s*\{",
        re.MULTILINE,
    )
    with open(cpp_path, "r", encoding="utf-8", errors="replace") as f:
        src = f.read()
    # Strip out commented lines first (preserve line structure for the regex)
    lines = src.split("\n")
    clean_lines = []
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("//"):
            clean_lines.append("")
        else:
            clean_lines.append(line)
    clean_src = "\n".join(clean_lines)
    names = pattern.findall(clean_src)
    # Dedupe preserving order
    seen = set()
    out = []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def generate_harness(fn_names: list[str], cpp_include: str) -> str:
    """Generate a C++ program that calls every function and prints name=value."""
    body = []
    for name in fn_names:
        # Use %.17g for full double precision; cast int returns to double
        body.append(
            f'    {{ double v = static_cast<double>(uqff::{name}()); '
            f'std::printf("{name}=%.17g\\n", v); }}'
        )
    return f"""\
#include <cstdio>
#include <cmath>
#include "{cpp_include}"

int main() {{
{chr(10).join(body)}
    return 0;
}}
"""


def call_python_closure(name: str) -> tuple[str, Any]:
    """Try to call the matching Python closure for a C++ function name.

    Returns ('match-type', value-or-explanation).
    """
    if not hasattr(u, "PARADOX_TO_CLOSURE"):
        return ("no-dispatch", None)

    dispatch = u.PARADOX_TO_CLOSURE
    # Strip common C++ suffixes our generator added (_v2, _v3, ...)
    base = re.sub(r"_v\d+$", "", name)
    candidates = [name, base, name.lower(), base.lower()]
    for cand in candidates:
        if cand in dispatch:
            try:
                result = u.calculate_paradox({"paradox": cand})
            except Exception as e:
                return ("python-error", repr(e)[:80])
            v = result.get("value") if isinstance(result, dict) else None
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                return ("scalar", float(v))
            if isinstance(v, tuple) and v and isinstance(v[0], (int, float)):
                return ("scalar", float(v[0]))
            if isinstance(v, dict) and "UQFF_formula_value" in v:
                try:
                    return ("scalar", float(v["UQFF_formula_value"]))
                except (TypeError, ValueError):
                    pass
            return ("dict", str(v)[:80])
    return ("not-found", None)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tol", type=float, default=1e-6,
                        help="relative tolerance for MATCH (default 1e-6)")
    parser.add_argument("--abs-tol", type=float, default=1e-12,
                        help="absolute tolerance floor (default 1e-12)")
    parser.add_argument("--cpp", default="uqff_exact_closures.cpp",
                        help="path to the C++ file (default: ./uqff_exact_closures.cpp)")
    parser.add_argument("--out", default="CPP_PYTHON_CROSSCHECK_REPORT.md",
                        help="output report path")
    parser.add_argument("--limit", type=int, default=0,
                        help="limit C++ function count for quick testing (default: 0 = all)")
    args = parser.parse_args()

    cpp_path = os.path.abspath(args.cpp)
    if not os.path.isfile(cpp_path):
        print(f"ERROR: {cpp_path} not found", file=sys.stderr)
        return 1

    print(f"Scanning {cpp_path} ...")
    fn_names = extract_cpp_functions(cpp_path)
    if args.limit:
        fn_names = fn_names[: args.limit]
    print(f"  Extracted {len(fn_names)} C++ function names")

    # Generate harness in same dir as cpp file so the #include works
    cpp_dir = os.path.dirname(cpp_path) or "."
    cpp_basename = os.path.basename(cpp_path)
    harness_path = os.path.join(cpp_dir, "_crosscheck_harness.cpp")
    binary_path = os.path.join(tempfile.gettempdir(), "uqff_crosscheck_bin")

    print(f"Generating harness at {harness_path} ...")
    with open(harness_path, "w", encoding="utf-8") as f:
        f.write(generate_harness(fn_names, cpp_basename))

    print(f"Compiling ...")
    compile_cmd = ["g++", "-std=c++17", "-O0", "-w", harness_path, "-o", binary_path]
    cp = subprocess.run(compile_cmd, capture_output=True, text=True)
    if cp.returncode != 0:
        print(f"  COMPILE FAILED:\n{cp.stderr[:3000]}", file=sys.stderr)
        os.remove(harness_path)
        return 2
    print(f"  OK -> {binary_path}")

    print(f"Running C++ harness to capture all {len(fn_names)} values ...")
    rp = subprocess.run([binary_path], capture_output=True, text=True, timeout=60)
    cpp_values: dict[str, float] = {}
    for line in rp.stdout.splitlines():
        if "=" in line:
            n, _, vstr = line.partition("=")
            try:
                cpp_values[n.strip()] = float(vstr.strip())
            except ValueError:
                pass
    print(f"  Captured {len(cpp_values)} numeric values from C++")

    # Now compare each to Python
    print(f"Comparing to Python closures ...")
    results = {"MATCH": [], "DRIFT": [], "MISSING": [], "UNCALLABLE": []}
    for name in fn_names:
        if name not in cpp_values:
            continue
        cv = cpp_values[name]
        kind, pv = call_python_closure(name)
        if kind == "scalar":
            denom = max(abs(pv), args.abs_tol)
            rel_err = abs(cv - pv) / denom
            if rel_err <= args.tol:
                results["MATCH"].append((name, cv, pv, rel_err))
            else:
                results["DRIFT"].append((name, cv, pv, rel_err))
        elif kind == "not-found":
            results["MISSING"].append((name, cv))
        else:
            results["UNCALLABLE"].append((name, cv, kind, pv))

    # Cleanup
    try: os.remove(harness_path)
    except OSError: pass
    try: os.remove(binary_path)
    except OSError: pass

    # Write report
    print(f"Writing {args.out} ...")
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(f"# C++ vs Python cross-check report\n\n")
        f.write(f"**C++ source:** `{args.cpp}`\n")
        f.write(f"**Tolerance:** relative={args.tol}, absolute floor={args.abs_tol}\n")
        f.write(f"**Total C++ functions:** {len(fn_names)}\n")
        f.write(f"**Captured numeric values:** {len(cpp_values)}\n\n")
        f.write(f"## Summary\n\n")
        f.write(f"| Category | Count | Meaning |\n|---|---:|---|\n")
        f.write(f"| MATCH      | {len(results['MATCH']):>5} | Python and C++ agree within tolerance |\n")
        f.write(f"| DRIFT      | {len(results['DRIFT']):>5} | Both return scalars but they differ |\n")
        f.write(f"| MISSING    | {len(results['MISSING']):>5} | C++ function has no matching Python closure key |\n")
        f.write(f"| UNCALLABLE | {len(results['UNCALLABLE']):>5} | Python closure exists but doesn't return a comparable scalar |\n\n")

        if results["DRIFT"]:
            f.write(f"## DRIFT — {len(results['DRIFT'])} value mismatches\n\n")
            f.write(f"These need investigation. C++ and Python disagree for the same name.\n\n")
            f.write(f"| Function | C++ value | Python value | Relative error |\n|---|---|---|---:|\n")
            for name, cv, pv, err in sorted(results["DRIFT"], key=lambda x: -x[3])[:50]:
                f.write(f"| `{name}` | {cv:.6g} | {pv:.6g} | {err:.3e} |\n")
            if len(results["DRIFT"]) > 50:
                f.write(f"\n_({len(results['DRIFT']) - 50} more drift entries omitted)_\n")

        if results["MISSING"]:
            f.write(f"\n## MISSING — {len(results['MISSING'])} C++ functions with no Python counterpart\n\n")
            f.write(f"These are C++-only functions; no PARADOX_TO_CLOSURE key matches.\n\n")
            f.write(f"```\n")
            for name, cv in sorted(results["MISSING"])[:50]:
                f.write(f"  {name:50s}  = {cv:.6g}\n")
            if len(results["MISSING"]) > 50:
                f.write(f"  ... and {len(results['MISSING']) - 50} more\n")
            f.write(f"```\n")

        if results["UNCALLABLE"]:
            f.write(f"\n## UNCALLABLE — Python closure exists but isn't a scalar\n\n")
            f.write(f"{len(results['UNCALLABLE'])} entries. Python closures returning dicts/None/structured data don't have a single number to compare.\n\n")
            f.write(f"_First 20:_\n\n```\n")
            for entry in results["UNCALLABLE"][:20]:
                name, cv = entry[0], entry[1]
                kind = entry[2] if len(entry) > 2 else "?"
                f.write(f"  {name:40s}  C++={cv:<12.6g}  python={kind}\n")
            f.write(f"```\n")

        f.write(f"\n## Verdict\n\n")
        n_compared = len(results["MATCH"]) + len(results["DRIFT"])
        if n_compared:
            match_pct = 100.0 * len(results["MATCH"]) / n_compared
            f.write(f"Of {n_compared} comparable entries, **{match_pct:.1f}% match** within "
                    f"relative tolerance {args.tol}.\n\n")
        f.write(f"DRIFT entries are the high-priority audit targets. Each represents a "
                f"case where the C++ value encoded does not match the current Python closure "
                f"value. Likely causes: (a) the Python closure formula changed after the C++ "
                f"port was authored; (b) the C++ value was a target rather than a derived "
                f"value; (c) a typo. Verify each manually and either update the C++ value or "
                f"flag the closure as deprecated.\n")

    # Console summary
    print()
    print(f"  MATCH      = {len(results['MATCH'])}")
    print(f"  DRIFT      = {len(results['DRIFT'])}")
    print(f"  MISSING    = {len(results['MISSING'])}")
    print(f"  UNCALLABLE = {len(results['UNCALLABLE'])}")
    n_cmp = len(results["MATCH"]) + len(results["DRIFT"])
    if n_cmp:
        print(f"\n  Match rate: {100.0 * len(results['MATCH']) / n_cmp:.1f}% of {n_cmp} comparable entries")
    print(f"\n  Report saved to: {args.out}")
    return 0 if not results["DRIFT"] else 3


if __name__ == "__main__":
    sys.exit(main())
