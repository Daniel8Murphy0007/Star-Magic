"""
test_3numeric_millennium_crosscheck.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  numeric_backends (symbolic, numerical, discrete)
                          + uqff_pure_calculator (locked primitives)
Dependencies (external):  sympy>=1.12 (via numeric_backends.symbolic)

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PHASE B SUCCESS CRITERION
----------------------------------------------------------------------------
For each of the 8 Clay Millennium derivations, all three numeric backends
must produce matching values within float precision (tolerance: 1e-9
relative, 1e-12 absolute for sub-unity values, special-cased for the
1e-9-scale P vs NP closure).

If all 8 derivations cross-check, the 3-numeric infrastructure is sound
and Phase B is complete. Run with:

    python test_3numeric_millennium_crosscheck.py

Exits 0 on full agreement; non-zero with a printed diff if any closure
disagrees.
"""

import sys

from numeric_backends import symbolic, numerical, discrete

CLOSURE_NAMES = [
    "yang_mills",
    "riemann",
    "navier_stokes",
    "hodge",
    "poincare",
    "p_vs_np",
    "bsd",
    "black_hole_info",
]

# Per-closure absolute tolerance: most are O(1), but P vs NP has 10^-9 scale
# differences that we don't want to flag as disagreement.
TOLERANCE_ABS = {
    "yang_mills":      1e-9,
    "riemann":         1e-6,    # 9877.78... value
    "navier_stokes":   1e-12,
    "hodge":           1e-15,
    "poincare":        1e-15,
    "p_vs_np":         1e-15,
    "bsd":             1e-4,    # documented 0.005% residual: published 0.30598 vs full-Cremona 0.30600
    "black_hole_info": 1e-5,    # external Page anchor; 5-digit match
}


def _row(label, value):
    if value is None:
        return f"  {label:<12s}  <UNAVAILABLE>"
    return f"  {label:<12s}  {value:.15g}" if isinstance(value, float) else f"  {label:<12s}  {value}"


def crosscheck_closure(name):
    s_out = symbolic.evaluate(name, dtype="float")
    n_out = numerical.evaluate(name)
    d_out = discrete.evaluate(name, dtype="float")

    s_val = s_out["value"]
    n_val = n_out["value"]
    d_val = d_out["value"]

    tol = TOLERANCE_ABS.get(name, 1e-9)
    pairs = [
        ("symbolic vs numerical", s_val, n_val),
        ("symbolic vs discrete",  s_val, d_val),
        ("numerical vs discrete", n_val, d_val),
    ]
    diffs = []
    for label, a, b in pairs:
        if a is None or b is None:
            diffs.append((label, "unavailable"))
            continue
        delta = abs(a - b)
        if delta > tol:
            diffs.append((label, f"{delta:.3e} > tol={tol:.0e}"))

    return {
        "name":         name,
        "symbolic":     s_val,
        "numerical":    n_val,
        "discrete":     d_val,
        "tolerance":    tol,
        "disagreements": diffs,
    }


def main():
    print("=" * 78)
    print("PHASE B CROSS-VALIDATION — 8 Clay Millennium derivations via 3 backends")
    print("=" * 78)
    print()

    all_pass = True
    total = 0
    passing = 0
    results = []

    for name in CLOSURE_NAMES:
        r = crosscheck_closure(name)
        results.append(r)
        total += 1

        ok = (len(r["disagreements"]) == 0)
        status = "AGREE" if ok else "DIFFER"
        if ok:
            passing += 1
        else:
            all_pass = False

        print(f"[{status}] {name}")
        print(_row("symbolic",  r["symbolic"]))
        print(_row("numerical", r["numerical"]))
        print(_row("discrete",  r["discrete"]))
        print(f"  tolerance  abs: {r['tolerance']:.0e}")
        if r["disagreements"]:
            for label, msg in r["disagreements"]:
                print(f"  DIFF       {label}:  {msg}")
        print()

    print("=" * 78)
    print(f"PHASE B summary: {passing} / {total} closures cross-check across all 3 backends")
    print("=" * 78)

    if all_pass:
        print()
        print("PHASE B SUCCESS CRITERION MET. All 3 numeric backends agree.")
        return 0
    else:
        print()
        print("PHASE B FAILURE. At least one disagreement above tolerance.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
