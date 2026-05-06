# -*- coding: utf-8 -*-
"""
Fix all known LaTeX compilation errors across all whitepapers.

Errors fixed:
1. VDS Gompertz formula: \right\right) → \right)\right)   (132 papers)
2. exp!\left → \exp\!\left                                 (945 papers)
3. \argmin → \operatorname{argmin}                         (6 papers)
4. \Deltak → \Delta k                                      (7 papers)
5. \Lambdac → \Lambda_c                                    (5 papers)
6. \gammagamma → \gamma\gamma                              (3 papers)
7. \tauto → \tau_0                                         (1 paper)
8. \mumu → \mu\mu                                          (1 paper)
9. \nunu → \nu\nu                                          (1 paper)
10. \pir → \pi r                                           (1 paper)
11. \odotapprox → \odot\approx                             (1 paper)
12. \sumdelta → \sum\delta                                 (1 paper)
13. \epsilonsw → \epsilon_{\rm sw}                         (1 paper)
14. \langlevarepsilon → \langle\varepsilon                 (1 paper)
15. \nablaU → \nabla U                                     (4 papers)
"""

from pathlib import Path
import re

WHITEPAPER_DIR = Path("whitepapers")

# --- REPLACEMENT TABLE ---
# Each entry: (search_str, replace_str, description)
# Order matters — do VDS fix BEFORE exp!\left fix (they overlap)
FIXES = [
    # 1. VDS Gompertz: missing ) between the two \right
    #    BAD:  \exp!\left(-\exp!\left(-\frac{r-r_0}{\lambda_\mathrm{VDS}}\right\right)
    #    GOOD: \exp\!\left(-\exp\!\left(-\frac{r-r_0}{\lambda_\mathrm{VDS}}\right)\right)
    #    We handle this as a two-step: first add the missing ), then fix exp!
    (
        r'\lambda_{\mathrm{VDS}}}\right\right)',
        r'\lambda_{\mathrm{VDS}}}\right)\right)',
        "VDS Gompertz: add missing ) before \\right)",
    ),
    # 2. exp!\left → \exp\!\left  (spacing command — !\left is not LaTeX-valid spacing)
    (
        r'\exp!\left',
        r'\exp\!\left',
        "exp!\\left → \\exp\\!\\left",
    ),
    # 3. \argmin → \operatorname{argmin}
    (
        r'\argmin',
        r'\operatorname{argmin}',
        "\\argmin → \\operatorname{argmin}",
    ),
    # 4. \Deltak → \Delta k  (space needed between macro and k)
    (
        r'\Deltak',
        r'\Delta k',
        "\\Deltak → \\Delta k",
    ),
    # 5. \Lambdac → \Lambda_c
    (
        r'\Lambdac',
        r'\Lambda_c',
        "\\Lambdac → \\Lambda_c",
    ),
    # 6. \gammagamma → \gamma\gamma
    (
        r'\gammagamma',
        r'\gamma\gamma',
        "\\gammagamma → \\gamma\\gamma",
    ),
    # 7. \tauto → \tau_0
    (
        r'\tauto',
        r'\tau_0',
        "\\tauto → \\tau_0",
    ),
    # 8. \mumu → \mu\mu
    (
        r'\mumu',
        r'\mu\mu',
        "\\mumu → \\mu\\mu",
    ),
    # 9. \nunu → \nu\nu
    (
        r'\nunu',
        r'\nu\nu',
        "\\nunu → \\nu\\nu",
    ),
    # 10. \pir → \pi r  (must NOT match \pi r if already correct — check context)
    (
        r'\pir',
        r'\pi r',
        "\\pir → \\pi r",
    ),
    # 11. \odotapprox → \odot\approx
    (
        r'\odotapprox',
        r'\odot\approx',
        "\\odotapprox → \\odot\\approx",
    ),
    # 12. \sumdelta → \sum\delta
    (
        r'\sumdelta',
        r'\sum\delta',
        "\\sumdelta → \\sum\\delta",
    ),
    # 13. \epsilonsw → \epsilon_{\rm sw}
    (
        r'\epsilonsw',
        r'\epsilon_{\rm sw}',
        "\\epsilonsw → \\epsilon_{\\rm sw}",
    ),
    # 14. \langlevarepsilon → \langle\varepsilon
    (
        r'\langlevarepsilon',
        r'\langle\varepsilon',
        "\\langlevarepsilon → \\langle\\varepsilon",
    ),
    # 15. \nablaU → \nabla U
    (
        r'\nablaU',
        r'\nabla U',
        "\\nablaU → \\nabla U",
    ),
]


def fix_paper(md_path: Path) -> tuple[bool, list[str]]:
    """Apply all fixes to a single whitepaper. Returns (changed, list_of_applied_fixes)."""
    original = md_path.read_text("utf-8", errors="replace")
    text = original
    applied = []

    for search, replace, desc in FIXES:
        if search in text:
            count = text.count(search)
            text = text.replace(search, replace)
            applied.append(f"  [{count}x] {desc}")

    if text != original:
        md_path.write_text(text, encoding="utf-8")
        return True, applied
    return False, []


def main():
    papers = sorted(WHITEPAPER_DIR.glob("PAPER_*.md"))
    total_changed = 0
    fix_counts = {desc: 0 for _, _, desc in FIXES}

    for md in papers:
        changed, applied = fix_paper(md)
        if changed:
            total_changed += 1
            for line in applied:
                # Extract description from applied line
                for _, _, desc in FIXES:
                    if desc in line:
                        fix_counts[desc] += 1

    print(f"\nFixed {total_changed} / {len(papers)} whitepapers")
    print("\nFixes applied:")
    for _, _, desc in FIXES:
        count = fix_counts[desc]
        if count:
            print(f"  {count:4d} papers — {desc}")

    # Verification pass
    print("\nVerification (should all be 0):")
    remaining = {}
    for search, _, desc in FIXES:
        found = sum(1 for md in papers if search in md.read_text("utf-8", errors="replace"))
        if found:
            remaining[desc] = found
            print(f"  STILL PRESENT ({found} papers): {search!r}")
    if not remaining:
        print("  All patterns cleared ✓")


if __name__ == "__main__":
    main()
