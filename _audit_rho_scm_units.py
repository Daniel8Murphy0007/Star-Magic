"""Audit rho_SCm unit conventions (kg/m^3 vs J/m^3) across whitepapers.

For every PAPER_*.md, find lines mentioning 7.09e-37 / 7.09\\times10^{-37} / rho_SCm
and capture the unit label within +/- 80 chars.
"""
from __future__ import annotations
import re, os, json, sys
from collections import Counter, defaultdict

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "whitepapers")

# Regex for the 7.09e-37 magnitude in many notations
NUM_RE = re.compile(
    r"7\.09\s*(?:\\?times|\\cdot|\*|x|\\\\?times)?\s*10\s*[\^\{]?\s*-?\s*[-\u2212]?\s*37",
    re.IGNORECASE,
)
# Also catch 7.09e-37 / 7.09E-37
SCI_RE = re.compile(r"7\.09\s*[eE]\s*-\s*37")
# Unit detectors
UNIT_PATTERNS = [
    ("J/m^3", re.compile(r"J\s*/\s*m[}\s\$]*\^?\s*\{?\s*3|\\text\{J/m\}\^?3|\\mathrm\{J/m\^3\}|J\\,?m\^?\{?-3|J/m3|J\\cdot m\^?-3", re.IGNORECASE)),
    ("kg/m^3", re.compile(r"kg\s*/\s*m[}\s\$]*\^?\s*\{?\s*3|\\text\{kg/m\}\^?3|\\mathrm\{kg/m\^3\}|kg\\,?m\^?\{?-3|kg/m3", re.IGNORECASE)),
    ("eV/m^3", re.compile(r"eV\s*/\s*m\s*\^?\s*\{?\s*3", re.IGNORECASE)),
    ("dimensionless", re.compile(r"dimensionless", re.IGNORECASE)),
]

def classify_units(context: str):
    hits = []
    for name, rx in UNIT_PATTERNS:
        if rx.search(context):
            hits.append(name)
    return hits

def scan_file(path: str):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            text = f.read()
    except OSError:
        return []
    results = []
    for m in list(NUM_RE.finditer(text)) + list(SCI_RE.finditer(text)):
        start = max(0, m.start() - 100)
        end = min(len(text), m.end() + 120)
        ctx = text[start:end].replace("\n", " ")
        units = classify_units(ctx)
        results.append({"pos": m.start(), "ctx": ctx, "units": units})
    return results

def main():
    files = sorted(f for f in os.listdir(ROOT) if f.startswith("PAPER_") and f.endswith(".md"))
    per_unit = Counter()
    per_paper_unit = defaultdict(set)
    no_unit_count = 0
    detail = []
    for fname in files:
        path = os.path.join(ROOT, fname)
        hits = scan_file(path)
        if not hits:
            continue
        for h in hits:
            if h["units"]:
                for u in h["units"]:
                    per_unit[u] += 1
                    per_paper_unit[fname].add(u)
            else:
                no_unit_count += 1
        detail.append({
            "file": fname,
            "n_hits": len(hits),
            "units_seen": sorted(set().union(*[set(h["units"]) for h in hits])) if hits else [],
        })

    print("=" * 70)
    print("RHO_SCM UNIT AUDIT (7.09e-37 magnitude in whitepapers/PAPER_*.md)")
    print("=" * 70)
    print(f"Papers with at least one 7.09e-37 hit: {len(detail)}")
    print(f"Total hits with explicit units near them: {sum(per_unit.values())}")
    print(f"Total hits WITHOUT an explicit unit label nearby: {no_unit_count}")
    print()
    print("Unit label frequency (across all hits, summed):")
    for u, n in per_unit.most_common():
        print(f"  {u:15s} {n}")
    print()

    # Papers where BOTH kg/m^3 and J/m^3 appear (internal inconsistency)
    both = [p for p, us in per_paper_unit.items() if {"kg/m^3", "J/m^3"} <= us]
    print(f"Papers containing BOTH kg/m^3 AND J/m^3 for 7.09e-37: {len(both)}")
    for p in both:
        print(f"  - {p}")
    print()

    only_kg = sorted(p for p, us in per_paper_unit.items() if us == {"kg/m^3"})
    only_J = sorted(p for p, us in per_paper_unit.items() if us == {"J/m^3"})
    print(f"Papers using ONLY kg/m^3 ({len(only_kg)}):")
    for p in only_kg[:20]:
        print(f"  - {p}")
    if len(only_kg) > 20:
        print(f"  ... +{len(only_kg)-20} more")
    print()
    print(f"Papers using ONLY J/m^3 ({len(only_J)}):")
    for p in only_J[:20]:
        print(f"  - {p}")
    if len(only_J) > 20:
        print(f"  ... +{len(only_J)-20} more")

    # Save JSON
    out = os.path.join(os.path.dirname(__file__), "_audit_rho_scm_units.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump({
            "summary": {
                "papers_with_hits": len(detail),
                "hits_with_units": sum(per_unit.values()),
                "hits_without_units": no_unit_count,
                "unit_counts": dict(per_unit),
                "papers_with_both_units": both,
                "papers_kg_only": only_kg,
                "papers_J_only": only_J,
            },
        }, f, indent=2)
    print(f"\nSaved: {out}")

if __name__ == "__main__":
    main()
