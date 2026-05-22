"""COMPLETE inventory: gold standard, CP1-CP4, all derivation/audit files."""
from pathlib import Path
import re

ROOT = Path('.').resolve()

PATTERNS = [
    '*gold*standard*', '*Gold*Standard*', '*GOLD_STANDARD*',
    'CondensedPhysics*.py', 'CP1*', 'CP2*', 'CP3*', 'CP4*', 'CP5*',
    '*cp4*', '*CP4*',
    '*audit*gold*', '*audit*GOLD*', 'audit_gold*',
    'master_closures*', 'unified_closure*', 'sigma_table*',
    'UQFF*.py', 'UQFF*.md', 'UQFF*.json',
    'uqff_*.py', 'uqff_*.json',
    '*manifest*.json', '*registry*.json',
]

results = {}
for pat in PATTERNS:
    for p in ROOT.rglob(pat):
        if p.is_file() and '.git' not in p.parts and 'node_modules' not in p.parts:
            results[str(p)] = p.stat().st_size

# group by directory
from collections import defaultdict
by_dir = defaultdict(list)
for path, sz in results.items():
    pp = Path(path)
    rel = pp.relative_to(ROOT)
    by_dir[str(rel.parent)].append((rel.name, sz))

print(f"=== COMPLETE FILE INVENTORY ({len(results)} files) ===\n")
for d in sorted(by_dir.keys()):
    files = sorted(by_dir[d], key=lambda x: -x[1])
    print(f"[{d}/]  ({len(files)} files)")
    for name, sz in files:
        sz_str = f"{sz/1024:.1f} KB" if sz < 1e6 else f"{sz/1e6:.2f} MB"
        print(f"   {sz_str:>10}  {name}")
    print()

# Specific high-value lookups
print("=" * 60)
print("KEY FILES STATUS:")
print("=" * 60)
KEY = [
    'CondensedPhysics.py', 'CondensedPhysics2.py', 'CondensedPhysics3.py',
    'CondensedPhysics4.py', 'CP4.py', 'CP4',
    'audit_gold_standard_results.json',
    'gold_standard.json', 'gold_standard_results.json',
    'master_closures.csv', 'unified_closure_audit.json',
    'sigma_table.csv',
    'UQFF_UNIFIED_CLOSURE_DERIVATIONS.py',
    'uqff_lagrangian_derivation.py',
]
for k in KEY:
    hits = list(ROOT.rglob(k))
    hits = [h for h in hits if h.is_file() and '.git' not in h.parts]
    if hits:
        for h in hits:
            sz = h.stat().st_size
            sz_str = f"{sz/1024:.1f} KB" if sz < 1e6 else f"{sz/1e6:.2f} MB"
            print(f"  ✓ {k}: {sz_str}  ({h.relative_to(ROOT)})")
    else:
        print(f"  ✗ {k}: NOT FOUND")
