"""Forensic audit: what's actually in the 478 _session*.py files?"""
import csv
import os
from pathlib import Path

ROOT = Path('.')

# 1) Count + size distribution of _session*.py
sess = sorted(ROOT.glob('_session*.py'))
print(f"=== _session*.py FORENSICS ===")
print(f"Total files: {len(sess)}")

sizes = [(p.name, p.stat().st_size) for p in sess]
tiny = [s for s in sizes if s[1] < 500]
small = [s for s in sizes if 500 <= s[1] < 2000]
medium = [s for s in sizes if 2000 <= s[1] < 10000]
large = [s for s in sizes if s[1] >= 10000]
print(f"  <500 bytes  (likely stubs/empty): {len(tiny)}")
print(f"  500-2K  (one-liner derivations):  {len(small)}")
print(f"  2-10K   (real script):            {len(medium)}")
print(f"  >=10K   (substantial):            {len(large)}")

print(f"\nSmallest 5 examples:")
for n, s in sorted(sizes, key=lambda x: x[1])[:5]:
    print(f"  {s:>6}B  {n}")

print(f"\nLargest 5 examples:")
for n, s in sorted(sizes, key=lambda x: x[1], reverse=True)[:5]:
    print(f"  {s:>6}B  {n}")

# 2) Quick content check: do they contain real Python or are they stubs?
print(f"\n=== CONTENT SAMPLING ===")
import random
random.seed(42)
samples = random.sample(sess, 5)
for p in samples:
    txt = p.read_text(encoding='utf-8', errors='ignore')
    has_print = 'print(' in txt
    has_math = 'math.' in txt or 'numpy' in txt or 'Fraction' in txt
    has_def = 'def ' in txt
    line_count = len(txt.splitlines())
    print(f"  {p.name}: {line_count} lines, print={has_print}, math={has_math}, def={has_def}")

# 3) master_closures.csv - does each row reference a script that EXISTS?
print(f"\n=== master_closures.csv -> script reference audit ===")
mc = ROOT / 'master_closures.csv'
exist_count = 0
missing_count = 0
no_script_field = 0
script_examples_missing = []
with open(mc, encoding='utf-8', newline='') as f:
    reader = csv.DictReader(f)
    rows = list(reader)
print(f"  Total rows: {len(rows)}")
print(f"  Columns: {reader.fieldnames}")
for r in rows:
    script_ref = (r.get('script') or '').strip()
    if not script_ref:
        no_script_field += 1
        continue
    # Try to resolve as a file path
    if Path(script_ref).exists():
        exist_count += 1
    else:
        missing_count += 1
        if len(script_examples_missing) < 5:
            script_examples_missing.append((r['ID'], script_ref))
print(f"  Rows whose 'script' field points to existing file: {exist_count}")
print(f"  Rows whose 'script' field points to MISSING file: {missing_count}")
print(f"  Rows with empty 'script' field:                    {no_script_field}")
if script_examples_missing:
    print(f"  Missing examples:")
    for i, s in script_examples_missing:
        print(f"    ID={i}  script={s!r}")
