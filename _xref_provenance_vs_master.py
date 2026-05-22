"""
Cross-reference provenance_audit_v2.csv against master_closures.csv to
quantify how many audit rows derive from curve-fit scripts vs honest
structural / primitive-derivation scripts.
"""
from __future__ import annotations
import csv
from pathlib import Path
from collections import Counter

ROOT = Path(__file__).parent

prov = {}
with (ROOT / "provenance_audit_v2.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        prov[row["script"]] = row["classification"]

# Map each master_closures row to its script's class
counts_by_class = Counter()
counts_by_class_and_status = Counter()
totals_by_status = Counter()
script_to_status = {}

with (ROOT / "master_closures.csv").open(encoding="utf-8") as f:
    for row in csv.DictReader(f):
        script = row.get("script", "")
        status = row.get("status", "")
        cls = prov.get(script, "EXTRA_or_UNKNOWN")
        counts_by_class[cls] += 1
        counts_by_class_and_status[(cls, status)] += 1
        totals_by_status[status] += 1

print("=" * 70)
print("master_closures.csv rows by script provenance class:")
print("=" * 70)
for cls, n in counts_by_class.most_common():
    print(f"  {cls:24s} {n:4d}")

print(f"\nTotal rows: {sum(counts_by_class.values())}")
print(f"\nStatus totals: {dict(totals_by_status)}")

print("\n" + "=" * 70)
print("Cross-tab: provenance class x status")
print("=" * 70)
classes = sorted(counts_by_class.keys())
statuses = sorted(totals_by_status.keys())
print(f"  {'class':<24s}  " + "  ".join(f"{s:>14s}" for s in statuses))
for cls in classes:
    row = "  ".join(f"{counts_by_class_and_status[(cls, s)]:>14d}" for s in statuses)
    print(f"  {cls:<24s}  {row}")

print("\n" + "=" * 70)
print("OF THE 'EXACT' (residual=0) ROWS, BY CLASS:")
print("=" * 70)
exact_by_class = Counter()
for (cls, status), n in counts_by_class_and_status.items():
    if status in ("EXACT", "candidate-EXACT", "candidate_EXACT"):
        exact_by_class[cls] += n
for cls, n in exact_by_class.most_common():
    print(f"  {cls:24s} {n:4d}")
