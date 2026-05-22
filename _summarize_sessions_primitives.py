"""Summarize SESSIONS_PRIMITIVES_NORMALIZED.csv into two views:
   1. SESSIONS_PRIMITIVES_VALUE_FREQ.csv  — per primitive, value frequency across sessions
   2. SESSIONS_PRIMITIVES_WIDE.csv        — wide grid: file × primitive (most-common value per file)
"""
from __future__ import annotations
import csv
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).parent
SRC = ROOT / "SESSIONS_PRIMITIVES_NORMALIZED.csv"

def main() -> None:
    by_prim_val: dict[str, Counter[str]] = defaultdict(Counter)
    by_prim_val_files: dict[tuple[str, str], list[str]] = defaultdict(list)
    by_file_prim: dict[str, dict[str, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    primitives_seen: set[str] = set()
    files_seen: set[str] = set()

    with SRC.open(encoding="utf-8") as f:
        for row in csv.DictReader(f):
            prim = row["primitive"]
            if prim == "READ_ERROR":
                continue
            val = row["literal_form"]
            file = row["file"]
            by_prim_val[prim][val] += 1
            by_prim_val_files[(prim, val)].append(f"{file}#L{row['line']}")
            by_file_prim[file][prim][val] += 1
            primitives_seen.add(prim)
            files_seen.add(file)

    # 1. Value frequency
    freq_path = ROOT / "SESSIONS_PRIMITIVES_VALUE_FREQ.csv"
    with freq_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["primitive", "value", "count", "n_files", "example_sites"])
        for prim in sorted(primitives_seen):
            for val, cnt in by_prim_val[prim].most_common():
                sites = by_prim_val_files[(prim, val)]
                file_set = sorted({s.split("#")[0] for s in sites})
                examples = "; ".join(sites[:3])
                w.writerow([prim, val, cnt, len(file_set), examples])
    print(f"Wrote {freq_path.name}")

    # 2. Wide grid: file × primitive, dominant value (with conflict marker)
    primitives_sorted = sorted(primitives_seen)
    wide_path = ROOT / "SESSIONS_PRIMITIVES_WIDE.csv"
    with wide_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file"] + primitives_sorted)
        for file in sorted(files_seen):
            row = [file]
            for prim in primitives_sorted:
                vals = by_file_prim[file].get(prim)
                if not vals:
                    row.append("")
                elif len(vals) == 1:
                    row.append(next(iter(vals)))
                else:
                    # Multiple values for same primitive in same file -> mark
                    parts = "|".join(f"{v}x{c}" for v, c in vals.most_common())
                    row.append(f"CONFLICT[{parts}]")
            w.writerow(row)
    print(f"Wrote {wide_path.name}")

    # Console preview of dominant values per primitive
    print("\n== Dominant value per primitive (top across all sessions) ==")
    for prim in sorted(primitives_seen):
        top_val, top_cnt = by_prim_val[prim].most_common(1)[0]
        total = sum(by_prim_val[prim].values())
        n_variants = len(by_prim_val[prim])
        print(f"  {prim:12s} {top_val:>30s}  ({top_cnt}/{total} = {top_cnt/total:.0%}, {n_variants} variants)")

if __name__ == "__main__":
    main()
