"""
_phase_e7_merge_dispatch_into_master_closures.py
Phase E7 — Re-run Phase 2 merge with 3 new schema columns:
  geometry_used, numeric_system, assimilation_status

Extends master_closures.csv from 13 -> 16 columns. Backfills the 3 new
columns from assimilation_dispatch.py for every row that matches a dispatch
entry. Untagged rows keep empty values (they remain valid historical
closures awaiting future solver-bus routing).

Match strategy (preferred order):
  1. exact (csv.script == dispatch.session_script) AND closure-token match
  2. csv.closure or csv.label contains the dispatch observable name (normalized)
  3. otherwise: leave untagged

Assimilation status mapping:
  - notes contain "OPEN_QUESTION" -> TENSION
  - residual_pct < 1e-9              -> EXACT
  - otherwise                         -> OK
"""
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import assimilation_dispatch as ad

MCSV = ROOT / "master_closures.csv"

NEW_COLS = ["geometry_used", "numeric_system", "assimilation_status"]


def _norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", (s or "").lower())


def _status_for(rec) -> str:
    notes = rec.get("notes", "") or ""
    if "OPEN_QUESTION" in notes:
        return "TENSION"
    rp = rec.get("residual_pct", 0.0)
    try:
        if abs(float(rp)) < 1e-9:
            return "EXACT"
    except (TypeError, ValueError):
        pass
    return "OK"


def main():
    # Build dispatch lookups
    by_script = {}
    by_norm_name = {}
    for obs, rec in ad.DISPATCH.items():
        s = (rec.get("session_script") or "").strip()
        if s:
            by_script.setdefault(s, []).append((obs, rec))
        by_norm_name[_norm(obs)] = (obs, rec)

    # Read existing CSV
    with MCSV.open(encoding="utf-8", newline="") as f:
        rdr = csv.DictReader(f)
        existing_fields = list(rdr.fieldnames)
        rows = list(rdr)

    if any(c in existing_fields for c in NEW_COLS):
        print(f"ERROR: one of {NEW_COLS} already present in CSV. Aborting.")
        return 2

    new_fields = existing_fields + NEW_COLS
    print(f"Existing schema: {len(existing_fields)} cols")
    print(f"New schema:      {len(new_fields)} cols")
    print(f"Rows to process: {len(rows)}")

    tagged_via_script = 0
    tagged_via_name = 0
    tag_records = []  # (row_idx, observable, geometry, status, match_method)

    for i, row in enumerate(rows):
        matched = None
        method = None

        # Strategy 1: exact script match
        s = (row.get("script") or "").strip()
        if s in by_script:
            # If multiple dispatch entries on same script, try closure-token match
            candidates = by_script[s]
            closure_norm = _norm(row.get("closure", "") + " " + row.get("label", ""))
            for obs, rec in candidates:
                if _norm(obs) in closure_norm:
                    matched = (obs, rec)
                    method = "script+name"
                    break
            if matched is None:
                # Pick the only candidate, or the first
                matched = candidates[0]
                method = "script_only"

        # Strategy 2: normalized name appears in closure or label
        if matched is None:
            closure_norm = _norm(row.get("closure", ""))
            label_norm = _norm(row.get("label", ""))
            for obs_norm, (obs, rec) in by_norm_name.items():
                if len(obs_norm) < 6:
                    continue  # too short, risks false positives
                if obs_norm in closure_norm or obs_norm in label_norm:
                    matched = (obs, rec)
                    method = "name_substring"
                    break

        if matched is None:
            # untagged
            for c in NEW_COLS:
                row[c] = ""
            continue

        obs, rec = matched
        row["geometry_used"] = rec["owner_geometry"]
        row["numeric_system"] = "numerical"
        row["assimilation_status"] = _status_for(rec)
        if method.startswith("script"):
            tagged_via_script += 1
        else:
            tagged_via_name += 1
        tag_records.append((i, obs, rec["owner_geometry"], row["assimilation_status"], method))

    # Write back with extended schema
    with (MCSV.with_suffix(".csv.tmp")).open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=new_fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    total_tagged = tagged_via_script + tagged_via_name
    untagged = len(rows) - total_tagged

    print()
    print("=" * 70)
    print("PHASE E7 MERGE RESULT")
    print("=" * 70)
    print(f"Tagged via script match:    {tagged_via_script}")
    print(f"Tagged via name substring:  {tagged_via_name}")
    print(f"TOTAL tagged:               {total_tagged}")
    print(f"Untagged (historical):      {untagged}")
    print(f"Dispatch catalog size:      {len(ad.DISPATCH)}")
    print()

    # Status distribution
    status_counts = {}
    geom_counts = {}
    for _, _, g, s, _ in tag_records:
        status_counts[s] = status_counts.get(s, 0) + 1
        geom_counts[g] = geom_counts.get(g, 0) + 1
    print(f"Tagged status:    {dict(sorted(status_counts.items()))}")
    print(f"Tagged geometry:  {dict(sorted(geom_counts.items()))}")
    print()

    # Write a tag audit log for transparency
    audit = ROOT / "phase_e7_tag_audit.csv"
    with audit.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row_index", "csv_closure", "dispatch_observable", "geometry_used",
                    "assimilation_status", "match_method"])
        for i, obs, g, s, m in tag_records:
            w.writerow([i, rows[i].get("closure", ""), obs, g, s, m])
    print(f"Audit log: {audit.name} ({len(tag_records)} rows)")

    import shutil
    shutil.copyfile(str(MCSV.with_suffix(".csv.tmp")), str(MCSV))
    try:
        MCSV.with_suffix(".csv.tmp").unlink()
    except Exception:
        pass
    return 0


if __name__ == "__main__":
    sys.exit(main())
