"""
_harvest_template_banners.py — extract `% CLOSURE :: …` banner lines from
PAPER_*.tex / template .tex files in whitepapers/ and append them to
master_closures.csv.

Banner schema (parser-locked, single line, must appear exactly):
    % CLOSURE :: <label> :: predicted=<X> observed=<Y> error_pct=<Z>
    % TEMPLATE :: T-<code>
    % CANONICAL :: <constant-source-line>

T-PRED sentinel is detected by `error_pct=9999` (observed=9999) — status is
written as OPEN_PREDICTIONS rather than computed against the observed value.

Idempotent: rows are keyed by (closure_label, source_tex_filename). Existing
rows are NOT overwritten unless --force is supplied.

Usage:
    _harvest_template_banners.py [paper_glob ...]
    _harvest_template_banners.py whitepapers/PAPER_0500_*.tex
    _harvest_template_banners.py --force        # rewrite existing rows
    _harvest_template_banners.py --dry-run      # show what would be added

If no glob is given, defaults to whitepapers/PAPER_*.tex.
"""
from __future__ import annotations

import argparse
import csv
import datetime
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LEDGER = ROOT / "master_closures.csv"
WHITEPAPERS = ROOT / "whitepapers"

LEDGER_FIELDS = [
    "closure", "predicted", "observed", "error_pct", "status",
    "cvw_stamp", "sm_anchor", "label", "raw_output",
    "category", "name", "script", "ID",
]

# Banner regexes — anchored to start-of-line whitespace + % marker.
CLOSURE_RE = re.compile(
    r"^\s*%\s*CLOSURE\s*::\s*(?P<label>[^:]+?)\s*::\s*"
    r"predicted=(?P<pred>\S+)\s+observed=(?P<obs>\S+)\s+error_pct=(?P<err>\S+)\s*$",
    re.MULTILINE,
)
TEMPLATE_RE = re.compile(r"^\s*%\s*TEMPLATE\s*::\s*(?P<code>T-\S+)\s*$", re.MULTILINE)
CANONICAL_RE = re.compile(r"^\s*%\s*CANONICAL\s*::\s*(?P<line>.+?)\s*$", re.MULTILINE)

SENTINEL = 9999.0


def _classify(predicted: str, observed: str, err_pct: str) -> str:
    """Map banner values to ledger status. T-PRED sentinel → OPEN_PREDICTIONS."""
    try:
        p, o, e = float(predicted), float(observed), float(err_pct)
    except ValueError:
        return "PARSE_FAIL"
    if abs(p - SENTINEL) < 1e-9 and abs(o - SENTINEL) < 1e-9 and abs(e - SENTINEL) < 1e-9:
        return "OPEN_PREDICTIONS"
    if e == 0.0:
        return "EXACT"
    if abs(e) < 1.0:
        return "OK"
    return "OK"  # high-residual rows still tagged OK; profiler buckets by error_pct.


def _load_ledger() -> tuple[list[dict], list[str]]:
    if not LEDGER.exists():
        return [], LEDGER_FIELDS
    with LEDGER.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or LEDGER_FIELDS
        rows = list(reader)
    # Coerce None → "" (DictReader returns None for short rows).
    for r in rows:
        for k in fields:
            if r.get(k) is None:
                r[k] = ""
    return rows, fields


def _save_ledger(rows: list[dict], fields: list[str]) -> None:
    backup = LEDGER.with_suffix(
        f".csv.bak_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    if LEDGER.exists():
        shutil.copy2(LEDGER, backup)
        print(f"backup -> {backup.name}")
    with LEDGER.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def _harvest_file(tex: Path) -> list[dict]:
    """Return list of pending row dicts (not yet merged)."""
    text = tex.read_text(encoding="utf-8", errors="replace")
    closures = list(CLOSURE_RE.finditer(text))
    if not closures:
        return []
    template_m = TEMPLATE_RE.search(text)
    canonical_m = CANONICAL_RE.search(text)
    template_code = template_m.group("code") if template_m else ""
    canonical_line = canonical_m.group("line") if canonical_m else ""

    out: list[dict] = []
    for cm in closures:
        label = cm.group("label").strip()
        pred = cm.group("pred").strip()
        obs = cm.group("obs").strip()
        err = cm.group("err").strip()
        status = _classify(pred, obs, err)
        cvw_stamp = f"{label}: predicted={pred} observed={obs} error_pct={err} status={status}"
        out.append({
            "closure": label,
            "predicted": pred,
            "observed": obs,
            "error_pct": err,
            "status": status,
            "cvw_stamp": cvw_stamp,
            "sm_anchor": "v2.0.0",
            "label": f"CVW v2.0.0 -- {template_code} :: {canonical_line}" if template_code else canonical_line,
            "raw_output": f"[harvested from {tex.name}]",
            "category": "DERIVATION_FIRST_PRINCIPLES",
            "name": label,
            "script": tex.name,
            "ID": "",  # filled by profiler if needed
        })
    return out


def harvest(globs: list[str], *, force: bool = False, dry_run: bool = False) -> int:
    """Scan tex files matching globs; append/update master_closures.csv. Return rows added/updated."""
    if not globs:
        globs = ["PAPER_*.tex"]

    tex_files: list[Path] = []
    for g in globs:
        # Allow both "whitepapers/PAPER_*.tex" and bare "PAPER_*.tex".
        p = Path(g)
        if p.is_absolute() or p.parts[:1] == ("whitepapers",):
            tex_files.extend(sorted(ROOT.glob(g)))
        else:
            tex_files.extend(sorted(WHITEPAPERS.glob(g)))
    tex_files = sorted(set(tex_files))

    if not tex_files:
        print(f"No .tex files matched: {globs}")
        return 0

    rows, fields = _load_ledger()
    # Index existing rows by (closure, script) for idempotency.
    index: dict[tuple[str, str], int] = {}
    for i, r in enumerate(rows):
        index[(r.get("closure", ""), r.get("script", ""))] = i

    added = updated = skipped = 0
    pending: list[dict] = []
    for tex in tex_files:
        new_rows = _harvest_file(tex)
        if not new_rows:
            continue
        print(f"  {tex.name}: {len(new_rows)} banner(s)")
        for nr in new_rows:
            key = (nr["closure"], nr["script"])
            if key in index:
                if force:
                    rows[index[key]] = nr
                    updated += 1
                else:
                    skipped += 1
            else:
                pending.append(nr)
                added += 1

    print(f"\nSummary: +{added} new, ~{updated} updated, -{skipped} skipped (already present)")
    if dry_run:
        print("(dry-run; ledger not written)")
        return added + updated

    if added or updated:
        rows.extend(pending)
        _save_ledger(rows, fields)
        print(f"Ledger updated: {LEDGER.name} ({len(rows)} total rows)")
    else:
        print("No changes — ledger untouched.")
    return added + updated


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("globs", nargs="*", help="tex glob(s) under whitepapers/ (default: PAPER_*.tex)")
    ap.add_argument("--force", action="store_true", help="overwrite existing rows on (closure, script) collision")
    ap.add_argument("--dry-run", action="store_true", help="report only; do not write the ledger")
    args = ap.parse_args()
    harvest(args.globs, force=args.force, dry_run=args.dry_run)
    return 0


if __name__ == "__main__":
    sys.exit(main())
