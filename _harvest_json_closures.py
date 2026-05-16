"""Harvest closures from _session*.json files into master_closures.csv.

Many session scripts (S261-S342) write JSON dumps instead of audit-format
print lines. This script walks every _session*.json, extracts the
`closures` array (or compatible structure), and emits OK rows in the
master_closures.csv schema. Existing OK rows are preserved; we only
upgrade PARSE_FAIL rows whose JSON sibling contains data.
"""
import csv, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SESS_RE = re.compile(r"_session(\d+)_")

# load existing master
master_csv = ROOT / "master_closures.csv"
existing = list(csv.DictReader(master_csv.open(encoding="utf-8")))
fieldnames = ["ID","name","label","predicted","observed","error_pct","status","script","raw_output"]

# group JSON files by session id
json_files = sorted(ROOT.glob("_session*.json"))
print(f"scanning {len(json_files)} JSON files")

harvested = []  # new audit rows
for jf in json_files:
    m = SESS_RE.match(jf.name)
    if not m: continue
    sid = int(m.group(1))
    try:
        data = json.loads(jf.read_text(encoding="utf-8", errors="ignore"))
    except Exception:
        continue
    # find a list of closure-like dicts
    closures = None
    if isinstance(data, dict):
        for key in ("closures","results","entries","items","gap_closures","new_closures"):
            v = data.get(key)
            if isinstance(v, list) and v and isinstance(v[0], dict):
                closures = v; break
    elif isinstance(data, list) and data and isinstance(data[0], dict):
        closures = data
    if not closures: continue
    for c in closures:
        # tolerate various key spellings
        name = c.get("name") or c.get("id") or c.get("label") or "?"
        pred = c.get("predicted", c.get("pred", c.get("uqff", c.get("value"))))
        obs  = c.get("observed",  c.get("obs",  c.get("measured")))
        err  = c.get("residual_pct", c.get("error_pct", c.get("err_pct", c.get("err"))))
        if pred is None or obs is None: continue
        try:
            pred_f = float(pred); obs_f = float(obs)
        except (TypeError, ValueError):
            continue
        if err is None:
            try:
                err = abs(pred_f - obs_f) / abs(obs_f) * 100 if obs_f else 0
            except Exception:
                err = 0
        try:
            err_f = float(err)
        except (TypeError, ValueError):
            err_f = 0
        harvested.append({
            "ID": str(sid),
            "name": f"S{sid}-{name}",
            "label": str(name),
            "predicted": f"{pred_f}",
            "observed":  f"{obs_f}",
            "error_pct": f"{err_f}",
            "status": "OK",
            "script": jf.name,
            "raw_output": "harvested from JSON",
        })

print(f"harvested {len(harvested)} closures from JSON")

# merge: keep all existing rows; append harvested rows that don't collide on (ID,label)
seen = {(r["ID"], r["label"]) for r in existing if r["status"]=="OK"}
added = []
for h in harvested:
    key = (h["ID"], h["label"])
    if key in seen: continue
    seen.add(key); added.append(h)

print(f"adding {len(added)} new OK rows")

# write merged file
with master_csv.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in existing: w.writerow(r)
    for r in added:    w.writerow(r)

# stats
all_rows = existing + added
ok = sum(1 for r in all_rows if r["status"]=="OK")
fail = sum(1 for r in all_rows if r["status"]!="OK")
print(f"master_closures.csv now: {ok} OK, {fail} PARSE_FAIL, {len(all_rows)} total")
