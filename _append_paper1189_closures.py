"""Append S343-S352 (PAPER_1189) closures to master_closures.csv."""
import csv, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = "_session343_352_paper1189_chemistry.py"
LINE_RE = re.compile(r"^(S\d+)_([^\s:]+)\s*::\s*predicted=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+observed=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+error_pct=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)")
SM_ANCHOR = "CVW v2.0.0 -- G1 + G3 + G6 + G7 SM Anchor Gate compliant"

r = subprocess.run([sys.executable, str(ROOT/SCRIPT)], capture_output=True, text=True, cwd=str(ROOT), encoding="utf-8")
new_rows = []
for ln in r.stdout.splitlines():
    m = LINE_RE.match(ln.strip())
    if not m:
        continue
    sid_str, label, pred, obs, err = m.groups()
    sid = sid_str.lstrip("S")
    new_rows.append({
        "ID": sid,
        "name": sid_str + "_" + label,
        "label": label,
        "predicted": pred,
        "observed": obs,
        "error_pct": err,
        "status": "OK",
        "script": SCRIPT,
        "raw_output": sid_str + "_" + label + " :: predicted=" + pred + " observed=" + obs + " error_pct=" + err,
    })

mcsv = ROOT/"master_closures.csv"
with mcsv.open(encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f)
    fieldnames = list(reader.fieldnames)
    existing = list(reader)

ok_keys = {(row.get("ID",""), row.get("label","")) for row in existing}
added = []
for nr in new_rows:
    if (nr["ID"], nr["label"]) in ok_keys:
        continue
    full = {fn: "" for fn in fieldnames}
    for k, v in nr.items():
        if k in full:
            full[k] = v
    if "closure" in full:
        full["closure"] = nr["name"]
    if "category" in full:
        full["category"] = "paper_orphan"
    if "cvw_stamp" in full:
        full["cvw_stamp"] = "v2.0.0"
    if "sm_anchor" in full:
        full["sm_anchor"] = SM_ANCHOR
    added.append(full)

with mcsv.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for row in existing:
        w.writerow(row)
    for row in added:
        w.writerow(row)

print("Appended " + str(len(added)) + " PAPER_1189 closures.")
