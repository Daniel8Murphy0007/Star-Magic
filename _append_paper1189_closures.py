"""Append S343-S352 (PAPER_1189) closures to master_closures.csv."""
import csv, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = "_session343_352_paper1189_chemistry.py"
LINE_RE = re.compile(r"^(S\d+)_([^\s:]+)\s*::\s*predicted=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+observed=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+error_pct=([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)")

r = subprocess.run([sys.executable, str(ROOT/SCRIPT)], capture_output=True, text=True, cwd=str(ROOT), encoding="utf-8")
new_rows = []
for ln in r.stdout.splitlines():
    m = LINE_RE.match(ln.strip())
    if not m: continue
    sid_str, label, pred, obs, err = m.groups()
    sid = sid_str.lstrip("S")
    new_rows.append({
        "ID": sid, "name": f"{sid_str}_{label}", "label": label,
        "predicted": pred, "observed": obs, "error_pct": err,
        "status": "OK", "script": SCRIPT,
        "raw_output": f"{sid_str}_{label} :: predicted={pred} observed={obs} error_pct={err}",
    })

mcsv = ROOT/"master_closures.csv"
existing = list(csv.DictReader(mcsv.open(encoding="utf-8")))
fieldnames = ["ID","name","label","predicted","observed","error_pct","status","script","raw_output"]
ok_keys = {(r["ID"], r["label"]) for r in existing}
added = [r for r in new_rows if (r["ID"], r["label"]) not in ok_keys]

with mcsv.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames); w.writeheader()
    for r in existing: w.writerow(r)
    for r in added:    w.writerow(r)

print(f"Appended {len(added)} PAPER_1189 closures.")
