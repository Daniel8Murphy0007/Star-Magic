"""Append S280-S286 frontier closures to master_closures.csv with dedup."""
import csv, re, subprocess, sys
from pathlib import Path

LINE_RE = re.compile(
    r"^(?P<name>\S+)\s*::\s*predicted\s*=\s*(?P<pred>[-+0-9eE\.\/]+)"
    r"\s+observed\s*=\s*(?P<obs>[-+0-9eE\.\/]+)"
    r"\s+error_pct\s*=\s*(?P<err>[-+0-9eE\.\/]+)"
)
SCRIPT = "_session280_286_frontier_harvester.py"
CSV    = Path("master_closures.csv")

PARENT = {
    "S280": "_session280_quark_closure.py",
    "S281": "_session281_neutrino_hunt.py",
    "S283": "_session283_finetune.py",
    "S285": "_session285_CKM_closure.py",
    "S286": "_session286_PMNS_closure.py",
}

def parent_for(name):
    return PARENT.get(name.split("_",1)[0].upper(), SCRIPT)

def parent_id(name):
    p = name.split("_",1)[0].upper()
    return int(p[1:]) if p.startswith("S") else 0

def main():
    rows = list(csv.DictReader(CSV.open(encoding="utf-8")))
    field = list(rows[0].keys())
    existing = {(r["ID"], r["label"]) for r in rows}

    r = subprocess.run([sys.executable, SCRIPT], capture_output=True, text=True, timeout=180)
    added = 0
    for line in r.stdout.splitlines():
        m = LINE_RE.match(line.strip())
        if not m: continue
        name = m.group("name"); sid = parent_id(name)
        key  = (str(sid), name)
        if key in existing: continue
        rows.append({
            "ID": str(sid), "name": name, "label": name,
            "predicted": f"{float(m.group('pred')):.10g}",
            "observed":  f"{float(m.group('obs')):.10g}",
            "error_pct": f"{float(m.group('err')):.6f}",
            "status": "OK", "script": parent_for(name),
            "raw_output": line.strip(),
        })
        existing.add(key); added += 1

    with CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=field)
        w.writeheader(); w.writerows(rows)
    print(f"Appended {added} frontier closures.")
    print(f"Total rows: {len(rows)} | OK: {sum(1 for r in rows if r['status']=='OK')}")

if __name__ == "__main__":
    main()
