"""Append millennium-cluster closures (S296-S302) to master_closures.csv with dedup."""
import csv, re, subprocess, sys
from pathlib import Path

LINE_RE = re.compile(
    r"^(?P<name>\S+)\s*::\s*predicted\s*=\s*(?P<pred>[-+0-9eE\.\/]+)"
    r"\s+observed\s*=\s*(?P<obs>[-+0-9eE\.\/]+)"
    r"\s+error_pct\s*=\s*(?P<err>[-+0-9eE\.\/]+)"
)
SCRIPT = "_session296_302_millennium_harvester.py"
CSV    = Path("master_closures.csv")

# session ID -> originating script (for traceability)
SID_FROM_NAME = {
    "S296": "_session296_millennium_poincare_perelman.py",
    "S297": "_session297_millennium_riemann.py",
    "S298": "_session298_millennium_p_vs_np.py",
    "S299": "_session299_millennium_yang_mills.py",
    "S300": "_session300_millennium_navier_stokes.py",
    "S301": "_session301_millennium_hodge.py",
    "S302": "_session302_millennium_bsd.py",
}

def parent_script_for(name):
    prefix = name.split("_", 1)[0].upper()
    return SID_FROM_NAME.get(prefix, SCRIPT)

def parent_id_for(name):
    prefix = name.split("_", 1)[0].upper()
    return int(prefix[1:]) if prefix.startswith("S") else 0

def run_harvester():
    r = subprocess.run([sys.executable, SCRIPT], capture_output=True, text=True, timeout=30)
    if r.returncode != 0:
        sys.stderr.write(r.stderr); sys.exit(r.returncode)
    return r.stdout.splitlines()

def main():
    rows = list(csv.DictReader(CSV.open(encoding="utf-8")))
    field = list(rows[0].keys()) if rows else ["ID","name","label","predicted","observed","error_pct","status","script","raw_output"]
    existing = {(r["ID"], r["label"]) for r in rows}

    added = 0
    for line in run_harvester():
        m = LINE_RE.match(line.strip())
        if not m: continue
        name = m.group("name")
        sid  = parent_id_for(name)
        key  = (str(sid), name)
        if key in existing: continue
        pred = float(m.group("pred"))
        obs  = float(m.group("obs"))
        err  = float(m.group("err"))
        rows.append({
            "ID": str(sid),
            "name": name,
            "label": name,
            "predicted": f"{pred:.10g}",
            "observed":  f"{obs:.10g}",
            "error_pct": f"{err:.6f}",
            "status": "OK",
            "script": parent_script_for(name),
            "raw_output": line.strip(),
        })
        existing.add(key); added += 1

    with CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=field)
        w.writeheader(); w.writerows(rows)
    print(f"Appended {added} millennium-cluster closures.")
    print(f"Total rows: {len(rows)} | OK: {sum(1 for r in rows if r['status']=='OK')}")

if __name__ == "__main__":
    main()
