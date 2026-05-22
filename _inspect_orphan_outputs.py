"""Run the 19 parse-fail orphans and dump the last ~15 lines of each so we
can see what closure format they emit.  Read-only — no file writes."""
import csv, subprocess, sys, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = ROOT / ".venv_py314_backup" / "Scripts" / "python.exe"
if not PY.exists():
    PY = Path(sys.executable)

fails = []
with open("master_closures.csv", encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if int(r["ID"]) >= 800 and r["status"] != "OK" and r["status"] != "OK_JSON":
            fails.append(r["script"])

env = os.environ.copy(); env["PYTHONIOENCODING"] = "utf-8"
print(f"Inspecting {len(fails)} parse-fail orphans\n")
for i, fname in enumerate(fails, 1):
    p = ROOT / fname
    print("=" * 78)
    print(f"[{i}/{len(fails)}] {fname}")
    print("=" * 78)
    try:
        r = subprocess.run([str(PY), str(p)], capture_output=True, text=True,
                           timeout=30, env=env)
        out = (r.stdout or "") + (r.stderr or "")
    except Exception as e:
        out = f"ERROR: {e}"
    lines = [l.rstrip() for l in out.splitlines() if l.strip()]
    if not lines:
        print("  (no output)")
        continue
    # Print last 15 lines (closures usually summary near end)
    for ln in lines[-15:]:
        print(f"  {ln[:160]}")
    print()
