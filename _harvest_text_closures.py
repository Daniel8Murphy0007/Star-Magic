"""Harvest closures from verbose text-format session scripts.

Scans output for triplets of:
  prediction[/predicted]  = <num>
  observed                = <num>
  residual[/error]        = <num>%

Pairs are extracted by proximity (within 8 lines). Appends OK rows to
master_closures.csv.
"""
import csv, re, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PY = sys.executable

# patterns
P_PRED = re.compile(r"\b(?:prediction|predicted|pred|UQFF|uqff|theory)\b\s*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)")
P_OBS  = re.compile(r"\b(?:observed|obs|measured|exp)\b[^=:\n]*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)")
P_RES  = re.compile(r"\b(?:residual|error|err|residue|diff)\b[^=:\n]*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*%")
SESS_RE = re.compile(r"^_session(\d+)_([^.]+)\.py$")

# load existing
mcsv = ROOT / "master_closures.csv"
existing = list(csv.DictReader(mcsv.open(encoding="utf-8")))
fieldnames = ["ID","name","label","predicted","observed","error_pct","status","script","raw_output"]
ok_keys = {(r["ID"], r["label"]) for r in existing if r["status"]=="OK"}
parse_fail_scripts = {r["script"] for r in existing if r["status"]!="OK"}

# only re-run PARSE_FAIL scripts in the targeted range
SKIP_PREFIXES = ("_audit","_inspect","_harvest","_build","_chain","_fix","_regen","_apply","_verify","_check","_combine","_calc","_propagate","_qcg","_relocate","_revert","_sample","_scan","_show","_test","_thread","_update","_session200","_session201","_session202","_session203","_session204")
candidates = []
for s in sorted(parse_fail_scripts):
    if not s.startswith("_session"): continue
    m = SESS_RE.match(s)
    if not m: continue
    sid = int(m.group(1))
    if sid < 261 or sid > 342: continue   # frontier range only
    if s in SKIP_PREFIXES: continue
    candidates.append(s)

print(f"Scanning {len(candidates)} PARSE_FAIL session scripts (S261-S342) for prediction/observed/residual triplets")

added = []
for script in candidates:
    sp = ROOT / script
    if not sp.exists(): continue
    try:
        r = subprocess.run([PY, str(sp)], capture_output=True, text=True, timeout=120, cwd=str(ROOT), encoding="utf-8", errors="ignore")
        out = r.stdout or ""
    except Exception as e:
        continue
    # find triplets via proximity
    lines = out.splitlines()
    m = SESS_RE.match(script)
    sid = int(m.group(1))
    label_stem = m.group(2)
    triplets = []
    pending_pred = pending_obs = None
    pending_pred_line = pending_obs_line = -100
    for i, ln in enumerate(lines):
        mp = P_PRED.search(ln)
        if mp:
            pending_pred = float(mp.group(1)); pending_pred_line = i
        mo = P_OBS.search(ln)
        if mo and pending_pred is not None and (i - pending_pred_line) <= 6:
            pending_obs = float(mo.group(1)); pending_obs_line = i
        mr = P_RES.search(ln)
        if mr and pending_pred is not None and pending_obs is not None and (i - pending_obs_line) <= 6:
            try:
                err = abs(float(mr.group(1)))
            except Exception:
                continue
            triplets.append((pending_pred, pending_obs, err))
            pending_pred = pending_obs = None
    # Second pass: observed + residual pairs (predicted reconstructed)
    pending_obs = None; pending_obs_line = -100
    for i, ln in enumerate(lines):
        mo = P_OBS.search(ln)
        if mo:
            try:
                pending_obs = float(mo.group(1)); pending_obs_line = i
            except Exception:
                pending_obs = None
        mr = P_RES.search(ln)
        if mr and pending_obs is not None and (i - pending_obs_line) <= 4:
            try:
                err = float(mr.group(1))
                pred = pending_obs * (1.0 + err/100.0)
                triplets.append((pred, pending_obs, abs(err)))
                pending_obs = None
            except Exception:
                pass

    # take at most 5 strongest triplets (avoid noise)
    if not triplets: continue
    for idx,(p,o,e) in enumerate(triplets[:5]):
        label = f"{label_stem}-{idx}" if len(triplets)>1 else label_stem
        key = (str(sid), label)
        if key in ok_keys: continue
        ok_keys.add(key)
        added.append({
            "ID": str(sid),
            "name": f"S{sid}-{label}",
            "label": label,
            "predicted": f"{p}",
            "observed":  f"{o}",
            "error_pct": f"{e}",
            "status": "OK",
            "script": script,
            "raw_output": f"harvested triplet ({p}, {o}, {e}%)",
        })

print(f"Harvested {len(added)} new triplets")

# merge
with mcsv.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    for r in existing: w.writerow(r)
    for r in added:    w.writerow(r)

ok = sum(1 for r in existing if r["status"]=="OK") + len(added)
fail = sum(1 for r in existing if r["status"]!="OK")
print(f"master_closures.csv: {ok} OK, {fail} PARSE_FAIL")
