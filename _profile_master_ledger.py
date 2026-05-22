"""Read-only profile of master_closures.csv + cross-ref with session primitives.

Outputs:
  - stdout: summary tables (category, status, error distribution, top-N scripts)
  - MASTER_LEDGER_BY_CATEGORY.csv
  - MASTER_LEDGER_BY_STATUS.csv
  - MASTER_LEDGER_BY_SCRIPT.csv
  - LEDGER_VS_PRIMITIVES_XREF.csv  (which ledger scripts also declare primitives)
"""
from __future__ import annotations
import csv, math, re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
LEDGER = ROOT / "master_closures.csv"
NORMALIZED = ROOT / "SESSIONS_PRIMITIVES_NORMALIZED.csv"

# ---------- load ledger ----------
rows = list(csv.DictReader(LEDGER.open(encoding="utf-8")))
# Normalize: DictReader returns None for missing columns on short rows; coerce to "".
for _r in rows:
    for _k, _v in list(_r.items()):
        if _v is None:
            _r[_k] = ""
print(f"[ledger] {LEDGER.name}: {len(rows)} rows\n")

def safe_float(s):
    try: return float(s)
    except (TypeError, ValueError): return None

# ---------- (a1) by category ----------
cat_counter = Counter(r["category"] for r in rows)
cat_status = defaultdict(Counter)
for r in rows:
    cat_status[r["category"]][r["status"]] += 1

print("=== BY CATEGORY ===")
with (ROOT / "MASTER_LEDGER_BY_CATEGORY.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["category","count","status_breakdown"])
    for cat, n in cat_counter.most_common():
        sb = ";".join(f"{s}={c}" for s,c in cat_status[cat].most_common())
        cat_disp = cat if cat else "<uncategorized>"
        w.writerow([cat_disp, n, sb])
        print(f"  {n:5d}  {cat_disp:45s}  {sb}")

# ---------- (a2) by status ----------
status_counter = Counter(r["status"] for r in rows)
print("\n=== BY STATUS ===")
with (ROOT / "MASTER_LEDGER_BY_STATUS.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["status","count","pct"])
    total = len(rows)
    for st, n in status_counter.most_common():
        pct = 100.0*n/total
        w.writerow([st, n, f"{pct:.2f}"])
        print(f"  {n:5d}  ({pct:5.1f}%)  {st}")

# ---------- (a3) error distribution ----------
buckets = [("EXACT (0)",0,0), ("<1e-6",0,1e-6), ("<1e-3",1e-6,1e-3),
           ("<1%",1e-3,1.0), ("<10%",1.0,10.0), (">=10%",10.0,math.inf)]
hist = Counter()
worst = []
SENTINEL_PCT = 9999.0
for r in rows:
    e = safe_float(r["error_pct"])
    if e is None: hist["NaN"] += 1; continue
    a = abs(e)
    # Skip T-PRED sentinels (predicted=observed=error_pct=9999, status OPEN_PREDICTIONS)
    # from both the error-distribution histogram and the worst-residual ranking.
    if r.get("status","").strip() == "OPEN_PREDICTIONS" and abs(a - SENTINEL_PCT) < 1e-6:
        hist["OPEN_PREDICTIONS (T-PRED sentinel)"] += 1
        continue
    for name, lo, hi in buckets:
        if name == "EXACT (0)" and a == 0: hist[name]+=1; break
        if lo < a <= hi: hist[name]+=1; break
    worst.append((a, int(r["ID"]) if str(r.get("ID","")).strip().isdigit() else -1, r))
worst.sort(key=lambda x: (-x[0], x[1]))

print("\n=== ERROR DISTRIBUTION (|error_pct|) ===")
for name,_,_ in buckets:
    print(f"  {hist[name]:5d}  {name}")
if hist["NaN"]: print(f"  {hist['NaN']:5d}  NaN/non-numeric")
if hist["OPEN_PREDICTIONS (T-PRED sentinel)"]:
    print(f"  {hist['OPEN_PREDICTIONS (T-PRED sentinel)']:5d}  OPEN_PREDICTIONS (T-PRED sentinel, excluded from ranking)")

print("\n=== TOP 15 WORST RESIDUALS (sentinels excluded) ===")
for a, _id, r in worst[:15]:
    print(f"  {a:12.4e}%  ID={r['ID']:>4}  {r['status']:<10}  {r['label'][:55]}")

# ---------- (a4) by script ----------
script_counter = Counter(r["script"] for r in rows)
script_status = defaultdict(Counter)
for r in rows:
    script_status[r["script"]][r["status"]] += 1
with (ROOT / "MASTER_LEDGER_BY_SCRIPT.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f); w.writerow(["script","count","status_breakdown"])
    for sc, n in script_counter.most_common():
        sb = ";".join(f"{s}={c}" for s,c in script_status[sc].most_common())
        w.writerow([sc, n, sb])

print(f"\n=== SCRIPTS ===  total={len(script_counter)} unique scripts contributing")
print("Top 10 contributors:")
for sc, n in script_counter.most_common(10):
    print(f"  {n:4d}  {sc}")

# ---------- (b) cross-ref with normalized primitives ----------
SESSION_RE = re.compile(r"_session(\d+)_")
ledger_sessions = Counter()
ledger_scripts_by_session = defaultdict(set)
for r in rows:
    m = SESSION_RE.search(r["script"])
    if m:
        sid = int(m.group(1))
        ledger_sessions[sid] += 1
        ledger_scripts_by_session[sid].add(r["script"])

prim_sessions = Counter()
prim_files_by_session = defaultdict(set)
prim_by_session_primitive = defaultdict(set)  # sid -> {primitives}
if NORMALIZED.exists():
    for r in csv.DictReader(NORMALIZED.open(encoding="utf-8")):
        m = SESSION_RE.search(r["file"])
        if not m: continue
        sid = int(m.group(1))
        prim_sessions[sid] += 1
        prim_files_by_session[sid].add(r["file"])
        prim_by_session_primitive[sid].add(r["primitive"])
else:
    print(f"\n[warn] {NORMALIZED.name} missing — skipping cross-ref")

all_sids = sorted(set(ledger_sessions) | set(prim_sessions))
xref_path = ROOT / "LEDGER_VS_PRIMITIVES_XREF.csv"
with xref_path.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["session","ledger_closures","primitive_decls",
                "in_ledger_only","in_primitives_only","both",
                "primitives_declared","ledger_scripts"])
    only_ledger = only_prim = both = 0
    for sid in all_sids:
        L = ledger_sessions.get(sid, 0)
        P = prim_sessions.get(sid, 0)
        if L and P: both += 1; flag = "both"
        elif L:     only_ledger += 1; flag = "ledger_only"
        else:       only_prim += 1; flag = "primitives_only"
        w.writerow([sid, L, P,
                    1 if flag=="ledger_only" else 0,
                    1 if flag=="primitives_only" else 0,
                    1 if flag=="both" else 0,
                    ";".join(sorted(prim_by_session_primitive.get(sid, []))),
                    ";".join(sorted(ledger_scripts_by_session.get(sid, [])))])

print("\n=== CROSS-REF: LEDGER vs SESSION PRIMITIVES ===")
print(f"  Sessions in ledger only       : {only_ledger}")
print(f"  Sessions in primitives only   : {only_prim}")
print(f"  Sessions in BOTH              : {both}")
print(f"  Ledger session range          : {min(ledger_sessions):d} ... {max(ledger_sessions):d}")
if prim_sessions:
    print(f"  Primitives session range      : {min(prim_sessions):d} ... {max(prim_sessions):d}")
print(f"\n  Output: {xref_path.name}")
