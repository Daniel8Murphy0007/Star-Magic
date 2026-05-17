"""
_uqff_program.py — UQFF closure program master entry point.

Usage:
  python _uqff_program.py --audit
        Scan all _session*.py scripts, build master_closures.csv with
        {ID, name, predicted, observed, error%, expr_present, script_path}.

  python _uqff_program.py --search --tier JJ --targets targets.csv
        Brute-force search closures for a CSV of (name,value) targets.
        Writes _tier_<letter>_search.py-style report.

  python _uqff_program.py --status
        Print cumulative-counter + locking + tier summary.

LOCKED PRIMITIVES (frozen May 2026):
  F_TRZ=1/10  Phi_res=5/6  SSq=57/100  K_Mex=25/12  beta_i=6029/10000
  D_phys=4    D_BSFG=6     D_crit=26   N_ch=9       SO5=10    A_5=60
"""
from __future__ import annotations
import argparse, csv, os, re, subprocess, sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
VENV_PY = ROOT / ".venv" / "Scripts" / "python.exe"

F = Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12)
Dp=4; DB=6; Dc=26; N=9; SO=10; A=60; beta=Fraction(6029,10000)


# ---------------------------------------------------------------------------
# AUDIT MODE — build master_closures.csv from all _session*.py scripts
# ---------------------------------------------------------------------------
SESSION_RE = re.compile(r"^_session(\d+)_([^.]+)\.py$")
# Pattern A: "label: PRED vs OBS -> ERR%" or "-> EXACT"
OUTPUT_RE_A = re.compile(r"([\w\-+/ ()^.]+?):\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*vs\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*->\s*(EXACT|[\d.]+%)", re.I)
# Pattern B: "... = PRED ...; obs (...) = OBS; match ERR%"
OUTPUT_RE_B = re.compile(r"=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^;]*?;\s*obs[^=]*?=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^;]*?;\s*match\s*([\d.]+)\s*%", re.I)
# Pattern C: "label = PRED, obs = OBS, err = ERR%"  /  "label: PRED (obs OBS, err ERR%)"
OUTPUT_RE_C = re.compile(r"([\w\-+/ ()^.]+?)\s*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^,]*?,\s*obs[^=:]*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^,]*?,\s*(?:err|error|match)[^=:]*[:=]\s*([\d.]+)\s*%", re.I)

def _normalize_err(err_str: str) -> str:
    """Promote machine-zero (< 1e-9 %) errors to EXACT.

    Many closure scripts compute err_pct = abs(d-t)/t*100 which can yield
    values like 1e-14 for IEEE-754-equal results.  Treat anything below
    1 part in 10^11 (i.e. 1e-9 %) as algebraically EXACT.
    """
    try:
        v = float(err_str)
    except (TypeError, ValueError):
        return err_str
    if abs(v) < 1e-9:
        return "0"
    return err_str

def _parse_line(line, fallback_name):
    m = OUTPUT_RE_A.search(line)
    if m:
        raw = m.group(4)
        err = "0" if raw.upper()=="EXACT" else _normalize_err(raw.rstrip("%"))
        return m.group(1).strip(), m.group(2), m.group(3), err
    m = OUTPUT_RE_B.search(line)
    if m:
        return fallback_name, m.group(1), m.group(2), _normalize_err(m.group(3))
    m = OUTPUT_RE_C.search(line)
    if m:
        return m.group(1).strip(), m.group(2), m.group(3), _normalize_err(m.group(4))
    return None

def audit():
    py = str(VENV_PY) if VENV_PY.exists() else sys.executable
    env = os.environ.copy(); env["PYTHONIOENCODING"] = "utf-8"
    rows = []
    scripts = sorted(ROOT.glob("_session*.py"), key=lambda p: int(SESSION_RE.match(p.name).group(1)) if SESSION_RE.match(p.name) else 0)
    print(f"Found {len(scripts)} session scripts. Running...")
    for sp in scripts:
        m = SESSION_RE.match(sp.name)
        if not m: continue
        sid, name = m.group(1), m.group(2)
        try:
            r = subprocess.run([py, str(sp)], capture_output=True, text=True, timeout=30, env=env)
            out = (r.stdout or "") + (r.stderr or "")
        except Exception as e:
            out = f"ERROR: {e}"
        lines = [l.strip() for l in out.splitlines() if l.strip()]
        parsed = None
        for ln in reversed(lines):
            parsed = _parse_line(ln, name)
            if parsed:
                break
        if parsed:
            label, predicted, observed, err_pct = parsed
            status = "OK"
        else:
            label = name
            predicted = observed = err_pct = ""
            status    = "PARSE_FAIL"
        rows.append((sid, name, label, predicted, observed, err_pct, status, sp.name, (lines[-1] if lines else "")))
    csv_path = ROOT / "master_closures.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID","name","label","predicted","observed","error_pct","status","script","raw_output"])
        w.writerows(rows)
    ok    = sum(1 for r in rows if r[6]=="OK")
    fails = sum(1 for r in rows if r[6]!="OK")
    print(f"\nMaster registry written: {csv_path}")
    print(f"  OK:           {ok}")
    print(f"  parse fails:  {fails}")
    print(f"  total scripts: {len(rows)}")
    # sigma table — only when err_pct numeric and < 0.005 we flag candidate EXACT
    cands = [r for r in rows if r[6]=="OK" and r[5] and float(r[5])<0.0005]
    print(f"  candidate-EXACT (err<0.0005%): {len(cands)}")
    exacts = [r for r in rows if r[6]=="OK" and r[5]=="0"]
    print(f"  exact (=0%): {len(exacts)}")
    return csv_path


# ---------------------------------------------------------------------------
# SEARCH MODE — brute-force closure search
# ---------------------------------------------------------------------------
def build_pool():
    pool=[]; seen=set()
    def add(v,l):
        v=Fraction(v)
        if v not in seen: seen.add(v); pool.append((v,l))
    for base,bl in [(Dp,'Dp'),(DB,'DB'),(Dc,'Dc'),(N,'N'),(SO,'SO'),(A,'A'),
                     (SSq,'SSq'),(Phi,'Phi'),(K,'K'),(beta,'beta')]:
        for fexp,fl in [(Fraction(1),''),(F,'F'),(F*F,'F2'),(F*F*F,'F3')]:
            for p,pl in [(1,''),(2,'2'),(3,'3'),(4,'4'),(5,'5')]:
                if bl in ('SSq','Phi','beta','K') or p==1:
                    v = fexp * (Fraction(base)**p)
                    add(v, (fl+'.' if fl else '')+bl+pl)
    for k in [1,2,3,4,5]: add(k,str(k))
    return pool

def search_one(target, pool, max_terms=5, tol=0.0003):
    target=Fraction(target).limit_denominator(10**8)
    best=None
    n=len(pool)
    def rec(idx,cur,terms,uc):
        nonlocal best
        if uc>0:
            err=abs(float(cur-target)/float(target)) if target!=0 else abs(float(cur))
            if best is None or err<best[0]:
                best=(err,float(cur),list(terms))
                if err<tol: return True
        if uc>=max_terms or idx>=n: return False
        if rec(idx+1,cur,terms,uc): return True
        v,l=pool[idx]
        terms.append('+'+l)
        if rec(idx+1,cur+v,terms,uc+1): return True
        terms.pop()
        terms.append('-'+l)
        if rec(idx+1,cur-v,terms,uc+1): return True
        terms.pop()
        return False
    rec(0,Fraction(0),[],0)
    return best

def search(tier, targets_path, max_terms=5, tol=0.0003):
    pool = build_pool()
    with open(targets_path, newline="") as f:
        targets = [(row[0], float(row[1])) for row in csv.reader(f) if row and not row[0].startswith("#")]
    out_path = ROOT / f"_tier_{tier}_search.py"
    print(f"Tier {tier}: searching {len(targets)} targets (pool={len(pool)}, max_terms={max_terms})")
    print(f"Results -> {out_path}")
    lines = [f"# Auto-generated tier {tier} search results"]
    for name, t in targets:
        b = search_one(t, pool, max_terms, tol)
        line = f"{name:24s} target={t}  best={b[1]:.6g}  err={b[0]*100:.4f}%  expr={' '.join(b[2])}"
        print(line); lines.append("# " + line)
    out_path.write_text("\n".join(lines), encoding="utf-8")


# ---------------------------------------------------------------------------
# SIGMA MODE — build falsifiability table from master_closures.csv
# ---------------------------------------------------------------------------
def sigma():
    csv_p = ROOT / "master_closures.csv"
    if not csv_p.exists():
        print("master_closures.csv missing — run --audit first."); return
    rows = list(csv.DictReader(csv_p.open(encoding="utf-8")))
    ok = [r for r in rows if r["status"]=="OK" and r["error_pct"]]
    # bucket by sigma proxy: |err%|/100 vs typical experimental uncertainty
    # crude tiers (no measurement uncertainty available here):
    tiers = {"EXACT":[], "<0.01%":[], "0.01-0.1%":[], "0.1-1%":[], ">=1%":[]}
    for r in ok:
        e = float(r["error_pct"])
        if e == 0:        tiers["EXACT"].append(r)
        elif e < 0.01:    tiers["<0.01%"].append(r)
        elif e < 0.1:     tiers["0.01-0.1%"].append(r)
        elif e < 1:       tiers["0.1-1%"].append(r)
        else:             tiers[">=1%"].append(r)
    out = ROOT / "sigma_table.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["tier","ID","name","label","predicted","observed","error_pct"])
        for t,items in tiers.items():
            for r in sorted(items, key=lambda x:int(x["ID"])):
                w.writerow([t, r["ID"], r["name"], r["label"], r["predicted"], r["observed"], r["error_pct"]])
    print(f"Sigma table written: {out}")
    for t,items in tiers.items():
        print(f"  {t:12s} {len(items):4d}")
    print(f"  TOTAL OK:    {len(ok):4d}")


def predict():
    p = ROOT / "_predictive_targets.csv"
    if not p.exists():
        print("_predictive_targets.csv missing."); return
    rows = [r for r in csv.reader(p.open(encoding="utf-8")) if r and not r[0].startswith("#")]
    print(f"{'experiment':14s} {'quantity':28s} {'predicted':45s} year status")
    print("-"*120)
    for r in rows:
        if len(r) < 5: continue
        print(f"{r[0]:14s} {r[1]:28s} {r[2]:45s} {r[3]} {r[4]}")
    print(f"\nTotal predictive targets: {len(rows)}")


def status():
    scripts = list(ROOT.glob("_session*.py"))
    ids = sorted(int(SESSION_RE.match(s.name).group(1)) for s in scripts if SESSION_RE.match(s.name))
    print(f"Total session scripts: {len(ids)}")
    if ids:
        print(f"  ID range: S{ids[0]:03d} .. S{ids[-1]:03d}")
    papers = list((ROOT/"whitepapers").glob("PAPER_1209*.tex"))
    print(f"Tier papers (PAPER_1209*): {len(papers)}")
    pdfs = list((ROOT/"pdf").glob("PAPER_*.pdf"))
    print(f"PDFs in pdf/: {len(pdfs)}")
    csv_p = ROOT / "master_closures.csv"
    if csv_p.exists():
        print(f"master_closures.csv: {csv_p} ({csv_p.stat().st_size} bytes)")
    else:
        print("master_closures.csv: NOT YET BUILT (run --audit)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit", action="store_true", help="run all _session*.py and build master_closures.csv")
    ap.add_argument("--search", action="store_true", help="brute-force closure search")
    ap.add_argument("--status", action="store_true", help="print program status")
    ap.add_argument("--sigma",  action="store_true", help="build sigma_table.csv from master_closures.csv")
    ap.add_argument("--predict",action="store_true", help="print predictive-validation tracker")
    ap.add_argument("--tier", help="tier letter (e.g. JJ) for --search")
    ap.add_argument("--targets", help="path to targets CSV (name,value) for --search")
    ap.add_argument("--max-terms", type=int, default=5)
    ap.add_argument("--tol", type=float, default=0.0003)
    args = ap.parse_args()
    if args.audit:    audit()
    elif args.search:
        if not args.tier or not args.targets:
            ap.error("--search requires --tier and --targets")
        search(args.tier, args.targets, args.max_terms, args.tol)
    elif args.status: status()
    elif args.sigma:  sigma()
    elif args.predict: predict()
    else: ap.print_help()

if __name__ == "__main__":
    main()
