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
OUTPUT_RE  = re.compile(r"([\w\-+/ ]+?):?\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*(?:vs|->|=)\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*->\s*([\d.]+)%", re.I)

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
        last = lines[-1] if lines else ""
        mm = OUTPUT_RE.search(last)
        if mm:
            label = mm.group(1).strip()
            predicted = mm.group(2)
            observed  = mm.group(3)
            err_pct   = mm.group(4)
            status    = "OK"
        else:
            label = name
            predicted = observed = err_pct = ""
            status    = "PARSE_FAIL"
        rows.append((sid, name, label, predicted, observed, err_pct, status, sp.name, last))
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
# STATUS MODE
# ---------------------------------------------------------------------------
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
    else: ap.print_help()

if __name__ == "__main__":
    main()
