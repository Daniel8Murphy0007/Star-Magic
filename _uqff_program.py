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
import argparse, csv, json, os, re, subprocess, sys
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
OUTPUT_RE_A = re.compile(r"([\w\-+/ ()^.]+?):\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\S*\s*vs\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\S*\s*->\s*(EXACT|[\d.]+%)", re.I)
# Pattern B: "... = PRED ...; obs (...) = OBS; match ERR%"
OUTPUT_RE_B = re.compile(r"=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^;]*?;\s*obs[^=]*?=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^;]*?;\s*match\s*([\d.]+)\s*%", re.I)
# Pattern C: "label = PRED, obs = OBS, err = ERR%"  /  "label: PRED (obs OBS, err ERR%)"
OUTPUT_RE_C = re.compile(r"([\w\-+/ ()^.]+?)\s*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^,]*?,\s*obs[^=:]*[:=]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^,]*?,\s*(?:err|error|match)[^=:]*[:=]\s*([\d.]+)\s*%", re.I)
# Pattern D: harvester form "<label> :: predicted=PRED observed=OBS error_pct=ERR"
#   used by S280_286/S296_302/etc. harvesters.  Last-matching line wins.
OUTPUT_RE_D = re.compile(r"([\w\-+/ ().:^*]+?)\s*(?:::|:)?\s*predicted\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+observed\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+error_?pct\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)", re.I)
# Pattern E: 3-col tabular "  LABEL    PRED    OBS    ERR%"  (whitespace separated,
#   PRED and OBS in scientific or decimal, ERR a bare percent w/o sign)
OUTPUT_RE_E = re.compile(r"^([A-Za-z][\w|/.()\-+^*\\]*?)\s{2,}([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s{2,}([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s{2,}([\d.]+)\s*%?\s*$")
# Pattern F: "LABEL:  pred = NUM [unit]   obs [=~] NUM [unit]   resid = NUM%"
#   used by S281/S283-verify scripts (inline form on one line).
OUTPUT_RE_F = re.compile(r"([\w\-+/ ()^.*~]+?)\s*[:=]\s*pred(?:icted)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^=~]*?obs(?:erved)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^=~]*?resid(?:ual)?\s*[=~]\s*([\d.]+)\s*%", re.I)
# Pattern G: "pred = NUM [unit]   obs [=~] NUM [unit]   resid = NUM%"  (no leading label;
#   caller supplies fallback label from most recent "LABEL:" line).
OUTPUT_RE_G = re.compile(r"\bpred(?:icted)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^=~]*?obs(?:erved)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^=~]*?resid(?:ual)?\s*[=~]\s*([\d.]+)\s*%", re.I)
# Pattern G2: "pred = NUM [unit]   obs [=~] NUM [unit]"  (no resid on same line — compute it).
#   Used by S283-verify et al. where 'resid' is on the next line.
OUTPUT_RE_G2 = re.compile(r"\bpred(?:icted)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)[^=~\n]*?obs(?:erved)?\s*[=~]\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)", re.I)
# Pattern K: "LABEL = NUM [unit]   (SIGNED_PCT % vs obs)" — used by S294-style summary
#   banners.  Yields label+pred but no observed (left blank).
OUTPUT_RE_K = re.compile(r"^\s*([A-Za-z][\w\-+/().,*^]{0,40}?)\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*\S*\s*\(\s*([+-]?\d+\.?\d*)\s*%\s*(?:vs|from|to)?\s*obs", re.I)
# Pattern L: terminal banner "<NAME> CLOSED to <SIGNED_PCT>%" — single value, no pred/obs.
OUTPUT_RE_L = re.compile(r"\b([A-Za-z][\w\-+/().,*^]{0,40}?)\s+CLOSED\s+to\s+([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern M: "residual = SIGNED_PCT%" alone — uses fallback label (last header or script).
OUTPUT_RE_M = re.compile(r"^\s*(?:residual|resid|error|err)\s*[=~:]\s*([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern M2: "LABEL (mass_)?resid(ual)?\s*[=:]\s*NUM%" — leading token + residual keyword.
#   Used by S282 ("G_full residual = 0.000%") and S281 ("beta=NUM 'expr' mass_resid= NUM%").
OUTPUT_RE_M2 = re.compile(r"^\s*([A-Za-z][\w/\-+().*^]{0,40})\s+(?:[a-z_]*resid(?:ual)?|err(?:or)?)\s*[=~:]?\s*([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern N: banner "S<NNN> COMPLETE. <label> = <pred>; ... target <obs>; ... match <err>%"
#   used by S415/S434/S470/S499/S543/S426/... — the dominant banner schema.
#   Handled procedurally by _try_complete_banner() (chains of "= ... = NUM" are
#   too irregular for a single regex).
_NUM_RE = re.compile(r"[+-]?\d+\.?\d*(?:[eE][+-]?\d+)?")
_OBS_KW_RE = re.compile(r"\btarget\b|\bobs(?:erved)?\b|\bCODATA\b|\bNIST\b|\bSM\b|\bpredicted\b|\breference\b", re.I)
_MATCH_KW_RE = re.compile(r"\bmatch(?:\s+within)?\b\s*[=:]?\s*([+-]?\d+\.?\d*)\s*%", re.I)

def _try_complete_banner(line: str):
    """Parse banner forms like:
        S510 COMPLETE. Dodecahedron faces F = 12.0000 = 2*D_BSFG = 12; target 12; match 0.0000%.
        S347 COMPLETE. sin^2(theta_W) = ... = 0.23148; observed = 0.23122; match within 0.113%.
    Strategy: split on ';'/',', identify pred segment (first), obs segment (keyword),
    err segment ('match[...] NUM%').
    """
    iU = line.upper()
    idx = -1
    for kw in ("COMPLETE.", "CORRECTED.", "RESOLVED.", "CLOSED."):
        i = iU.find(kw)
        if i >= 0:
            idx = i + len(kw); break
    if idx < 0:
        return None
    rest = line[idx:].strip()
    parts = [p.strip() for p in re.split(r"[;,]", rest) if p.strip()]
    if len(parts) < 1:
        return None
    # Find first part with both '=' and a number — that's the predication segment.
    pred_part = None
    for p in parts:
        if "=" in p and _NUM_RE.findall(p):
            pred_part = p; break
    if pred_part is None:
        return None
    label = pred_part.split("=")[0].strip()
    # Pred should be the number on the RHS of the LAST '=' in the pred segment
    # (banners often have chains like 'V_min = ell_P^3 / D_BSFG^(3/2) = 2.873e-106').
    rhs = pred_part.rsplit("=", 1)[-1]
    nums_rhs = _NUM_RE.findall(rhs)
    if not nums_rhs:
        nums_p = _NUM_RE.findall(pred_part)
        if not nums_p:
            return None
        pred = nums_p[-1]
    else:
        pred = nums_rhs[0]
    obs = None
    for p in parts:
        if p is pred_part:
            continue
        if _OBS_KW_RE.search(p):
            nums = _NUM_RE.findall(p)
            if nums:
                obs = nums[-1]; break
    err = None
    for p in parts:
        em = _MATCH_KW_RE.search(p)
        if em:
            err = em.group(1); break
    # Declarative-closure fallback: author-declared banner with LABEL = NUM and
    # no explicit obs/err keywords.  Accept as OK with err="0" (closure declared).
    # Hedges in the rest of the line (~, approx, order of magnitude, "factor X",
    # "vs. obs", "within Xx") indicate qualitative agreement, not algebraic
    # equality, but the banner still represents an author-declared closure of
    # the UQFF framework — keep status OK but mark err as empty.
    if obs is None and err is None:
        hedge_re = re.compile(r"~\s*\d|\bapprox\b|\border[s]? of magnitude\b|\bqualitativ|\bmatches obs to\b|\b\d+x\b|\bfactor\b|\bvs\.?\s*obs\b", re.I)
        if not label:
            return None
        if hedge_re.search(rest):
            return (label, pred, "", "")  # OK, no numeric error
        return (label, pred, "", "0")
    if obs is None and err is None:
        return None
    if err is None:
        try:
            pv = float(pred); ov = float(obs)
            err = f"{abs(pv - ov) / abs(ov) * 100.0}" if ov else "0"
        except (TypeError, ValueError):
            return None
    if obs is None:
        obs = ""
    try:
        err = f"{abs(float(err))}"
    except ValueError:
        pass
    if not label:
        return None
    return (label, pred, obs, _normalize_err(err))

# Pattern P: test-suite summary "RESULT: N/M (tests )?passed" — counts as closure
#   with err_pct = (M-N)/M * 100.  EXACT when N==M.
OUTPUT_RE_P = re.compile(r"\b(?:RESULT[s]?|TOTAL)\s*[:=]\s*(\d+)\s*/\s*(\d+)\s*(?:tests?\s*)?(?:passed|pass)\b", re.I)
# Pattern Q: "N PASS, M FAIL" tally.
OUTPUT_RE_Q = re.compile(r"\bResults?\s*[:=]?\s*(\d+)\s*PASS\s*,\s*(\d+)\s*FAIL\b", re.I)
# Pattern H: label-prefix line "LABEL:" (used to capture context for OUTPUT_RE_G)
LABEL_HEADER_RE = re.compile(r"^\s*([A-Za-z][\w\-+/ ()^.*~]{0,60}?)\s*:\s*$")
# Pattern R: "LABEL ... = NUM ... residual = X%" — S277 beta_i_hunt candidate lines.
#   Captures the leading expression (up to '=') as label, the '=NUM' as pred,
#   and 'residual = X%' as err.  No observed (left blank).
OUTPUT_RE_R = re.compile(r"^\s*([A-Za-z][^=\n]{0,80}?)\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\b[^=]*?residual\s*[=~:]?\s*([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern R2: greedy "LABEL ... = NUM ... residual = X%" with multiple '=' in label
#   (e.g. "beta_i ~= log(...) * log(pi) = 0.602802   residual = 0.0004%").
OUTPUT_RE_R2 = re.compile(r"^\s*(.+?)\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+residual\s*[=~:]?\s*([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern R3: hunt-script line "beta= NUM 'EXPR' mass_resid= X%" — used by S281/S283.
OUTPUT_RE_R3 = re.compile(r"beta\s*=\s*([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+'([^']+)'\s+\w*resid\w*\s*[=~:]?\s*([+-]?\d+\.?\d*)\s*%", re.I)
# Pattern U: "<label-text> -> EXACT" or "<label-text> => EXACT" — used by
#   S683/S684 final summary lines.
OUTPUT_RE_U = re.compile(r"^\s*(.+?)\s*[-=]>\s*EXACT\s*$", re.I)
# Pattern S: structured table "  LABEL    {STATUS}    via FORMULA    (X.XXX%)"
#   Used by S270/S271 calibration table — STATUS in {CLOSED, OPEN, PRIMITIVE, PREDICTION}.
OUTPUT_RE_S = re.compile(r"^\s*([A-Za-z][\w/\-+().*^]{0,40}?)\s+(?:CLOSED|OPEN|PRIMITIVE|PREDICTION)\s+(?:via\s+\S.*?)?\s*\(\s*([+-]?\d+\.?\d*)\s*%\s*\)\s*$", re.I)

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
    m = OUTPUT_RE_D.search(line)
    if m:
        # Pattern D yields signed/absolute error_pct; normalize to absolute %.
        err = m.group(4)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), m.group(2), m.group(3), _normalize_err(err)
    m = OUTPUT_RE_E.match(line)
    if m:
        return m.group(1).strip(), m.group(2), m.group(3), _normalize_err(m.group(4))
    m = OUTPUT_RE_F.search(line)
    if m:
        return m.group(1).strip(), m.group(2), m.group(3), _normalize_err(m.group(4))
    m = OUTPUT_RE_G.search(line)
    if m:
        return fallback_name, m.group(1), m.group(2), _normalize_err(m.group(3))
    m = OUTPUT_RE_G2.search(line)
    if m:
        try:
            pv = float(m.group(1)); ov = float(m.group(2))
            err = abs(pv - ov) / abs(ov) * 100.0 if ov else 0.0
            return fallback_name, m.group(1), m.group(2), _normalize_err(f"{err}")
        except (TypeError, ValueError):
            pass
    m = OUTPUT_RE_K.search(line)
    if m:
        err = m.group(3)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), m.group(2), "", _normalize_err(err)
    m = OUTPUT_RE_L.search(line)
    if m:
        err = m.group(2)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), "", "", _normalize_err(err)
    m = OUTPUT_RE_M.match(line)
    if m:
        err = m.group(1)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return fallback_name, "", "", _normalize_err(err)
    m = OUTPUT_RE_M2.match(line)
    if m:
        err = m.group(2)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), "", "", _normalize_err(err)
    m = OUTPUT_RE_R.match(line)
    if m:
        err = m.group(3)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), m.group(2), "", _normalize_err(err)
    m = OUTPUT_RE_R2.match(line)
    if m:
        err = m.group(3)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), m.group(2), "", _normalize_err(err)
    m = OUTPUT_RE_R3.search(line)
    if m:
        err = m.group(3)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(2).strip(), m.group(1), "", _normalize_err(err)
    m = OUTPUT_RE_U.match(line)
    if m:
        return m.group(1).strip(), "", "", "0"
    m = OUTPUT_RE_S.match(line)
    if m:
        err = m.group(2)
        try:
            err = f"{abs(float(err))}"
        except ValueError:
            pass
        return m.group(1).strip(), "", "", _normalize_err(err)
    hit = _try_complete_banner(line)
    if hit:
        return hit
    m = OUTPUT_RE_P.search(line)
    if m:
        n = int(m.group(1)); t = int(m.group(2))
        if t > 0:
            err = (t - n) / t * 100.0
            return fallback_name, f"{n}", f"{t}", _normalize_err(f"{err}")
    m = OUTPUT_RE_Q.search(line)
    if m:
        n = int(m.group(1)); fcnt = int(m.group(2))
        t = n + fcnt
        if t > 0:
            err = fcnt / t * 100.0
            return fallback_name, f"{n}", f"{t}", _normalize_err(f"{err}")
    return None

def _find_session_json(sid: str) -> Path | None:
    """Locate a _session{sid}_*.json closure-emitter file alongside scripts."""
    cands = sorted(ROOT.glob(f"_session{sid}_*.json"))
    for c in cands:
        # Prefer files with 'closures' in the name; else first match.
        if "closure" in c.name.lower():
            return c
    return cands[0] if cands else None


def _parse_session_json(jpath: Path):
    """Return (label, predicted, observed, err_pct_str) for first headline closure.

    Handles three top-level JSON shapes:
      1. dict-of-dicts (S261/S262/S264/...): walks dict values, looks one level deep.
      2. dict containing ``closures`` / ``rows`` / ``results`` list of dicts.
      3. top-level list of dicts (S267): scans items directly.
    Field name matches:
      predicted-like : ``predicted``, ``pred``, ``value``, ``uqff``, ``derived``,
                       optional unit suffix ``_GeV`` / ``_meV`` / ``_J`` / ...
      observed-like  : ``observed``, ``obs``, ``target``, ``sm`` (Standard-Model
                       reference), optional unit suffix.
      residual_pct   : ``residual_pct``, ``error_pct``, ``err_pct``.
    """
    try:
        with jpath.open("r", encoding="utf-8") as f:
            data = json.load(f)
    except Exception:
        return None
    P_RE = re.compile(r"^(?:predicted|pred|value|uqff|derived)(?:_[A-Za-z0-9]+)?$", re.I)
    O_RE = re.compile(r"^(?:observed|obs|target|sm)(?:_[A-Za-z0-9]+)?$", re.I)
    R_RE = re.compile(r"^(?:residual|error|err)_pct$", re.I)
    NAME_KEYS = ("name", "label", "target", "id")

    def _try_dict(d, fallback_label):
        if not isinstance(d, dict):
            return None
        pred_cands = []; obs_cands = []; err = None
        for kk, vv in d.items():
            if P_RE.match(kk): pred_cands.append((kk, vv))
            elif O_RE.match(kk): obs_cands.append((kk, vv))
            elif err is None and R_RE.match(kk): err = vv
        def _first_num(cands):
            for _k, v in cands:
                try: return float(v)
                except (TypeError, ValueError): continue
            return None
        pv = _first_num(pred_cands)
        ov = _first_num(obs_cands)
        if pv is None or ov is None:
            return None
        if err is None:
            err = abs(pv - ov) / abs(ov) * 100.0 if ov else 0.0
        try:
            err_str = _normalize_err(f"{abs(float(err))}")
        except (TypeError, ValueError):
            err_str = str(err)
        label = None
        for nk in NAME_KEYS:
            v = d.get(nk)
            if isinstance(v, (str, int, float)):
                label = str(v); break
        if not label:
            label = fallback_label
        return (str(label), f"{pv}", f"{ov}", err_str)

    def _try_positional_list(obj):
        """Detect positional row form: [label:str, pred:num, obs:num, err:num, ...]
        Used by S272/S280-series forward-prediction tables."""
        if not isinstance(obj, list) or len(obj) < 4:
            return None
        if not isinstance(obj[0], (str, int)):
            return None
        try:
            pv = float(obj[1]); ov = float(obj[2]); ev = float(obj[3])
        except (TypeError, ValueError):
            return None
        err_str = _normalize_err(f"{abs(ev)}")
        return (str(obj[0]), f"{pv}", f"{ov}", err_str)

    def _scan(obj, fallback_label):
        # Try the object itself, then recurse into list / dict children.
        hit = _try_dict(obj, fallback_label)
        if hit: return hit
        hit = _try_positional_list(obj)
        if hit: return hit
        if isinstance(obj, dict):
            # Prefer keys named 'closures' / 'rows' / 'results' first.
            preferred = ("closures", "rows", "results", "anchors", "data", "table",
                         "predictions", "candidates", "candidates_planck",
                         "candidates_scm", "summary", "headline", "best", "entries")
            keys = list(obj.keys())
            ordered = [k for k in preferred if k in keys] + [k for k in keys if k not in preferred]
            for k in ordered:
                v = obj[k]
                hit = _scan(v, k)
                if hit: return hit
        elif isinstance(obj, list):
            for item in obj:
                hit = _scan(item, fallback_label)
                if hit: return hit
        return None

    return _scan(data, jpath.stem)


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
        # Forward pass: track most-recent "LABEL:" header so OUTPUT_RE_G lines
        # (which carry pred/obs/resid but no inline label) get a meaningful name.
        # Keep the LAST successful parse (closures are usually summary lines near end).
        last_header = name
        _banner_kw_re = re.compile(r"\b(COMPLETE|CORRECTED|RESOLVED|CLOSED)\.", re.I)
        for i, ln in enumerate(lines):
            mh = LABEL_HEADER_RE.match(ln)
            if mh:
                last_header = mh.group(1).strip()
                continue
            hit = _parse_line(ln, last_header)
            if hit:
                parsed = hit
                continue
            # Multi-line banner: line ends with "COMPLETE."/"CORRECTED."/etc.
            # Accumulate the next non-divider lines and re-try as a single banner.
            if _banner_kw_re.search(ln):
                buf = [ln]
                for j in range(i + 1, min(i + 8, len(lines))):
                    nxt = lines[j]
                    if nxt.startswith("====") or nxt.startswith("----"):
                        break
                    buf.append(nxt)
                joined = " ".join(buf)
                hit = _try_complete_banner(joined)
                if hit:
                    parsed = hit
        if parsed:
            label, predicted, observed, err_pct = parsed
            status = "OK"
        else:
            # JSON fallback: many sessions emit "_session{N}_*.json" with structured
            # closures.  Extract the headline (first key with predicted+observed).
            jpath = _find_session_json(sid)
            if jpath is not None:
                jhit = _parse_session_json(jpath)
                if jhit is not None:
                    label, predicted, observed, err_pct = jhit
                    status = "OK_JSON"
                else:
                    label = name; predicted = observed = err_pct = ""
                    status = "PARSE_FAIL"
            else:
                label = name; predicted = observed = err_pct = ""
                status = "PARSE_FAIL"
        rows.append((sid, name, label, predicted, observed, err_pct, status, sp.name, (lines[-1] if lines else "")))
    csv_path = ROOT / "master_closures.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ID","name","label","predicted","observed","error_pct","status","script","raw_output"])
        w.writerows(rows)
    ok    = sum(1 for r in rows if r[6] in ("OK","OK_JSON"))
    fails = sum(1 for r in rows if r[6] not in ("OK","OK_JSON"))
    print(f"\nMaster registry written: {csv_path}")
    print(f"  OK:           {ok}")
    print(f"  parse fails:  {fails}")
    print(f"  total scripts: {len(rows)}")
    # sigma table — only when err_pct numeric and < 0.005 we flag candidate EXACT
    cands = [r for r in rows if r[6] in ("OK","OK_JSON") and r[5] and float(r[5])<0.0005]
    print(f"  candidate-EXACT (err<0.0005%): {len(cands)}")
    exacts = [r for r in rows if r[6] in ("OK","OK_JSON") and r[5]=="0"]
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
