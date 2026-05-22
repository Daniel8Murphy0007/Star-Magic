"""
_verify_derivation_provenance.py
Read-only forensic survey of every _session*.py.

Classifies each script into one of:
  PRIMITIVE_ONLY  - predicted values come only from the 11 locked primitives
                    (and integer/rational combinations of them). No CODATA,
                    no observed value, no search.
  PRIMITIVE_PLUS  - uses primitives plus pi/e/sqrt or other dimensionless
                    pure-math constants, but no CODATA input.
  CODATA_INPUT    - imports a measured constant (alpha_em, mpme, hbar, G, ...)
                    and uses it INSIDE the predicted-value computation
                    (not just for comparison).
  SEARCH_FIT      - defines / calls a search2/search3/search4 brute-force
                    that scans atom dictionaries for a value matching an
                    observed target. predicted = best fit, not derivation.
  OBS_TAUTOLOGY   - assigns <name>_obs = NUMBER and the predicted value
                    references that obs variable (or derives from it).
  COMPARISON_ONLY - imports CODATA/PDG values but only uses them after the
                    prediction is computed (legitimate residual reporting).
  UNCLASSIFIED    - none of the above patterns matched cleanly.

Writes provenance_audit.csv with columns:
  script, classification, evidence_lines, predicted_print_count
"""
from __future__ import annotations
import re, csv, ast
from pathlib import Path

ROOT = Path(__file__).parent

# 11 locked primitives (and a few derived aliases that are pure rationals over them)
PRIMITIVE_NAMES = {
    "F_TRZ", "Phi_res", "SSq", "K_Mex", "D_phys", "D_BSFG", "D_crit",
    "N_ch", "SO5", "A_5", "beta_i",
    # pure-math constants we count as "PRIMITIVE_PLUS" not "CODATA_INPUT"
}
PURE_MATH_NAMES = {"pi", "e", "PI", "E"}

# Names that, if assigned a numeric literal in the script, count as CODATA input
CODATA_NAME_PATTERNS = [
    re.compile(r"\balpha_em\b"),
    re.compile(r"\bmpme\b"),
    re.compile(r"\binv_137\b"),
    re.compile(r"\bm_p_obs\b|\bm_p\b\s*=\s*[0-9]"),
    re.compile(r"\bm_e_obs\b|\bm_e\b\s*=\s*[0-9]"),
    re.compile(r"\bhbar\b\s*=\s*[0-9]"),
    re.compile(r"\bG_SI\b|\bG\b\s*=\s*6\.6743"),
    re.compile(r"\br_p_obs"),
    re.compile(r"\bLambda_obs"),
    re.compile(r"\beta_obs|\bY_p_obs|\bD_H_obs"),
    re.compile(r"\bH0_obs|\bH_0_obs"),
    re.compile(r"\bRHO_LAMBDA_OBS"),
    re.compile(r"\bsin.*theta.*_obs"),
    re.compile(r"_obs\s*=\s*[0-9]"),
]

CODATA_COMMENT_PATTERNS = [
    re.compile(r"CODATA", re.I),
    re.compile(r"NuFIT|NuFit", re.I),
    re.compile(r"Planck 2018", re.I),
    re.compile(r"\bPDG\b"),
]

SEARCH_DEF_RE = re.compile(r"def\s+search[234abc]?\s*\(", re.I)
SEARCH_CALL_RE = re.compile(r"\bsearch[234abc]?\s*\(", re.I)

# Pattern: predicted=X observed=Y, with variable names
PRED_PRINT_RE = re.compile(
    r"f?[\"']\s*[^\"']*?predicted\s*=\s*\{?\s*(\w[\w\.\[\]]*?)\s*[:\}\s].*?observed",
    re.I,
)
# Looser: any "predicted=" line that references an obs var
PRED_USES_OBS_RE = re.compile(r"predicted\s*=\s*[^,]*?_obs\b", re.I)


def classify(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    evidence: list[str] = []

    has_search_def = False
    has_search_call = False
    for i, ln in enumerate(lines, 1):
        if SEARCH_DEF_RE.search(ln):
            has_search_def = True
            evidence.append(f"L{i}: search-def: {ln.strip()[:90]}")
        elif SEARCH_CALL_RE.search(ln) and "def " not in ln:
            # crude but enough to confirm a call site
            if not has_search_call:
                evidence.append(f"L{i}: search-call: {ln.strip()[:90]}")
            has_search_call = True

    # CODATA name assignments (not just comments)
    codata_assignment_lines = []
    for i, ln in enumerate(lines, 1):
        stripped = ln.strip()
        if stripped.startswith("#"):
            continue
        if "=" not in stripped:
            continue
        for pat in CODATA_NAME_PATTERNS:
            if pat.search(stripped):
                codata_assignment_lines.append((i, stripped))
                break

    # Was any CODATA assignment used in a 'predicted=' computation?
    # Heuristic: look at f-string print lines containing 'predicted=' and see
    # if any CODATA name appears in the expression.
    codata_used_in_prediction = False
    obs_used_in_prediction = False
    pred_print_count = 0
    # Scan multi-line: simple per-line scan
    for i, ln in enumerate(lines, 1):
        if "predicted" not in ln.lower():
            continue
        if "=" not in ln:
            continue
        pred_print_count += 1
        # Look at the *expression* fed to the predicted= slot.
        # Common pattern: f"...predicted={expr:.4e} observed={obs_var:.4e}..."
        m = re.search(r"predicted\s*=\s*\{([^:}]+)[:\}]", ln)
        if m:
            expr = m.group(1)
            # Does the expression contain a CODATA-tagged name?
            for pat in CODATA_NAME_PATTERNS:
                if pat.search(expr):
                    codata_used_in_prediction = True
                    if len(evidence) < 12:
                        evidence.append(f"L{i}: pred-uses-CODATA: {ln.strip()[:100]}")
                    break
            if re.search(r"_obs\b", expr):
                obs_used_in_prediction = True
                if len(evidence) < 12:
                    evidence.append(f"L{i}: pred-uses-_obs: {ln.strip()[:100]}")

    # Also check direct assignment: pred = ... _obs ...
    for i, ln in enumerate(lines, 1):
        s = ln.strip()
        if s.startswith("#"):
            continue
        if re.match(r"\s*(pred|prediction|p|y_pred)\b\s*=", s) and re.search(r"_obs\b", s):
            obs_used_in_prediction = True
            if len(evidence) < 14:
                evidence.append(f"L{i}: pred-assigned-from-_obs: {s[:100]}")

    has_codata_assignment = len(codata_assignment_lines) > 0
    if has_codata_assignment and not codata_used_in_prediction and not obs_used_in_prediction:
        # Could still be comparison-only.
        if codata_assignment_lines:
            i, s = codata_assignment_lines[0]
            evidence.append(f"L{i}: codata-assign-only: {s[:100]}")

    # Decision tree (priority order)
    if has_search_def or has_search_call:
        cls = "SEARCH_FIT"
    elif obs_used_in_prediction:
        cls = "OBS_TAUTOLOGY"
    elif codata_used_in_prediction:
        cls = "CODATA_INPUT"
    elif has_codata_assignment:
        cls = "COMPARISON_ONLY"
    else:
        # No CODATA, no search, no obs in prediction.
        # Check for ANY numeric literal in a pred=... line that isn't from primitives.
        # For now mark as PRIMITIVE_ONLY (best case) — caller can refine.
        cls = "PRIMITIVE_ONLY" if pred_print_count > 0 else "UNCLASSIFIED"

    return cls, evidence[:15], pred_print_count


def main():
    sessions = sorted(ROOT.glob("_session*.py"))
    rows = []
    counts = {}
    for sp in sessions:
        cls, ev, npred = classify(sp)
        counts[cls] = counts.get(cls, 0) + 1
        rows.append((sp.name, cls, npred, " | ".join(ev)))

    out = ROOT / "provenance_audit.csv"
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["script", "classification", "predicted_print_count", "evidence"])
        for r in rows:
            w.writerow(r)

    print(f"Total session scripts scanned: {len(sessions)}")
    print(f"Output: {out.name}")
    print()
    print("Classification breakdown:")
    for k in sorted(counts, key=lambda x: -counts[x]):
        print(f"  {k:18s} {counts[k]:4d}")
    print()
    print("Sample of each class (first 3):")
    seen = {}
    for name, cls, npred, ev in rows:
        seen.setdefault(cls, []).append((name, npred, ev))
    for cls, items in seen.items():
        print(f"\n[{cls}]")
        for name, npred, ev in items[:3]:
            print(f"  {name}  (pred-prints={npred})")
            if ev:
                for line in ev.split(" | ")[:3]:
                    print(f"      {line}")


if __name__ == "__main__":
    main()
