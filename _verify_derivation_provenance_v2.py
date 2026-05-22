"""
_verify_derivation_provenance_v2.py
Rigorous per-script forensic classification.

Distinguishes:
  STRUCTURAL_IDENTITY  - script computes finite-math facts: set sizes, file
                         line counts, modular arithmetic, Laplacian eigenvalues,
                         polylogarithm sums, combinatorial identities.
                         Both 'pred' and 'obs' are derived from the same
                         math; no CODATA involved. (e.g. S201, S202)
  PRIMITIVE_DERIVATION - 'pred' is a single closed expression in the 11 locked
                         primitives (and pure-math constants pi/e/sqrt); 'obs'
                         is loaded from CODATA/PDG for comparison only.
  SHORTLIST_PICK       - script tries 2..30 hand-written candidate expressions
                         and picks the best-fitting one. Forward in spirit
                         but allows curve-picking. (e.g. S272)
  ATOM_SEARCH_FIT      - script defines an ATOMS dictionary and runs an
                         exhaustive search2/search3/search4 brute-force over
                         it. 'pred' = best fit. ATOMS may include CODATA
                         constants. (e.g. S747, S785)
  OBS_TAUTOLOGY        - 'pred' assignment references '_obs' variable.
  NO_CLOSURE_OUTPUT    - script does not emit any predicted/observed line.

Output: provenance_audit_v2.csv with detailed evidence per script.
"""
from __future__ import annotations
import re, csv
from pathlib import Path

ROOT = Path(__file__).parent

# True brute-force search definitions (ATOMS-based, multi-arg combinatorial).
# Must look like:  def searchN(target, ...): ... for ... in ATOMS / items / VALS ...
SEARCH_DEF_RE = re.compile(
    r"^\s*def\s+search[234](?:_[a-z]+)?\s*\(\s*target",
    re.MULTILINE,
)
# Same call site
SEARCH_CALL_RE = re.compile(r"\bsearch[234](?:_[a-z]+)?\s*\(\s*[A-Za-z_]")

# Atom dictionary marker
ATOMS_DICT_RE = re.compile(r"^\s*ATOMS\s*=\s*\{", re.MULTILINE)

# Candidate-list pattern (S272-style)
CANDIDATES_DICT_RE = re.compile(r"^\s*candidates\s*=\s*\{", re.MULTILINE)
BEST_PICK_RE = re.compile(r"\bmin\s*\(\s*candidates\.items\(\)")

# Closure output: any of these printed forms
CLOSURE_PRINT_PATTERNS = [
    re.compile(r"predicted\s*=", re.I),
    re.compile(r"pred\s*=", re.I),
    re.compile(r"->.*EXACT", re.I),
    re.compile(r"vs\s+[-\d.e+]+\s*->", re.I),
]

# CODATA / external observed value names
CODATA_NAME_RE = re.compile(
    r"\b(?:alpha_em|alpha_inv|mpme|mp_over_me|"
    r"m_p(?:_obs)?\s*=\s*[\d.e\-+]|"
    r"m_e(?:_obs)?\s*=\s*[\d.e\-+]|"
    r"hbar\s*=\s*[\d.e\-+]|"
    r"G_SI|G\s*=\s*6\.6[47]|"
    r"r_p_obs|Lambda_obs|eta_obs|Y_p_obs|D_H_obs|H0_obs|H_0_obs|"
    r"RHO_LAMBDA_OBS|sin2_thW|J_per_(?:GeV|MeV|eV))"
)
OBS_ASSIGN_RE = re.compile(r"^\s*\w*_obs\s*=\s*[-+]?[\d.eE+\-]+", re.MULTILINE)

# Structural/identity hints
STRUCTURAL_HINTS = [
    re.compile(r"definitional", re.I),
    re.compile(r"set cardinal", re.I),
    re.compile(r"modular", re.I),
    re.compile(r"mod\s+\d", re.I),
    re.compile(r"line count", re.I),
    re.compile(r"Laplacian|eigenvalue|polylog|Li_\d+", re.I),
    re.compile(r"closure footprint", re.I),
    re.compile(r"NULL extraction", re.I),
    re.compile(r"number-theoretic", re.I),
]


def classify(path: Path):
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    # Strip comments/docstrings for code-only checks
    # Quick & dirty: drop lines starting with #
    code_text = "\n".join(ln for ln in lines if not ln.lstrip().startswith("#"))

    has_atoms_dict = bool(ATOMS_DICT_RE.search(code_text))
    has_search_def = bool(SEARCH_DEF_RE.search(code_text))
    has_search_call = bool(SEARCH_CALL_RE.search(code_text))
    has_candidates_dict = bool(CANDIDATES_DICT_RE.search(code_text))
    has_best_pick = bool(BEST_PICK_RE.search(code_text))

    has_closure_output = any(p.search(text) for p in CLOSURE_PRINT_PATTERNS)

    codata_assigns = list(CODATA_NAME_RE.finditer(code_text))
    obs_assigns = list(OBS_ASSIGN_RE.finditer(code_text))

    structural_score = sum(1 for p in STRUCTURAL_HINTS if p.search(text))

    # Check for pred = ... _obs ... in actual code
    obs_in_pred = False
    for ln in lines:
        s = ln.strip()
        if s.startswith("#") or "=" not in s:
            continue
        if re.match(r"(pred(?:icted)?|p)\s*=", s) and re.search(r"_obs\b", s):
            obs_in_pred = True
            break

    # Decision tree
    if not has_closure_output:
        cls = "NO_CLOSURE_OUTPUT"
    elif obs_in_pred:
        cls = "OBS_TAUTOLOGY"
    elif has_atoms_dict and (has_search_def or has_search_call):
        cls = "ATOM_SEARCH_FIT"
    elif has_candidates_dict and has_best_pick:
        cls = "SHORTLIST_PICK"
    elif structural_score >= 1 and not codata_assigns:
        cls = "STRUCTURAL_IDENTITY"
    elif codata_assigns or obs_assigns:
        cls = "PRIMITIVE_DERIVATION"  # provisional: CODATA present but used for comparison
    else:
        cls = "STRUCTURAL_IDENTITY"  # no CODATA, has closure output

    return {
        "classification": cls,
        "has_atoms_dict": has_atoms_dict,
        "has_search_def": has_search_def,
        "has_search_call": has_search_call,
        "has_candidates_dict": has_candidates_dict,
        "has_best_pick": has_best_pick,
        "has_closure_output": has_closure_output,
        "codata_count": len(codata_assigns),
        "obs_assign_count": len(obs_assigns),
        "structural_score": structural_score,
        "obs_in_pred": obs_in_pred,
    }


def main():
    sessions = sorted(ROOT.glob("_session*.py"))
    rows = []
    counts = {}
    for sp in sessions:
        d = classify(sp)
        counts[d["classification"]] = counts.get(d["classification"], 0) + 1
        rows.append({"script": sp.name, **d})

    out = ROOT / "provenance_audit_v2.csv"
    fields = ["script", "classification", "has_atoms_dict", "has_search_def",
              "has_search_call", "has_candidates_dict", "has_best_pick",
              "has_closure_output", "codata_count", "obs_assign_count",
              "structural_score", "obs_in_pred"]
    with out.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"Total session scripts: {len(sessions)}")
    print(f"Output: {out.name}\n")
    print("Classification breakdown:")
    for k in sorted(counts, key=lambda x: -counts[x]):
        print(f"  {k:24s} {counts[k]:4d}")

    print("\nSample of each class (first 5):")
    by_class = {}
    for r in rows:
        by_class.setdefault(r["classification"], []).append(r["script"])
    for cls, names in by_class.items():
        print(f"\n[{cls}]  ({len(names)} scripts)")
        for n in names[:5]:
            print(f"  {n}")


if __name__ == "__main__":
    main()
