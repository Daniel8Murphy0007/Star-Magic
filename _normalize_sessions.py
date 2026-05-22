"""Read-only normalizer: extract primitive declarations from every _session*.py.

Output: SESSIONS_PRIMITIVES_NORMALIZED.csv  (long-form, one row per declaration)
Columns: file,line,primitive,raw_value,literal_form,scope
"""
from __future__ import annotations
import csv
import re
from pathlib import Path

ROOT = Path(__file__).parent

# Primitive name -> regex of accepted LHS aliases.
# Match either at top of line (module scope) or single-indent (class scope) or
# deeper (function scope) — we record scope by indent depth.
PRIMITIVES = {
    "D_PHYS":      r"D_PHYS|D_phys",
    "D_BSFG":      r"D_BSFG|D_bsfg",
    "D_CRIT":      r"D_CRIT|D_crit",
    "N_CH":        r"N_CH|N_ch",
    "ORD_SO5":     r"ORD_SO5|SO5_ORDER|SO5_ord|dim_SO5",
    "PHI_RES":     r"PHI_RES|Phi_res|phi_res",
    "K_MEX":       r"K_MEX|K_Mex|KMEX",
    "F_TRZ":       r"F_TRZ|f_TRZ|fTRZ|F_trz",
    "SSQ":         r"SSQ|SSq|\[SSq\]",
    "KAPPA":       r"KAPPA|KAPPA_PER_DAY",
    "H_SCM":       r"H_SCM|H_SCm|HSCm",
    "V_SCM":       r"V_SCM|v_SCm|V_SCm",
    "V_UA":        r"V_UA|v_UA",
    "BETA_I":      r"BETA_I|beta_i|BETA_i|Beta_i",
    "RHO_VAC_SCM": r"RHO_VAC_SCM|RHO_SCM|rho_vac_SCm|rho_SCm",
    "RHO_VAC_UA":  r"RHO_VAC_UA|RHO_UA|rho_vac_UA|rho_UA",
}

# Combined LHS pattern with named groups
LHS_PARTS = [f"(?P<{k}>{v})" for k, v in PRIMITIVES.items()]
LHS_RE = re.compile(
    r"^(?P<indent>\s*)(?:" + "|".join(LHS_PARTS) + r")\s*=\s*(?P<rhs>[^#\n]+?)\s*(?:#.*)?$"
)

def scope_of(indent: str) -> str:
    n = len(indent.expandtabs(4))
    if n == 0:
        return "module"
    if n == 4:
        return "class_or_func"
    return f"nested_{n}"

def normalize_value(raw: str) -> str:
    """Best-effort numeric normalization; falls back to raw expression."""
    expr = raw.strip().rstrip(",;")
    # Common patterns we accept as-is (Fraction, Rational, sympy)
    if any(tok in expr for tok in ("Fraction(", "Rational(", "sympify(", "sqrt(", "pi", "math.")):
        return expr
    try:
        # Only allow simple arithmetic / floats / ints
        if re.fullmatch(r"[\d\.\-+eE_\*/\s\(\)]+", expr):
            val = eval(expr, {"__builtins__": {}}, {})
            if isinstance(val, (int, float)):
                return repr(val)
    except Exception:
        pass
    return expr

def main() -> None:
    files = sorted(ROOT.glob("_session*.py"))
    out_path = ROOT / "SESSIONS_PRIMITIVES_NORMALIZED.csv"
    rows = 0
    file_hits: dict[str, int] = {}
    with out_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "line", "primitive", "raw_value", "literal_form", "scope"])
        for fp in files:
            try:
                text = fp.read_text(encoding="utf-8", errors="replace")
            except Exception as exc:
                w.writerow([fp.name, 0, "READ_ERROR", str(exc), "", ""])
                continue
            n_in_file = 0
            for lineno, line in enumerate(text.splitlines(), 1):
                m = LHS_RE.match(line)
                if not m:
                    continue
                # Identify which primitive matched
                matched = next(
                    (k for k in PRIMITIVES if m.group(k) is not None), None
                )
                if not matched:
                    continue
                rhs = m.group("rhs")
                w.writerow([
                    fp.name,
                    lineno,
                    matched,
                    rhs.strip(),
                    normalize_value(rhs),
                    scope_of(m.group("indent")),
                ])
                rows += 1
                n_in_file += 1
            if n_in_file:
                file_hits[fp.name] = n_in_file
    print(f"Files scanned:       {len(files)}")
    print(f"Files with hits:     {len(file_hits)}")
    print(f"Total declarations:  {rows}")
    print(f"Output:              {out_path.name}")

if __name__ == "__main__":
    main()
