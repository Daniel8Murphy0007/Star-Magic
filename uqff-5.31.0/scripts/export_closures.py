#!/usr/bin/env python3
"""Export UQFF closures to LaTeX, BibTeX, CSV, or HTML (Tier-3 I5).

Iterates the calculator's PARADOX_TO_CLOSURE dispatch table, evaluates each
closure, and writes the results in the requested format.

Usage:
    python scripts/export_closures.py --format latex   --out closures.tex
    python scripts/export_closures.py --format bibtex  --out closures.bib
    python scripts/export_closures.py --format csv     --out closures.csv
    python scripts/export_closures.py --format html    --out closures.html
    python scripts/export_closures.py --format markdown --out closures.md

Filter:
    python scripts/export_closures.py --format latex --filter hubble --out hubble.tex
"""
import argparse
import csv
import os
import sys
from typing import Any

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE) if os.path.basename(_HERE) == "scripts" else _HERE
for _p in (os.getcwd(), _REPO_ROOT):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import uqff_pure_calculator as u


def _extract_row(name: str) -> dict[str, Any]:
    """Pull a flat record for one closure: name, uqff_value, target, residual_pct, source."""
    out: dict[str, Any] = {"name": name, "uqff_value": None, "target": None,
                            "residual_pct": None, "primary_source": None,
                            "status_tier": None}
    try:
        result = u.calculate_paradox({"paradox": name})
        v = result.get("value") if isinstance(result, dict) else None
        if isinstance(v, dict):
            out["uqff_value"] = v.get("UQFF_formula_value")
            out["residual_pct"] = v.get("residual_pct")
            out["status_tier"] = v.get("status_tier")
            out["primary_source"] = v.get("primary_source")
            for k in v:
                if k.endswith("_target"):
                    out["target"] = v[k]
                    break
            if out["uqff_value"] is None:
                for fld in ("primary_result", "value"):
                    if fld in v and isinstance(v[fld], (int, float)):
                        out["uqff_value"] = v[fld]
                        break
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            out["uqff_value"] = float(v)
        elif isinstance(v, tuple) and v and isinstance(v[0], (int, float)):
            out["uqff_value"] = float(v[0])
    except Exception as e:
        out["error"] = repr(e)[:60]
    return out


def export_csv(rows: list[dict], outfile: str) -> None:
    cols = ["name", "uqff_value", "target", "residual_pct", "status_tier", "primary_source"]
    with open(outfile, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for r in rows:
            writer.writerow({c: r.get(c, "") for c in cols})


def export_latex(rows: list[dict], outfile: str) -> None:
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("% UQFF Star-Magic closures (auto-generated)\n")
        f.write(f"% {len(rows)} closures exported on auto-run\n\n")
        f.write("\\begin{table}[h!]\n\\centering\n")
        f.write("\\caption{UQFF Star-Magic — selected closure values vs measured targets.}\n")
        f.write("\\label{tab:uqff-closures}\n")
        f.write("\\begin{tabular}{lrrl}\n")
        f.write("\\hline\n")
        f.write("Closure & UQFF value & Target & Residual (\\%) \\\\\n")
        f.write("\\hline\n")
        for r in rows[:200]:  # cap to avoid pathologically large tables
            n = str(r.get("name", "")).replace("_", "\\_")[:50]
            u_val = r.get("uqff_value")
            t_val = r.get("target")
            res = r.get("residual_pct")
            u_str = f"{u_val:.4g}" if isinstance(u_val, (int, float)) else "—"
            t_str = f"{t_val:.4g}" if isinstance(t_val, (int, float)) else "—"
            r_str = f"{res:.4f}" if isinstance(res, (int, float)) else "—"
            f.write(f"\\texttt{{{n}}} & {u_str} & {t_str} & {r_str} \\\\\n")
        if len(rows) > 200:
            f.write(f"\\multicolumn{{4}}{{c}}{{... and {len(rows) - 200} more}} \\\\\n")
        f.write("\\hline\n\\end{tabular}\n\\end{table}\n")


def export_bibtex(rows: list[dict], outfile: str) -> None:
    """One @misc entry per closure, citing the backing PAPER_XXXX."""
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("% UQFF Star-Magic closures (auto-generated BibTeX)\n\n")
        for r in rows:
            name = r.get("name", "")
            paper = r.get("primary_source", "") or "UQFF_general"
            uqff_val = r.get("uqff_value")
            val_str = f"{uqff_val}" if uqff_val is not None else "(structured)"
            f.write(f"@misc{{uqff_{name},\n")
            f.write(f"  title = {{UQFF closure: {name}}},\n")
            f.write(f"  author = {{Murphy, Daniel T.}},\n")
            f.write(f"  year = 2026,\n")
            f.write(f"  note = {{Derivation in {paper}; value = {val_str}}},\n")
            f.write(f"  url = {{https://github.com/Daniel8Murphy0007/Star-Magic/blob/main/whitepapers/{paper}.md}},\n")
            f.write(f"}}\n\n")


def export_html(rows: list[dict], outfile: str) -> None:
    with open(outfile, "w", encoding="utf-8") as f:
        f.write("<!DOCTYPE html><html><head><meta charset='utf-8'>\n")
        f.write("<title>UQFF Closure Catalog</title>\n")
        f.write("<style>body{font-family:system-ui;margin:20px}")
        f.write("table{border-collapse:collapse;width:100%}")
        f.write("th,td{padding:6px 10px;border-bottom:1px solid #ccc;text-align:left}")
        f.write("th{background:#ddf;color:#003}")
        f.write("td.num{text-align:right;font-variant-numeric:tabular-nums}")
        f.write("</style></head><body>\n")
        f.write(f"<h1>UQFF Star-Magic Closure Catalog</h1>\n")
        f.write(f"<p>{len(rows)} closures auto-exported.</p>\n")
        f.write("<table><thead><tr><th>Closure</th><th>UQFF value</th><th>Target</th>")
        f.write("<th>Residual %</th><th>Status</th><th>Paper</th></tr></thead><tbody>\n")
        for r in rows:
            n = r.get("name", "")
            u_val = r.get("uqff_value")
            t_val = r.get("target")
            res = r.get("residual_pct")
            status = r.get("status_tier", "") or ""
            paper = r.get("primary_source", "") or ""
            u_str = f"{u_val:.6g}" if isinstance(u_val, (int, float)) else "—"
            t_str = f"{t_val:.6g}" if isinstance(t_val, (int, float)) else "—"
            r_str = f"{res:.4f}" if isinstance(res, (int, float)) else "—"
            f.write(f"<tr><td><code>{n}</code></td>"
                    f"<td class='num'>{u_str}</td><td class='num'>{t_str}</td>"
                    f"<td class='num'>{r_str}</td><td>{status}</td><td>{paper}</td></tr>\n")
        f.write("</tbody></table></body></html>\n")


def export_markdown(rows: list[dict], outfile: str) -> None:
    with open(outfile, "w", encoding="utf-8") as f:
        f.write(f"# UQFF Star-Magic Closure Catalog\n\n")
        f.write(f"{len(rows)} closures auto-exported.\n\n")
        f.write("| Closure | UQFF value | Target | Residual (%) | Status | Paper |\n")
        f.write("|---|---:|---:|---:|---|---|\n")
        for r in rows:
            n = r.get("name", "")
            u_val = r.get("uqff_value")
            t_val = r.get("target")
            res = r.get("residual_pct")
            status = r.get("status_tier", "") or ""
            paper = r.get("primary_source", "") or ""
            u_str = f"{u_val:.6g}" if isinstance(u_val, (int, float)) else "—"
            t_str = f"{t_val:.6g}" if isinstance(t_val, (int, float)) else "—"
            r_str = f"{res:.4f}" if isinstance(res, (int, float)) else "—"
            f.write(f"| `{n}` | {u_str} | {t_str} | {r_str} | {status} | {paper} |\n")


_FORMATS = {
    "csv": export_csv,
    "latex": export_latex,
    "bibtex": export_bibtex,
    "html": export_html,
    "markdown": export_markdown,
}


def main() -> int:
    parser = argparse.ArgumentParser(description="Export UQFF closures to multiple formats")
    parser.add_argument("--format", choices=list(_FORMATS), required=True,
                        help="output format")
    parser.add_argument("--out", required=True, help="output file path")
    parser.add_argument("--filter", help="substring filter on closure names (optional)")
    parser.add_argument("--limit", type=int, default=0,
                        help="cap number of exported rows (0 = no cap)")
    parser.add_argument("--only-schema", action="store_true",
                        help="export only closures with full schema (target, residual, paper)")
    args = parser.parse_args()

    keys = sorted(u.PARADOX_TO_CLOSURE.keys())
    if args.filter:
        flt = args.filter.lower()
        keys = [k for k in keys if flt in k.lower()]

    print(f"Evaluating {len(keys)} closures...")
    rows = []
    for k in keys:
        r = _extract_row(k)
        if args.only_schema and (r["target"] is None or r["primary_source"] is None):
            continue
        rows.append(r)

    if args.limit:
        rows = rows[: args.limit]

    print(f"Writing {len(rows)} rows to {args.out} as {args.format}...")
    _FORMATS[args.format](rows, args.out)
    size = os.path.getsize(args.out)
    print(f"  OK: {size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
