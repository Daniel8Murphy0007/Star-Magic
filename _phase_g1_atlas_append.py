"""
_phase_g1_atlas_append.py
Phase G1 — Append the assimilation overlay section to CLOSURE_ATLAS.md.

Builds Section 12 (Assimilation overlay) from assimilation_dispatch.py, including:
  - Per-domain rollup
  - Per-observable table (114 rows)
  - geometry_used / numeric_system / assimilation_status tagging columns
  - cross-refs to OVERDETERMINATION_MAP.csv and the Phase F public surfaces

Idempotent: detects if Section 12 already exists and replaces in place.
"""
import sys, re
from pathlib import Path
sys.path.insert(0, '.')
import assimilation_dispatch as ad

ROOT = Path(__file__).resolve().parent
ATLAS = ROOT / "CLOSURE_ATLAS.md"


def build_section():
    lines = []
    lines.append("---")
    lines.append("")
    lines.append("## 12. Assimilation overlay (Phase E + F)")
    lines.append("")
    lines.append("**Built:** 2026-06-29 (Round 670). Updated by `_phase_g1_atlas_append.py`.")
    lines.append("**Source of truth:** `assimilation_dispatch.py` (114 curated observables).")
    lines.append("**Public API:** 8 Phase F surfaces in `uqff_pure_calculator.py`:")
    lines.append("  `calculate_qcalcgeom_compute_{FUBi,FUBii,F_U,emergent_mass}`,")
    lines.append("  `calculate_qcalcgeom_solve_habitable_zone`,")
    lines.append("  `calculate_3numeric_decomposition`, `calculate_geometry_decomposition`,")
    lines.append("  `calculate_overdetermination`, plus the `qcalcgeom_solve` dispatch key in")
    lines.append("  `calculate_analytic_closures`.")
    lines.append("**Solver bus:** `qcalcgeom_solver.solve(observable, geometry, numeric, decompose)`.")
    lines.append("")

    # Per-domain rollup
    by_dom = {}
    for name, rec in ad.DISPATCH.items():
        by_dom.setdefault(rec['domain'], []).append((name, rec))

    lines.append("### 12.1 Per-domain rollup (114 observables, 10 domains)")
    lines.append("")
    lines.append("| Domain | Observables | Owner geometries | Worst residual |")
    lines.append("|---|---:|---|---:|")
    for dom in sorted(by_dom):
        obs = by_dom[dom]
        geoms = sorted({r['owner_geometry'] for _, r in obs})
        worst = max((float(r.get('residual_pct') or 0.0) for _, r in obs), default=0.0)
        lines.append(f"| {dom} | {len(obs)} | {', '.join(geoms)} | {worst:.4f}% |")
    lines.append("")

    # Per-observable table
    lines.append("### 12.2 Full 114-observable inventory")
    lines.append("")
    lines.append("Every row dispatchable via `qcalcgeom_solver.solve(observable)` or")
    lines.append("`calculate_analytic_closures({'qcalcgeom_solve': {'observable': name}})`.")
    lines.append("")
    lines.append("| Observable | Domain | Owner | Residual | Source |")
    lines.append("|---|---|---|---:|---|")
    for name in sorted(ad.DISPATCH.keys()):
        rec = ad.DISPATCH[name]
        rp = rec.get('residual_pct')
        rp_str = f"{rp:.4f}%" if isinstance(rp, (int, float)) else "n/a"
        src = rec.get('primary_source', '')
        lines.append(f"| `{name}` | {rec['domain']} | {rec['owner_geometry']} | {rp_str} | {src} |")
    lines.append("")

    # OVERDETERMINATION_MAP cross-ref
    lines.append("### 12.3 OVERDETERMINATION_MAP cross-references")
    lines.append("")
    lines.append("- `OVERDETERMINATION_MAP.csv` — long format, 1,368 rows = 114 obs × 4 geom × 3 numeric.")
    lines.append("- `OVERDETERMINATION_WIDE.csv` — wide format, 114 rows × 18 cols.")
    lines.append("- `OVERDETERMINATION_MAP.md` — peer-review summary with per-domain rollup, TENSION block,")
    lines.append("  and the Round 669 BAO multi-path closure section.")
    lines.append("- Round 669 closed the only OPEN_QUESTION (BAO) with dual closures; the framework now")
    lines.append("  ships with **0 TENSION cells**.")
    lines.append("")

    # Cross-refs
    lines.append("### 12.4 Discovery cheat sheet")
    lines.append("")
    lines.append("```")
    lines.append("# List all 114 observables programmatically:")
    lines.append("import assimilation_dispatch as ad")
    lines.append("for name in sorted(ad.DISPATCH):")
    lines.append("    print(name, ad.DISPATCH[name]['domain'])")
    lines.append("")
    lines.append("# Call any observable through the calculator:")
    lines.append("import uqff_pure_calculator as u")
    lines.append("u.calculate_analytic_closures({'qcalcgeom_solve': {'observable': 'alpha_inverse'}})")
    lines.append("# -> {'value': 137.0}")
    lines.append("")
    lines.append("# Decomposed view with full provenance (8-field result dict):")
    lines.append("u.calculate_analytic_closures(")
    lines.append("    {'qcalcgeom_solve': {'observable': 'LCDM_BAO_rd_H0_over_c_primary', 'decompose': True}})")
    lines.append("```")
    lines.append("")

    return "\n".join(lines)


def main():
    text = ATLAS.read_text(encoding='utf-8')
    section_anchor = "## 12. Assimilation overlay (Phase E + F)"
    new_section = build_section()

    if section_anchor in text:
        # Idempotent: replace existing section (between anchor and next '## ' or EOF)
        start = text.find("---\n\n" + section_anchor)
        if start < 0:
            start = text.find(section_anchor)
            start = text.rfind("---", 0, start)
        # Find end: next '## N.' header or EOF
        rest = text[start + 10:]
        next_header = re.search(r"^---\n\n## \d+\.", rest, re.MULTILINE)
        end = (start + 10 + next_header.start()) if next_header else len(text)
        new_text = text[:start] + new_section + ("\n" + text[end:] if end < len(text) else "\n")
        action = "replaced"
    else:
        new_text = text.rstrip() + "\n\n" + new_section + "\n"
        action = "appended"

    ATLAS.write_text(new_text, encoding='utf-8')
    print(f"CLOSURE_ATLAS.md {action}: section 12 ({len(new_section)} chars)")
    print(f"File now: {len(new_text)} bytes, {new_text.count(chr(10))} lines")


if __name__ == "__main__":
    main()
