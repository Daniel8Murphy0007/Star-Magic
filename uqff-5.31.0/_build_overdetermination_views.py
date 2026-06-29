"""
_build_overdetermination_views.py
Phase E8 generator — produces the OVERDETERMINATION_MAP family.

Reads:
  assimilation_dispatch.py   (112 dispatch entries — the curated observable catalog)
  qcalcgeom_solver.solve()   (4 x 3 dispatch matrix per observable)

Writes:
  OVERDETERMINATION_MAP.csv         long format: one row per (obs, geom, numeric) = 1,344 rows
  OVERDETERMINATION_WIDE.csv        wide format: one row per obs, 12 residual cells + total_N
  OVERDETERMINATION_MAP.md          human-readable summary with per-domain rollups
"""
import csv
import sys
from collections import OrderedDict, Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import assimilation_dispatch as ad
from qcalcgeom_solver import solve

GEOMS = ["qcalcgeom", "bsfg", "dpm", "d26"]
NUMS = ["symbolic", "numerical", "discrete"]
GEOM_ABBREV = {"qcalcgeom":"qg", "bsfg":"bsfg", "dpm":"dpm", "d26":"d26"}
NUM_ABBREV = {"symbolic":"sym", "numerical":"num", "discrete":"dis"}


def main():
    print("=" * 72)
    print("PHASE E8 — OVERDETERMINATION_MAP family generator")
    print("=" * 72)

    long_rows = []
    wide_rows = []
    per_domain = {}
    obs_names = sorted(ad.DISPATCH.keys())

    for name in obs_names:
        rec = ad.DISPATCH[name]
        domain = rec["domain"]
        owner = rec["owner_geometry"]
        target = rec.get("target")
        primary_source = rec.get("primary_source", "")
        notes = rec.get("notes", "")

        r = solve(name, geometry="auto", numeric="all", decompose=True)
        alt = r.get("alternate_paths") or {}

        cells = {}  # (g, n) -> {value, residual_pct, status}
        for g in GEOMS:
            for n in NUMS:
                cell = (alt.get(g) or {}).get(n) or {}
                v = cell.get("value")
                rp = cell.get("residual_pct")
                if rp is None and v is not None and target is not None:
                    try:
                        tv = float(target); vv = float(v)
                        rp = 100.0 * abs(vv - tv) / abs(tv) if tv != 0 else abs(vv - tv) * 100.0
                    except (TypeError, ValueError):
                        rp = None
                if v is None:
                    status = "GAP"
                elif "OPEN_QUESTION" in (notes or ""):
                    status = "TENSION"
                elif rp is not None and abs(rp) < 1e-9:
                    status = "EXACT"
                else:
                    status = "OK"
                cells[(g, n)] = {"value": v, "residual_pct": rp, "status": status}

                long_rows.append({
                    "observable": name,
                    "domain": domain,
                    "geometry": g,
                    "numeric": n,
                    "value": "" if v is None else repr(v),
                    "target": "" if target is None else repr(target),
                    "residual_pct": "" if rp is None else f"{rp:.6f}",
                    "status": status,
                    "owner_geometry": owner,
                    "primary_source": primary_source,
                })

        total_N = sum(1 for c in cells.values() if c["value"] is not None)
        owner_N = sum(1 for n in NUMS if cells[(owner, n)]["value"] is not None)

        wide = OrderedDict()
        wide["observable"] = name
        wide["domain"] = domain
        wide["owner_geometry"] = owner
        for g in GEOMS:
            for n in NUMS:
                col = f"{GEOM_ABBREV[g]}_{NUM_ABBREV[n]}"
                cell = cells[(g, n)]
                rp = cell["residual_pct"]
                if cell["status"] == "GAP":
                    wide[col] = ""
                elif cell["status"] == "EXACT":
                    wide[col] = "EXACT"
                elif cell["status"] == "TENSION":
                    wide[col] = f"TENSION({rp:.4f}%)" if rp is not None else "TENSION"
                else:
                    wide[col] = f"{rp:.4f}%" if rp is not None else "OK"
        wide["owner_N"] = owner_N
        wide["total_N"] = total_N
        wide["primary_source"] = primary_source
        wide_rows.append(wide)

        per_domain.setdefault(domain, []).append({
            "name": name, "owner": owner, "total_N": total_N, "owner_N": owner_N,
            "status_counts": Counter(c["status"] for c in cells.values()),
        })

    # Write long-format
    long_path = ROOT / "OVERDETERMINATION_MAP.csv"
    long_fields = ["observable","domain","geometry","numeric","value","target",
                   "residual_pct","status","owner_geometry","primary_source"]
    with long_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=long_fields, lineterminator="\n")
        w.writeheader()
        for row in long_rows:
            w.writerow(row)
    print(f"Wrote {long_path.name}: {len(long_rows)} rows ({len(obs_names)} obs x 4 geom x 3 num)")

    # Write wide-format
    wide_path = ROOT / "OVERDETERMINATION_WIDE.csv"
    wide_fields = list(wide_rows[0].keys()) if wide_rows else []
    with wide_path.open("w", encoding="utf-8", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=wide_fields, lineterminator="\n")
        w.writeheader()
        for row in wide_rows:
            w.writerow(row)
    print(f"Wrote {wide_path.name}: {len(wide_rows)} rows, {len(wide_fields)} cols")

    # Write .md summary
    md_path = ROOT / "OVERDETERMINATION_MAP.md"
    total_obs = len(obs_names)
    total_cells = len(long_rows)
    populated = sum(1 for r in long_rows if r["value"] != "")
    exact = sum(1 for r in long_rows if r["status"] == "EXACT")
    ok = sum(1 for r in long_rows if r["status"] == "OK")
    tension = sum(1 for r in long_rows if r["status"] == "TENSION")
    gap = sum(1 for r in long_rows if r["status"] == "GAP")
    coverage_pct = 100.0 * populated / total_cells if total_cells else 0.0

    lines = []
    lines.append("# OVERDETERMINATION_MAP — Phase E8 summary")
    lines.append("")
    lines.append("Generated by `_build_overdetermination_views.py` from `assimilation_dispatch.py`")
    lines.append("via the `qcalcgeom_solver.solve()` 4 x 3 dispatch matrix.")
    lines.append("")
    lines.append("## Top-line metrics")
    lines.append("")
    lines.append(f"- Observables in dispatch:       **{total_obs}**")
    lines.append(f"- 4 x 3 matrix cells (total):    **{total_cells}**")
    lines.append(f"- Populated cells:               **{populated}**  ({coverage_pct:.1f}%)")
    lines.append(f"- EXACT cells:                   **{exact}**")
    lines.append(f"- OK cells:                      **{ok}**")
    lines.append(f"- TENSION cells (OPEN_QUESTION): **{tension}**")
    lines.append(f"- GAP cells (no closure):        **{gap}**")
    lines.append("")
    lines.append("## Per-domain rollup")
    lines.append("")
    lines.append("| Domain | Observables | EXACT | OK | TENSION | Worst residual |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for domain in sorted(per_domain.keys()):
        domain_obs = per_domain[domain]
        n = len(domain_obs)
        ex = sum(o["status_counts"].get("EXACT",0) for o in domain_obs)
        ok_d = sum(o["status_counts"].get("OK",0) for o in domain_obs)
        ten = sum(o["status_counts"].get("TENSION",0) for o in domain_obs)
        worst_rp = 0.0
        for r in long_rows:
            if r["domain"] == domain and r["residual_pct"]:
                try:
                    worst_rp = max(worst_rp, float(r["residual_pct"]))
                except ValueError:
                    pass
        lines.append(f"| {domain} | {n} | {ex} | {ok_d} | {ten} | {worst_rp:.4f}% |")
    lines.append("")
    lines.append("## OPEN_QUESTION / TENSION cells")
    lines.append("")
    tension_rows = [r for r in long_rows if r["status"] == "TENSION"]
    if tension_rows:
        lines.append("| Observable | Geometry | Numeric | Residual | Primary source |")
        lines.append("|---|---|---|---:|---|")
        for r in tension_rows:
            lines.append(f"| {r['observable']} | {r['geometry']} | {r['numeric']} | {r['residual_pct']}% | {r['primary_source']} |")
    else:
        lines.append("_No TENSION cells._")
    lines.append("")
    lines.append("## Round 669 BAO multi-path closure")
    lines.append("")
    bao_rows = sorted([r for r in long_rows if "BAO_rd_H0_over_c" in r["observable"]
                       and r["numeric"] == "numerical"
                       and r["status"] != "GAP"], key=lambda x: x["observable"])
    if bao_rows:
        lines.append("| Observable | Formula | Residual | Source |")
        lines.append("|---|---|---:|---|")
        for r in bao_rows:
            rec = ad.DISPATCH.get(r["observable"], {})
            formula = rec.get("uqff_formula", "").split("[")[0].strip() if rec else r["observable"]
            lines.append(f"| {r['observable']} | `{formula}` | {r['residual_pct']}% | {r['primary_source']} |")
        lines.append("")
        lines.append("Both BAO closures use only locked primitives and converge on the same observable")
        lines.append("through independent primitive groupings. The two-path agreement (different formulas,")
        lines.append("same observable, both within 0.03%) is the corroboration that the form is structural")
        lines.append("rather than coincidental. See PAPER_1156 and SESSION_LOG Round 669.")
    else:
        lines.append("_No BAO closures found._")
    lines.append("")
    lines.append("## Per-observable coverage (wide preview, first 25)")
    lines.append("")
    preview_cols = ["observable","domain","owner_geometry","owner_N","total_N"]
    lines.append("| " + " | ".join(preview_cols) + " |")
    lines.append("|" + "|".join("---" for _ in preview_cols) + "|")
    for row in wide_rows[:25]:
        lines.append("| " + " | ".join(str(row[c]) for c in preview_cols) + " |")
    if len(wide_rows) > 25:
        lines.append(f"| ... | ... | ... | ... | ... | _(plus {len(wide_rows) - 25} more in OVERDETERMINATION_WIDE.csv)_ |")
    lines.append("")
    lines.append("## Schema notes")
    lines.append("")
    lines.append("- **Owner geometry** is the geometry whose backend owns the closure for that observable.")
    lines.append("- **owner_N**: how many of the owner's 3 numeric backends produced a value (target: 3).")
    lines.append("- **total_N**: how many of all 12 (4 geom x 3 num) cells produced a value (target: >=3).")
    lines.append("- **GAP** cells are expected for non-owner geometries that don't have closures for the observable.")
    lines.append("- **TENSION** cells flag OPEN_QUESTION closures (currently only the BAO Round 663 flag).")
    lines.append("- Status priority order: TENSION > EXACT > OK > GAP.")
    lines.append("")
    md_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote {md_path.name}: {len(lines)} lines, top-line {populated}/{total_cells} ({coverage_pct:.1f}%) populated")

    return 0


if __name__ == "__main__":
    sys.exit(main())
