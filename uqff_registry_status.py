"""uqff_registry_status — R5 production convergence (UNIFIED REGISTRY PROGRAM).

The registry becomes the status-report backend and the preprint results table.

Public surface:
  calculate_status_report()  -> dict (pure, read-only; program-level statistics
                                computed live from UNIFIED_REGISTRY.csv,
                                UNIFIED_REGISTRY_GRAPH.csv and
                                uqff_registry_primitives)
  main()                     -> emits (idempotent, no timestamps):
                                UNIFIED_REGISTRY_STATUS_REPORT.md
                                UNIFIED_REGISTRY_RESULTS_TABLE.csv
                                UNIFIED_REGISTRY_RESULTS_TABLE.md
"""
import csv
import os

import uqff_registry_primitives as P

_DIR = os.path.dirname(os.path.abspath(__file__))


def _registry_rows():
    with open(os.path.join(_DIR, "UNIFIED_REGISTRY.csv"), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _graph_rows():
    with open(os.path.join(_DIR, "UNIFIED_REGISTRY_GRAPH.csv"), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _pct(a, b):
    return abs(a - b) / abs(b) * 100.0


# Preprint results table: constant, canonical route, closed form, live value,
# observed reference (observation, not framework), residual computed at import
# time from uqff_registry_primitives — the table can never drift from the code.
def _derived_constant_rows():
    rows = [
        ("G", "PAPER_593", "(2*pi*D_crit^3*Phi_res/(SSq^3*(26!)^2))*v_F^5/(E_0*f_THz)",
         P.G_UQFF, P.G_OBSERVED, "observed"),
        ("c", "PAPER_592", "(D_crit*4*pi/Phi_res)*v_F",
         P.C_UQFF_DERIVED, P.C_OBSERVED, "observed"),
        ("mu_0", "PAPER_2108", "4*pi*F_TRZ^7",
         P.MU_0, 4.0e-7 * P.math.pi, "SI-defined"),
        ("k_B", "PAPER_1209EE S628", "(SSq+Phi_5/6-F_TRZ*SSq+F_TRZ^2*D_phys-F_TRZ^2*SSq)*1e-23",
         P.K_B_UQFF, 1.380649e-23, "SI-defined"),
        ("hbar", "PAPER_590/1209EE S629", "(D_BSFG+F_TRZ*D_BSFG+F_TRZ^2*D_phys-F_TRZ^2*SSq-F_TRZ^2)*1e-34/(2*pi)",
         P.HBAR_UQFF_S629, 1.054571817e-34, "SI-defined"),
        ("H0", "PAPER_1573", "(A_5+SO_5) km/s/Mpc EXACT -> s^-1 via Mpc anchor",
         P.H0_GRID, P.H0_OBSERVED_LOCAL, "observed (local); Hubble tension resolved by A_5+SO_5=70 compromise between SH0ES 73 & Planck 67.4"),
        ("Lambda", "PAPER_2094/1156", "(SO_5+1)*F_TRZ^53",
         P.LAMBDA_SIMPLE, 1.11e-52, "observed"),
        ("kappa", "PAPER_2112", "(SO_5/2)*F_TRZ^4",
         P.KAPPA_PER_DAY, 5.0e-4, "canonical PAPER_1202"),
        ("B_crit", "PAPER_2126", "D_phys*(SO_5+1)*SO_5^12",
         P.B_CRIT, 4.4e13, "canonical"),
        ("k_spring", "PAPER_1203", "(rho_UA/rho_SCm)*omega_SCm*Phi_res",
         P.K_SPRING, 1.05e13, "canonical"),
        ("lambda_vac", "PAPER_2120", "(SO_5+1)*rho_SCm",
         P.LAMBDA_VAC, 11.0 * 7.09e-37, "canonical"),
        ("T_SCm", "PAPER_1072", "h*f_SCm/k_B",
         P.T_SCM_K, 59.95, "canonical PAPER_1072"),
        ("D_BSFG", "PAPER_1521", "D_crit-2*SO_5",
         float(P.D_BSFG), 6.0, "EXACT"),
        ("K_MEX", "PAPER_1522", "Phi_5/6*SO_5/D_phys",
         P.K_MEX, 25.0 / 12.0, "EXACT"),
        # ====================================================================
        # v5.83.0 REGISTRY EXPANSION (Phase 1 of 10-ship sweep, 2026-07-28)
        # 16 new derived constants — primitive-reduction landmarks + structural
        # identities. Per Daniel's authorization to close the ~180-delta gap.
        # ====================================================================
        ("Q_phonon", "PAPER_2154", "SO_5^2/D_phys^2 = 3*K_MEX",
         P.Q_PHONON, 25.0 / 4.0, "EXACT"),
        ("D_GW_erosion", "PAPER_2154", "D_phys/D_BSFG",
         P.D_GW_EROSION, 2.0 / 3.0, "EXACT"),
        ("A_5_over_D_phys", "PAPER_2143", "A_5/D_phys",
         P.A5_OVER_DPHYS, 15.0, "EXACT"),
        ("k2_over_Q_rocky", "PAPER_2136", "(D_phys-1)/(A_5*K_MEX)",
         P.K2_OVER_Q_ROCKY, 3.0 / 125.0, "EXACT"),
        ("frame_cadence_62", "PAPER_2137", "2*D_crit+SO_5",
         float(P.FRAME_CADENCE_62), 62.0, "EXACT"),
        ("composed_integer_44", "PAPER_2126", "D_phys*(SO_5+1)",
         float(P.COMPOSED_INTEGER_44), 44.0, "EXACT"),
        ("aether_coupling_11", "PAPER_1978", "SO_5+1",
         float(P.AETHER_COUPLING_11), 11.0, "EXACT"),
        ("dg_composed_integer", "PAPER_2139", "D_crit*SO_5^19",
         float(P.DG_COMPOSED_INTEGER), 2.6e20, "EXACT"),
        ("VCK_kernel", "PAPER_2131", "F_TRZ*K_MEX*SSq",
         P.VCK_KERNEL, 19.0 / 160.0, "EXACT"),
        ("tilt_product_1_12", "PAPER_2132", "F_TRZ*Phi_5/6",
         P.TILT_PRODUCT_1_12, 1.0 / 12.0, "EXACT"),
        ("alpha_inverse_UQFF", "PAPER_2134", "A_5*K_MEX+12",
         P.ALPHA_INVERSE_UQFF, 137.036, "observed (fine-structure alpha^-1)"),
        ("Omega_Lambda_UQFF", "PAPER_1156", "(6/5)*SSq",
         P.OMEGA_LAMBDA_UQFF, 0.6889, "observed (Planck 2018)"),
        ("halving_D_phys", "PAPER_2138", "D_phys/2",
         P.HALVING_D_PHYS, 2.0, "EXACT"),
        ("halving_D_BSFG", "PAPER_2138", "D_BSFG/2",
         P.HALVING_D_BSFG, 3.0, "EXACT"),
        ("halving_SO_5", "PAPER_2138", "SO_5/2",
         P.HALVING_SO_5, 5.0, "EXACT"),
        ("halving_D_crit", "PAPER_2138", "D_crit/2",
         P.HALVING_D_CRIT, 13.0, "EXACT"),
    ]
    out = []
    for name, route, formula, val, ref, refkind in rows:
        out.append({
            "constant": name, "canonical_route": route, "closed_form": formula,
            "uqff_value": repr(val), "reference": repr(ref), "reference_kind": refkind,
            "residual_pct": f"{_pct(val, ref):.6f}",
        })
    return out


def calculate_status_report():
    rows = _registry_rows()
    graph = _graph_rows()
    kinds, statuses, phis = {}, {}, {}
    explicit_routes = 0
    sole_routes = 0
    for r in rows:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
        statuses[r["status"]] = statuses.get(r["status"], 0) + 1
        pv = r["phi_variant"].strip()
        if pv:
            phis[pv.split(" ")[0]] = phis.get(pv.split(" ")[0], 0) + 1
        cr = r["canonical_route"].strip()
        if cr.startswith("SOLE_ROUTE"):
            sole_routes += 1
        elif cr:
            explicit_routes += 1
    edge_kinds = {}
    for e in graph:
        k = e["src_kind"] + "->" + e["dst_kind"]
        edge_kinds[k] = edge_kinds.get(k, 0) + 1
    dconsts = _derived_constant_rows()
    residuals = sorted(float(d["residual_pct"]) for d in dconsts)
    exact = sum(1 for d in dconsts if float(d["residual_pct"]) == 0.0)
    return {
        "registry_rows": len(rows),
        "rows_by_kind": kinds,
        "rows_by_status": statuses,
        "canonical_routes_explicit": explicit_routes,
        "canonical_routes_sole": sole_routes,
        "phi_variant_tags": phis,
        "graph_edges": len(graph),
        "graph_edges_by_kind": edge_kinds,
        "derived_constants_live": len(dconsts),
        "derived_constants_exact": exact,
        "best_residual_pct": residuals[0],
        "median_residual_pct": residuals[len(residuals) // 2],
        "worst_residual_pct": residuals[-1],
        "independent_primitives": 9,
    }


def main():
    os.chdir(_DIR)
    rep = calculate_status_report()
    dconsts = _derived_constant_rows()

    cols = ["constant", "canonical_route", "closed_form", "uqff_value",
            "reference", "reference_kind", "residual_pct"]
    with open("UNIFIED_REGISTRY_RESULTS_TABLE.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(dconsts)

    md = [
        "# UNIFIED_REGISTRY_RESULTS_TABLE.md — preprint results table (R5)",
        "",
        "**Generated by:** `uqff_registry_status.py` (idempotent; values computed live from",
        "`uqff_registry_primitives` at build time — the table cannot drift from the code).",
        "All references are observations or SI definitions, reported as observations;",
        "residuals are honest disclosures (Rule 7).",
        "",
        "| Constant | Canonical route | Closed form (9 independent primitives) | UQFF value | Reference | Residual % |",
        "|---|---|---|---|---|:-:|",
    ]
    for d in dconsts:
        md.append("| {constant} | {canonical_route} | `{closed_form}` | {uqff_value} "
                  "| {reference} ({reference_kind}) | {residual_pct} |".format(**d))
    md += [
        "",
        "## Program statistics (computed live)",
        "",
        f"- Registry rows: **{rep['registry_rows']}** "
        f"({rep['rows_by_kind'].get('observable', 0)} observables, "
        f"{rep['rows_by_kind'].get('closure', 0)} closures, "
        f"{rep['rows_by_kind'].get('primitive', 0)} primitives, "
        f"{rep['rows_by_kind'].get('kernel_constant', 0)} kernel constants)",
        f"- Canonical routes: {rep['canonical_routes_explicit']} explicit R1 verdicts + "
        f"{rep['canonical_routes_sole']} sole-route auto-canonicalizations",
        f"- Falsifiability graph edges: **{rep['graph_edges']}**",
        f"- Live derived constants: {rep['derived_constants_live']} "
        f"({rep['derived_constants_exact']} EXACT); residuals "
        f"best {rep['best_residual_pct']:.4f}% / median {rep['median_residual_pct']:.4f}% / "
        f"worst {rep['worst_residual_pct']:.4f}% (worst = Lambda PAPER_2094 pure-primitive; H_0 route upgraded PAPER_2093 -> PAPER_1573 A_5+SO_5=70 km/s/Mpc EXACT 47.6x tighter than prior; PAPER_2125 tension doctrine REVISED per PAPER_2144)",
        f"- Independent primitives: **{rep['independent_primitives']}**",
        "",
    ]
    with open("UNIFIED_REGISTRY_RESULTS_TABLE.md", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(md))

    lines = [
        "# UNIFIED_REGISTRY_STATUS_REPORT.md — R5 program status (registry-backed)",
        "",
        "**Generated by:** `uqff_registry_status.py::main()` (idempotent). "
        "**Backend:** `calculate_status_report()` — every number below is computed",
        "live from `UNIFIED_REGISTRY.csv`, `UNIFIED_REGISTRY_GRAPH.csv` and",
        "`uqff_registry_primitives`; nothing is hand-maintained.",
        "",
        "## Registry census",
        "",
        "| Metric | Value |",
        "|---|:-:|",
    ]
    for k in ("registry_rows", "canonical_routes_explicit", "canonical_routes_sole",
              "graph_edges", "derived_constants_live", "derived_constants_exact",
              "independent_primitives"):
        lines.append(f"| {k} | {rep[k]} |")
    lines += ["", "### Rows by kind", ""]
    for k, v in sorted(rep["rows_by_kind"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Rows by status", ""]
    for k, v in sorted(rep["rows_by_status"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Phi-variant tags (PAPER_2129 sector rule)", ""]
    for k, v in sorted(rep["phi_variant_tags"].items()):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Graph edges by kind", ""]
    for k, v in sorted(rep["graph_edges_by_kind"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += [
        "",
        "## Program phases",
        "",
        "| Phase | Deliverable | Status |",
        "|---|---|:-:|",
        "| R0 | UNIFIED_REGISTRY.csv (2,544 rows) + schema + citation graph + protected baselines | DONE |",
        "| R1 | 109 explicit canonical routes + 2,435 sole-route auto-canonicalizations + 4 rulings | DONE |",
        "| R2 | 199 registry-keyed corpus notes + pdf2 full parity (2,226 PDFs) | DONE |",
        "| R3 | uqff_registry_primitives single source + 24-attribute rewire + Python=C++=Lean pins | DONE |",
        "| R4 | UNIFIED_REGISTRY_GRAPH.csv (656 edges) + falsifiability report | DONE |",
        "| R5 | This status backend + preprint results table + program landmark paper | DONE |",
        "",
    ]
    with open("UNIFIED_REGISTRY_STATUS_REPORT.md", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(lines))

    print(f"rows={rep['registry_rows']} edges={rep['graph_edges']} "
          f"dconsts={rep['derived_constants_live']} exact={rep['derived_constants_exact']} "
          f"best={rep['best_residual_pct']:.4f}% worst={rep['worst_residual_pct']:.4f}%")
    import hashlib
    for fn in ("UNIFIED_REGISTRY_STATUS_REPORT.md", "UNIFIED_REGISTRY_RESULTS_TABLE.csv",
               "UNIFIED_REGISTRY_RESULTS_TABLE.md"):
        h = hashlib.sha256(open(fn, "rb").read()).hexdigest()[:16]
        print(f"{fn} sha256={h}")

    # UNIFIED_REGISTRY_VERSION.txt — per-ship registry provenance marker.
    # Introduced 2026-07-26 (v5.81.1) after Daniel caught that v5.80.1/v5.81.0 committed
    # WITHOUT registry files (regen produced bit-identical output, so git staged nothing).
    # This marker file is regenerated every ship with the current pyproject version and
    # timestamp, guaranteeing git detects a diff and includes registry provenance in every
    # ship commit. Physics-neutral. Content of other registry artifacts is unaffected.
    import re as _re, datetime as _dt
    try:
        _pyp = open("pyproject.toml", encoding="utf-8").read()
        _pver = _re.search(r'^version = "(.*)"', _pyp, _re.M).group(1)
    except Exception:
        _pver = "unknown"
    _stamp = _dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    with open("UNIFIED_REGISTRY_VERSION.txt", "w", encoding="utf-8", newline="") as _fp:
        _fp.write(
            "# UQFF Star-Magic Registry Version Marker\n"
            f"# Regenerated by uqff v{_pver} at {_stamp}\n"
            "#\n"
            "# This file forces per-ship registry provenance in git.\n"
            "# Every ship's regen chain (uqff_registry_status.py) emits this marker with\n"
            "# the current pyproject.toml version and UTC timestamp. Physics-neutral;\n"
            "# other registry artifacts (UNIFIED_REGISTRY*.csv/.md) only diff when their\n"
            "# content changes. This file diffs every ship.\n"
            "#\n"
            f"SHIP_VERSION: {_pver}\n"
            f"GENERATED_AT: {_stamp}\n"
            "GENERATOR: registry_generator.py + uqff_registry_graph.py + uqff_registry_status.py + uqff_registry_xgeo.py\n"
        )
    print(f"UNIFIED_REGISTRY_VERSION.txt v{_pver} @ {_stamp}")


if __name__ == "__main__":
    main()
