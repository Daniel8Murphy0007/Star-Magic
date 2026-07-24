import csv
import hashlib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

REGISTRY_CSV = "UNIFIED_REGISTRY.csv"
GAPS_CSV = "UNIFIED_REGISTRY_GAPS.csv"
DUPES_CSV = "UNIFIED_REGISTRY_DUPLICATES.csv"

BASELINES = [
    "OVERDETERMINATION_MAP.csv",
    "OVERDETERMINATION_WIDE.csv",
    "MASTER_LEDGER_BY_CATEGORY.csv",
    "MASTER_LEDGER_BY_STATUS.csv",
    "MASTER_LEDGER_BY_SCRIPT.csv",
    "LEDGER_VS_PRIMITIVES_XREF.csv",
    "PRIMITIVES_RECONCILIATION.csv",
]

COLUMNS = [
    "quantity", "origin", "kind", "canonical_route", "formula", "value",
    "reference", "residual_pct", "sector", "phi_variant", "paper_source",
    "confirmations", "py_sites", "cpp_sites", "lean_sites",
    "corpus_citations", "status",
]


def norm(name):
    return re.sub(r"[^a-z0-9]+", "_", str(name).lower()).strip("_")


def paper_from_key(key):
    m = re.search(r"paper_?(\d{3,4}[a-z]{0,2})", key, re.I)
    return f"PAPER_{m.group(1).upper()}" if m else ""


def rows_from_dispatch(u):
    d = u.PARADOX_TO_CLOSURE
    fn_to_keys = {}
    for k, f in d.items():
        fn_to_keys.setdefault(f.__name__, []).append(k)
    rows, dupes = [], []
    for fname in sorted(fn_to_keys):
        keys = sorted(fn_to_keys[fname])
        primary = keys[0]
        try:
            out = d[primary]()
            val = out.get("primary_result", "") if isinstance(out, dict) else ""
            status = "derived"
        except Exception as e:
            val, status = "", f"ERROR:{type(e).__name__}"
        rows.append({
            "quantity": primary, "origin": "dispatch", "kind": "closure",
            "canonical_route": "", "formula": "", "value": repr(val),
            "reference": "", "residual_pct": "", "sector": "",
            "phi_variant": "", "paper_source": paper_from_key(primary + fname),
            "confirmations": len(keys) - 1, "py_sites": fname,
            "cpp_sites": "", "lean_sites": "", "corpus_citations": "",
            "status": status,
        })
        if len(keys) > 1:
            dupes.append({"quantity": primary, "kind": "dispatch_aliases",
                          "detail": ";".join(keys[1:])})
    return rows, dupes


def rows_from_odmap():
    rows, gaps = [], []
    with open("OVERDETERMINATION_MAP.csv", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            row = {
                "quantity": r.get("observable", ""), "origin": "odmap",
                "kind": "observable", "canonical_route": "",
                "formula": r.get("geometry", ""), "value": r.get("value", ""),
                "reference": r.get("target", ""),
                "residual_pct": r.get("residual_pct", ""),
                "sector": r.get("domain", ""), "phi_variant": "",
                "paper_source": r.get("primary_source", ""),
                "confirmations": "", "py_sites": r.get("owner_geometry", ""),
                "cpp_sites": "", "lean_sites": "", "corpus_citations": "",
                "status": r.get("status", ""),
            }
            rows.append(row)
            if r.get("status", "").upper() == "GAP":
                gaps.append({"quantity": r.get("observable", ""),
                             "kind": "odmap_gap",
                             "detail": r.get("primary_source", "")})
    return rows, gaps


def rows_from_primitives():
    rows = []
    with open("PRIMITIVES_RECONCILIATION.csv", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            rows.append({
                "quantity": r.get("primitive", ""), "origin": "primitives",
                "kind": "primitive", "canonical_route": "",
                "formula": "", "value": r.get("canonical_value", ""),
                "reference": "", "residual_pct": "",
                "sector": "", "phi_variant": "",
                "paper_source": r.get("canonical_source", ""),
                "confirmations": "", "py_sites": r.get("canonical_source", ""),
                "cpp_sites": "", "lean_sites": "", "corpus_citations": "",
                "status": "locked",
            })
    return rows


def rows_from_inventory():
    rows = []
    if not os.path.exists("CLOSED_CONSTANTS_INVENTORY.md"):
        return rows
    text = open("CLOSED_CONSTANTS_INVENTORY.md", encoding="utf-8").read()
    for m in re.finditer(r"^\*\*(\d+)\.\s+(.+?)\*\*", text, re.M):
        rows.append({
            "quantity": m.group(2).strip(), "origin": "inventory", "kind": "kernel_constant",
            "canonical_route": "", "formula": "", "value": "",
            "reference": "", "residual_pct": "", "sector": "",
            "phi_variant": "", "paper_source": "CLOSED_CONSTANTS_INVENTORY",
            "confirmations": "", "py_sites": "UQFFSimultaneousProofEngine",
            "cpp_sites": "", "lean_sites": "", "corpus_citations": "",
            "status": "CLOSED_flagged_see_2026_05_28_verification",
        })
    seen, out = set(), []
    for r in rows:
        k = norm(r["quantity"])
        if k not in seen:
            seen.add(k)
            out.append(r)
    return out


def load_citations():
    cites = {}
    if not os.path.exists("UNIFIED_REGISTRY_CORPUS_CITATIONS.csv"):
        return cites
    with open("UNIFIED_REGISTRY_CORPUS_CITATIONS.csv", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f):
            cites[r["paper"]] = r["cited_by_files"]
    return cites


def apply_citations(rows, cites):
    hit = 0
    for r in rows:
        m = re.search(r"PAPER[_ ]?(\d{3,4}[A-Z]{0,2})", str(r["paper_source"]), re.I)
        if m:
            key = f"PAPER_{m.group(1).upper()}"
            if key in cites:
                r["corpus_citations"] = cites[key]
                hit += 1
    return hit


def cross_origin_dupes(rows):
    by_norm = {}
    for r in rows:
        by_norm.setdefault(norm(r["quantity"]), set()).add(r["origin"])
    return [{"quantity": k, "kind": "cross_origin_match",
             "detail": ";".join(sorted(v))}
            for k, v in sorted(by_norm.items()) if len(v) > 1]


def load_cross_language_tokens():
    cpp_tokens, lean_tokens = set(), set()
    for cpp in ("uqff_exact_closures.cpp", "MAIN_1_CoAnQi.cpp"):
        if os.path.exists(cpp):
            src = open(cpp, encoding="utf-8", errors="ignore").read()
            for m in re.finditer(r"(?:^|\s)(?:static\s+)?(?:double|float|long double)\s+([A-Za-z_][A-Za-z0-9_]+)\s*\(", src):
                cpp_tokens.add(norm(m.group(1)))
    import glob as _g
    for lf in _g.glob("formal/**/*.lean", recursive=True) + _g.glob("formal/*.lean"):
        src = open(lf, encoding="utf-8", errors="ignore").read()
        for m in re.finditer(r"^\s*(?:theorem|def|lemma)\s+([A-Za-z_][A-Za-z0-9_'.]+)", src, re.M):
            lean_tokens.add(norm(m.group(1)))
    return cpp_tokens, lean_tokens


def apply_cross_language(rows, cpp_tokens, lean_tokens):
    cpp_hits = lean_hits = 0
    for r in rows:
        qn = norm(r["quantity"])
        qcore = qn.replace("_l96_uqff_", "").replace("_closure", "")
        for tok in sorted(cpp_tokens):
            if len(tok) > 8 and (tok in qcore or qcore in tok):
                r["cpp_sites"] = tok
                cpp_hits += 1
                break
        for tok in sorted(lean_tokens):
            if len(tok) > 6 and (tok in qcore or qcore in tok):
                r["lean_sites"] = tok
                lean_hits += 1
                break
    return cpp_hits, lean_hits


def compute_residuals(rows):
    hits = 0
    for r in rows:
        if r["residual_pct"]:
            continue
        try:
            v = float(str(r["value"]).strip("'\""))
            t = float(str(r["reference"]).strip("'\""))
            if t != 0.0:
                r["residual_pct"] = f"{abs(v - t) / abs(t) * 100.0:.6g}"
                hits += 1
        except (ValueError, TypeError):
            continue
    return hits


R1_CANONICAL_ROUTES = {
    "hbar": "PAPER_590 (R1 2026-07-22: physical-route standing rule; 1209EE S629 = confirmation)",
    "lambda": "PAPER_1156 (R228 precedent; PAPER_2094 = companion form)",
    "c_light": "PAPER_592 (1209EE S630 = confirmation)",
    "g_newton": "PAPER_593 (derive_G_newton crossing-form = confirmation)",
    "h0": "PAPER_2093 grid (2.27e-18 = observed local anchor; residual = Hubble tension per PAPER_2125)",
    "k_b": "1209EE S628 with Phi_5/6 (PAPER_2129)",
    "phi_res": "SECTOR RULE per PAPER_2129 (R1 2026-07-22: counting sectors 5/6, resonance sectors 0.84)",
}


PHI_COUNTING_TOKENS = ("magic", "shell", "nuclear", "binding", "thermo", "entropy",
                       "boltzmann", "k_b", "avogadro", "faraday", "occupancy")
PHI_RESONANCE_TOKENS = ("lenr", "holmlid", "630", "phonon", "k_spring", "quantum_chain",
                        "ker", "rossi", "parkhomov", "mizuno", "resonan", "s26", "s_26")


def apply_r1_verdicts(rows):
    route_hits = gap_marks = sole_marks = phi_marks = 0
    for r in rows:
        qn = norm(r["quantity"])
        matched = False
        for key, verdict in R1_CANONICAL_ROUTES.items():
            if key == qn or key in qn:
                r["canonical_route"] = verdict
                route_hits += 1
                matched = True
                break
        if not matched and not r["canonical_route"]:
            r["canonical_route"] = "SOLE_ROUTE (auto R1-completion 2026-07-22)"
            sole_marks += 1
        if r["origin"] == "odmap" and str(r["status"]).upper() == "GAP":
            r["status"] = "SYMBOLIC_PENDING_R1"
            gap_marks += 1
        if not r["phi_variant"]:
            hay = qn + "_" + norm(r.get("sector", ""))
            if any(t in hay for t in PHI_COUNTING_TOKENS):
                r["phi_variant"] = "5/6 (counting sector, PAPER_2129 rule)"
                phi_marks += 1
            elif any(t in hay for t in PHI_RESONANCE_TOKENS):
                r["phi_variant"] = "0.84 (resonance sector, PAPER_2129 rule)"
                phi_marks += 1
    return route_hits, gap_marks, sole_marks, phi_marks


MERGE_KEYS = {
    "G_newton": ("g_newton", "gravitational_constant", "newton_g"),
    "c_light": ("c_light", "speed_of_light", "lightspeed"),
    "hbar": ("hbar", "planck_reduced", "h_bar"),
    "planck_h": ("planck_h", "planck_constant"),
    "k_B": ("k_b", "boltzmann"),
    "Lambda": ("lambda", "cosmological_constant"),
    "H0": ("h0", "hubble"),
    "mu_0": ("mu_0", "permeability"),
    "beta_i": ("beta_i",),
    "rho_SCm": ("rho_scm", "rho_vac_scm"),
    "rho_UA": ("rho_ua", "rho_vac_ua", "rho_vac"),
    "kappa": ("kappa",),
    "Phi_res": ("phi_res", "phi_resonance"),
    "F_TRZ": ("f_trz", "trz"),
    "SSq": ("ssq",),
    "omega_SCm": ("omega_scm", "1p25_thz", "1_25_thz"),
}


def write_merged_view(rows):
    groups = {}
    for r in rows:
        qn = norm(r["quantity"])
        for canon, tokens in MERGE_KEYS.items():
            if any(t in qn for t in tokens):
                groups.setdefault(canon, []).append(r)
                break
    out = []
    for canon in sorted(groups):
        g = groups[canon]
        origins = sorted(set(r["origin"] for r in g))
        routes = sorted(set(r["canonical_route"] for r in g if r["canonical_route"] and "SOLE_ROUTE" not in r["canonical_route"]))
        out.append({
            "quantity_key": canon,
            "n_rows": len(g),
            "origins": ";".join(origins),
            "canonical_route": routes[0] if routes else "SOLE_ROUTE_family",
            "row_names_sample": ";".join(sorted(r["quantity"] for r in g)[:5]),
        })
    with open("UNIFIED_REGISTRY_MERGED.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["quantity_key", "n_rows", "origins", "canonical_route", "row_names_sample"])
        w.writeheader()
        w.writerows(out)
    return len(out), sum(g["n_rows"] for g in out)


def main():
    import uqff_pure_calculator as u
    d_rows, dupes = rows_from_dispatch(u)
    o_rows, gaps = rows_from_odmap()
    p_rows = rows_from_primitives()
    i_rows = rows_from_inventory()
    rows = p_rows + d_rows + o_rows + i_rows
    dupes += cross_origin_dupes(rows)
    cite_hits = apply_citations(rows, load_citations())
    cpp_hits, lean_hits = apply_cross_language(rows, *load_cross_language_tokens())
    resid_hits = compute_residuals(rows)
    route_hits, gap_marks, sole_marks, phi_marks = apply_r1_verdicts(rows)
    merged_groups, merged_rows = write_merged_view(rows)
    print(f"cross_language: cpp_hits={cpp_hits} lean_hits={lean_hits} | residuals_computed={resid_hits}")
    print(f"R1 verdicts applied: canonical_routes={route_hits} gap_rows_marked_SYMBOLIC_PENDING_R1={gap_marks}")
    print(f"R1-completion: sole_routes={sole_marks} phi_variants_assigned={phi_marks} | merged_view: {merged_groups} groups covering {merged_rows} rows")
    rows.sort(key=lambda r: (r["origin"], norm(r["quantity"])))

    with open(REGISTRY_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    with open(GAPS_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["quantity", "kind", "detail"])
        w.writeheader()
        w.writerows(gaps)
    with open(DUPES_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["quantity", "kind", "detail"])
        w.writeheader()
        w.writerows(dupes)

    h = hashlib.sha256(open(REGISTRY_CSV, "rb").read()).hexdigest()
    print(f"rows={len(rows)} (primitives={len(p_rows)} dispatch={len(d_rows)} "
          f"odmap={len(o_rows)} inventory={len(i_rows)}) citation_hits={cite_hits}")
    print(f"gaps={len(gaps)} duplicates={len(dupes)}")
    print(f"registry_sha256={h}")
    errors = [r for r in rows if str(r["status"]).startswith("ERROR")]
    print(f"closure_call_errors={len(errors)}")
    for r in errors[:5]:
        print("  ERR:", r["quantity"], r["status"])
    print("BASELINE_HASHES:")
    for b in BASELINES:
        bh = hashlib.sha256(open(b, "rb").read()).hexdigest()[:16]
        print(f"  {b}: {bh}")


if __name__ == "__main__":
    main()
