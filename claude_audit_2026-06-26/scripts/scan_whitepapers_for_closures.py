#!/usr/bin/env python3
"""scan_whitepapers_for_closures.py — extract closure formulas across all 1,877 whitepapers.

READ-ONLY scan of /sessions/vibrant-keen-bohr/mnt/Star-Magic/whitepapers/*.md
Writes a JSON catalog to MY sandbox.
"""
import re
import json
from pathlib import Path
from collections import defaultdict

WP_DIR = Path("/sessions/vibrant-keen-bohr/mnt/Star-Magic/whitepapers")
OUT    = Path("/sessions/vibrant-keen-bohr/mnt/outputs/uqff_recompute/whitepaper_closures_index.json")
SUMMARY= Path("/sessions/vibrant-keen-bohr/mnt/outputs/uqff_recompute/WHITEPAPER_INDEX_MY_NOTES.md")

md_files = sorted(WP_DIR.glob("PAPER_*.md"))
print(f"Found {len(md_files):,} PAPER_*.md files in {WP_DIR}")

# index: paper_id -> list of (line, formula_match, quantity_tag)
index = {}
quantity_tags = defaultdict(list)

# regex for closure lines (= ... EXACT, residual_pct, match, etc.)
EXACT_RE  = re.compile(r"^\s*([\w\d\\^\s\(\)\.\+\-\*\/]{3,80})\s*=\s*([\w\d\\^\s\(\)\.\+\-\*\/]{1,80})\s*(?:EXACT|exact|≈|=|→|—|--)", re.MULTILINE)
BOXED_RE  = re.compile(r"\\boxed\{[^}]+\}")
RESID_RE  = re.compile(r"(?:residual|Match|match|Error|diff)\s*[:=]?\s*([\d\.\-+e]+%)")
PRIMARY_RE= re.compile(r"primary_source.*?PAPER_(\d{3,4})")
TAG_RE    = re.compile(r"tags:\s*\[([^\]]+)\]")
TITLE_RE  = re.compile(r"^title:\s*[\"']?([^\"'\n]+)[\"']?", re.MULTILINE)
EQ_RE     = re.compile(r"\$\$([^\$]{1,400})\$\$", re.DOTALL)

# Quantity tag mappings — keywords -> tag
KEYWORDS = {
    "lambda_cc":      ["Lambda", "Λ", "cosmological constant", "dark energy"],
    "yang_mills":     ["Yang-Mills", "yang_mills", "glueball", "mass gap", "Λ_QCD", "Lambda_QCD"],
    "magic_numbers":  ["magic number", "Mayer-Jensen", "shell model", "shell-model"],
    "holmlid_lenr":   ["Holmlid", "630 eV", "KER", "ultra-dense"],
    "rossi_ecat":     ["Rossi", "E-Cat", "COP"],
    "ssq":            ["SSq", "[SSq]", "0.57"],
    "s_26":           ["S_26", "S26", "Ramanujan"],
    "phi_res":        ["Phi_res", "Φ_res", "5/6"],
    "f_trz":          ["F_TRZ", "TRZ", "1/10"],
    "k_mex":          ["K_MEX", "K_Mex", "Mexican-hat"],
    "beta_i":         ["beta_i", "β_i", "0.6029"],
    "proton_mass":    ["m_p", "proton mass", "nucleon mass"],
    "alpha_fs":       ["fine structure", "alpha_FS", "1/137"],
    "h_planck":       ["Planck constant", "h_Planck", "6.626"],
    "g_newton":       ["Newton G", "G_Newton", "6.6743"],
    "c_light":        ["c_light", "speed of light", "299792458"],
    "rho_scm":        ["rho_SCm", "ρ_SCm", "7.09e-37"],
    "hubble":         ["H_0", "Hubble", "tension"],
    "neutron":        ["neutron lifetime", "879", "bottle method"],
    "millennium":     ["Millennium", "Clay", "Riemann", "Poincaré", "Hodge", "BSD", "Navier", "P vs NP"],
    "gw_events":      ["GW170817", "GW190425", "GW150914", "LIGO", "NANOGrav"],
    "agn":            ["AGN", "M87", "Sgr A", "TON618", "M-sigma"],
    "higgs":          ["Higgs", "m_H", "125", "γγ"],
    "lagrangian":     ["Lagrangian", "L_EH", "L_YM", "L_SCm", "L_KK", "Mexican-hat"],
    "26d_projection": ["26D", "26-D", "bosonic", "string critical", "Polyakov"],
    "dpm":            ["DPM", "Di-Pseudo-Monopole", "vortex"],
    "buoyancy":       ["F_U_Bi", "FUBi", "buoyancy", "F_UBii"],
    "caduceus":       ["Caduceus", "26 pinch", "π decimal"],
}

def tag_text(text_lower):
    tags = []
    for tag, words in KEYWORDS.items():
        for w in words:
            if w.lower() in text_lower:
                tags.append(tag)
                break
    return tags

print("\nScanning whitepapers for closure formulas, equations, and quantity tags...")
print(f"{'paper':<45} {'eqs':>4} {'tags'}")

total_eqs = 0
total_resids = 0
total_boxed = 0
errors = 0

for i, fp in enumerate(md_files):
    try:
        txt = fp.read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        errors += 1
        continue
    paper_id_m = re.search(r"PAPER_(\d{3,4})", fp.stem)
    paper_id = paper_id_m.group(1) if paper_id_m else fp.stem

    eqs = EQ_RE.findall(txt)
    boxed = BOXED_RE.findall(txt)
    resids = RESID_RE.findall(txt)
    title_m = TITLE_RE.search(txt)
    title = title_m.group(1).strip() if title_m else fp.stem
    tags = tag_text(txt.lower()[:50000])  # scan first 50k chars

    total_eqs += len(eqs)
    total_resids += len(resids)
    total_boxed += len(boxed)

    index[paper_id] = {
        "title": title,
        "filename": fp.name,
        "size_bytes": len(txt),
        "n_equations": len(eqs),
        "n_boxed": len(boxed),
        "n_residuals_reported": len(resids),
        "tags": tags,
        "sample_residuals": resids[:5],
    }
    for tag in tags:
        quantity_tags[tag].append(paper_id)

# Summary stats
print(f"\nTotal whitepapers scanned: {len(md_files):,}")
print(f"Total $$equations$$:        {total_eqs:,}")
print(f"Total \\boxed expressions:   {total_boxed:,}")
print(f"Total residuals reported:   {total_resids:,}")
print(f"Read errors:                 {errors}")
print()
print(f"{'tag':<22} {'papers'}")
print("-" * 60)
for tag in sorted(quantity_tags, key=lambda k: -len(quantity_tags[k])):
    print(f"{tag:<22} {len(quantity_tags[tag]):>5}  examples: {', '.join(quantity_tags[tag][:3])}")

# write JSON
out_data = {
    "summary": {
        "papers_scanned": len(md_files),
        "total_equations": total_eqs,
        "total_boxed": total_boxed,
        "total_residuals": total_resids,
        "errors": errors,
    },
    "quantity_tag_counts": {k: len(v) for k, v in quantity_tags.items()},
    "quantity_tag_papers": {k: v for k, v in quantity_tags.items()},
    "papers": index,
}
OUT.write_text(json.dumps(out_data, indent=2))
print(f"\nWrote JSON index: {OUT}")

# Markdown summary
lines = [f"# Whitepaper closures index (READ-ONLY scan)\n"]
lines.append(f"**Generated by Claude in MY sandbox; Daniel's repo untouched.**\n")
lines.append(f"\nTotal whitepapers: **{len(md_files):,}**\n")
lines.append(f"\n## By quantity tag\n")
lines.append("| Tag | Papers | Sample |")
lines.append("|---|---:|---|")
for tag in sorted(quantity_tags, key=lambda k: -len(quantity_tags[k])):
    sample = ", ".join(f"PAPER_{p}" for p in quantity_tags[tag][:5])
    lines.append(f"| `{tag}` | {len(quantity_tags[tag])} | {sample} |")
lines.append("\n## Top 30 papers by equation count\n")
sorted_by_eq = sorted(index.items(), key=lambda x: -x[1]["n_equations"])[:30]
lines.append("| Paper | Equations | Title |")
lines.append("|---|---:|---|")
for pid, d in sorted_by_eq:
    lines.append(f"| PAPER_{pid} | {d['n_equations']} | {d['title'][:80]} |")
lines.append("\n## All scanned papers in order\n")
lines.append("Total scanned (covering computational closures): see JSON for full table.")

SUMMARY.write_text("\n".join(lines))
print(f"Wrote summary: {SUMMARY}")
