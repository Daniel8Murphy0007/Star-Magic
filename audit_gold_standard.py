#!/usr/bin/env python3
"""
Gold Standard Audit Scanner — PAPER_001 through PAPER_1012
Checks all 6 CVW v2.0.0 gates + §A Cosmogenesis + §B VDS/DVP/BSH appendices.
Outputs: audit_gold_standard_results.json with per-paper scores and deficiency lists.
"""

import os
import re
import json
import glob
from pathlib import Path

WHITEPAPER_DIR = os.path.join(os.path.dirname(__file__), "whitepapers")
RANGE_START = 1
RANGE_END = 1012

# --- Gate Check Functions ---

def check_g1_header(content: str) -> dict:
    """G1: Header Completeness — Author, Framework, Session, Date"""
    issues = []
    has_author = bool(re.search(r'\*\*Author\*\*.*Daniel|author:.*Daniel', content, re.I))
    has_framework = bool(re.search(r'Framework.*UQFF|framework.*UQFF|κ\s*=|kappa\s*=', content, re.I))
    has_session = bool(re.search(r'\*\*Session\*\*|session:', content, re.I))
    has_title = bool(re.search(r'^#\s+PAPER_\d+', content, re.M))
    if not has_author:
        issues.append("missing_author")
    if not has_framework:
        issues.append("missing_framework_constants")
    if not has_session:
        issues.append("missing_session")
    if not has_title:
        issues.append("missing_paper_title_header")
    score = sum([has_author, has_framework, has_session, has_title])
    return {"pass": score >= 3, "score": score, "max": 4, "issues": issues}


def check_g2_abstract(content: str) -> dict:
    """G2: Abstract — ≥3 sentences, 100+ words ideal"""
    issues = []
    abstract_match = re.search(
        r'(?:^##\s*Abstract|^Abstract)\s*\n(.*?)(?=\n##|\n---|\Z)',
        content, re.M | re.S | re.I
    )
    if not abstract_match:
        return {"pass": False, "score": 0, "max": 3, "issues": ["no_abstract_section"]}
    text = abstract_match.group(1).strip()
    sentences = [s.strip() for s in re.split(r'[.!?]+', text) if len(s.strip()) > 10]
    words = len(text.split())
    if len(sentences) < 1:
        issues.append("abstract_empty")
    elif len(sentences) < 3:
        issues.append(f"abstract_short_{len(sentences)}_sentences")
    if words < 50:
        issues.append(f"abstract_only_{words}_words")
    score = min(3, len(sentences))
    return {"pass": len(sentences) >= 2 and words >= 30, "score": score, "max": 3, "issues": issues}


def check_g3_equations(content: str) -> dict:
    """G3: Core Equations — ≥2 LaTeX $$ blocks"""
    issues = []
    # Match $$ ... $$ blocks
    eq_blocks = re.findall(r'\$\$.*?\$\$', content, re.S)
    # Also match ``` math blocks or indented equation patterns
    code_eqs = re.findall(r'```(?:math|latex)(.*?)```', content, re.S | re.I)
    total = len(eq_blocks) + len(code_eqs)
    if total < 1:
        issues.append("no_equations")
    elif total < 2:
        issues.append(f"only_{total}_equation_block")
    return {"pass": total >= 1, "score": min(3, total), "max": 3, "issues": issues}


def check_g4_numerical(content: str) -> dict:
    """G4: Numerical Results — e-notation values with units"""
    issues = []
    # Match patterns like 1.23e-22, 1.23×10⁻²², 1.23 × 10^{-22}, 10⁻⁴², etc.
    e_notation = re.findall(r'\d+\.?\d*[eE][+-]?\d+', content)
    sci_notation = re.findall(r'\d+\.?\d*\s*[×x]\s*10[\^⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹{}\-\d]+', content)
    power_notation = re.findall(r'10[\^{]\s*[-]?\d+', content)
    total = len(e_notation) + len(sci_notation) + len(power_notation)
    if total < 1:
        issues.append("no_numerical_results")
    elif total < 3:
        issues.append(f"only_{total}_numerical_values")
    # Check for unicode superscript issues (×10⁻²² should be e-notation)
    unicode_sci = re.findall(r'×10[⁻⁺]?[⁰¹²³⁴⁵⁶⁷⁸⁹]+', content)
    if unicode_sci:
        issues.append(f"unicode_superscripts_{len(unicode_sci)}")
    return {"pass": total >= 2, "score": min(3, total), "max": 3, "issues": issues}


def check_g5_crossrefs(content: str) -> dict:
    """G5: Anchor Cross-References — ≥2 PAPER_NNN citations"""
    issues = []
    refs = set(re.findall(r'PAPER_(\d+)', content))
    if len(refs) < 1:
        issues.append("no_cross_references")
    elif len(refs) < 2:
        issues.append(f"only_{len(refs)}_cross_reference")
    return {"pass": len(refs) >= 2, "score": min(3, len(refs)), "max": 3, "issues": issues}


def check_g6_sm_anchors(content: str, paper_num: int) -> dict:
    """G6: SM Anchor section — mandatory for PAPER_422+"""
    issues = []
    mandatory = paper_num >= 422
    has_sm = bool(re.search(r'§SM|SM Anchor|Standard Model Cross-Validation|sm_anchor', content, re.I))
    has_table = bool(re.search(r'\|\s*Observable\s*\||\|\s*UQFF\s*Prediction\s*\||\|\s*Alignment\s*\|', content, re.I))
    if not has_sm:
        if mandatory:
            issues.append("MISSING_SM_ANCHORS_MANDATORY")
        else:
            issues.append("no_sm_anchors_optional")
    elif not has_table:
        issues.append("sm_anchors_no_table")
    score = 2 if (has_sm and has_table) else (1 if has_sm else 0)
    passing = has_sm if mandatory else True
    return {"pass": passing, "score": score, "max": 2, "issues": issues, "mandatory": mandatory}


def check_cosmogenesis(content: str) -> dict:
    """§A Cosmogenesis appendix check"""
    issues = []
    has_cosmo = bool(re.search(r'§A|Cosmogenesis|Lagrangian.*Density|Euler-Lagrange|cosmogenesis', content, re.I))
    has_sectors = bool(re.search(r'Sector Classification|sector.*Lagrangian|L_\{|\\mathcal\{L\}', content, re.I))
    if not has_cosmo:
        issues.append("MISSING_COSMOGENESIS_APPENDIX")
    return {"present": has_cosmo, "has_sectors": has_sectors, "issues": issues}


def check_vds_dvp_bsh(content: str) -> dict:
    """§B VDS/DVP/BSH Deep Synthesis check"""
    issues = []
    has_vds = bool(re.search(r'VDS|Vacuum Density Series|VDS sub-ratio', content, re.I))
    has_dvp = bool(re.search(r'DVP|Dipole Vortex Prime', content, re.I))
    has_bsh = bool(re.search(r'BSH|Buoyancy.*Saturation.*Harmonic|BSH.*timescale', content, re.I))
    has_section = bool(re.search(r'§B|VDS/DVP/BSH', content))
    count = sum([has_vds, has_dvp, has_bsh])
    if count == 0:
        issues.append("MISSING_VDS_DVP_BSH")
    elif count < 3:
        missing = []
        if not has_vds: missing.append("VDS")
        if not has_dvp: missing.append("DVP")
        if not has_bsh: missing.append("BSH")
        issues.append(f"partial_VDS_DVP_BSH_missing_{'_'.join(missing)}")
    return {"present": has_section or count >= 2, "vds": has_vds, "dvp": has_dvp, "bsh": has_bsh, "issues": issues}


def check_content_depth(content: str) -> dict:
    """Check content depth — line count, word count, section count"""
    lines = content.count('\n') + 1
    words = len(content.split())
    sections = len(re.findall(r'^##\s', content, re.M))
    is_stub = lines < 100 or words < 500
    return {"lines": lines, "words": words, "sections": sections, "is_stub": is_stub}


def compute_score(g1, g2, g3, g4, g5, g6, cosmo, vds, depth) -> int:
    """Compute 0-100 quality score"""
    score = 0
    # G1 header: 15 pts
    score += int((g1["score"] / g1["max"]) * 15)
    # G2 abstract: 15 pts
    score += int((g2["score"] / g2["max"]) * 15)
    # G3 equations: 15 pts
    score += int((g3["score"] / g3["max"]) * 15)
    # G4 numerical: 10 pts
    score += int((g4["score"] / g4["max"]) * 10)
    # G5 cross-refs: 10 pts
    score += int((g5["score"] / g5["max"]) * 10)
    # G6 SM anchors: 10 pts
    score += int((g6["score"] / g6["max"]) * 10)
    # §A Cosmogenesis: 10 pts
    score += 10 if cosmo["present"] else 0
    # §B VDS/DVP/BSH: 10 pts
    score += 10 if vds["present"] else 0
    # Content depth: 5 pts
    score += 0 if depth["is_stub"] else 5
    return min(100, score)


def audit_paper(filepath: str, paper_num: int) -> dict:
    """Audit a single paper against gold standard"""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    g1 = check_g1_header(content)
    g2 = check_g2_abstract(content)
    g3 = check_g3_equations(content)
    g4 = check_g4_numerical(content)
    g5 = check_g5_crossrefs(content)
    g6 = check_g6_sm_anchors(content, paper_num)
    cosmo = check_cosmogenesis(content)
    vds = check_vds_dvp_bsh(content)
    depth = check_content_depth(content)
    score = compute_score(g1, g2, g3, g4, g5, g6, cosmo, vds, depth)

    # Collect all issues
    all_issues = []
    all_issues.extend(g1["issues"])
    all_issues.extend(g2["issues"])
    all_issues.extend(g3["issues"])
    all_issues.extend(g4["issues"])
    all_issues.extend(g5["issues"])
    all_issues.extend(g6["issues"])
    all_issues.extend(cosmo["issues"])
    all_issues.extend(vds["issues"])

    # Determine tier
    if score >= 90:
        tier = "GOLD"
    elif score >= 70:
        tier = "FULL"
    elif score >= 50:
        tier = "STANDARD"
    elif score >= 30:
        tier = "BRIEF"
    else:
        tier = "STUB"

    return {
        "paper": f"PAPER_{paper_num:03d}" if paper_num < 1000 else f"PAPER_{paper_num}",
        "score": score,
        "tier": tier,
        "lines": depth["lines"],
        "words": depth["words"],
        "sections": depth["sections"],
        "is_stub": depth["is_stub"],
        "gates": {
            "G1_header": g1["pass"],
            "G2_abstract": g2["pass"],
            "G3_equations": g3["pass"],
            "G4_numerical": g4["pass"],
            "G5_crossrefs": g5["pass"],
            "G6_sm_anchors": g6["pass"],
        },
        "appendices": {
            "cosmogenesis": cosmo["present"],
            "vds_dvp_bsh": vds["present"],
        },
        "issues": all_issues,
        "critical": [i for i in all_issues if i.startswith("MISSING_")],
    }


def main():
    results = []
    missing_files = []

    # Build a map of paper number -> filepath from actual directory listing
    paper_map = {}
    if os.path.isdir(WHITEPAPER_DIR):
        for fname in os.listdir(WHITEPAPER_DIR):
            if not fname.endswith('.md'):
                continue
            m = re.match(r'PAPER_(\d+)', fname)
            if m:
                num = int(m.group(1))
                paper_map[num] = os.path.join(WHITEPAPER_DIR, fname)

    for num in range(RANGE_START, RANGE_END + 1):
        filepath = paper_map.get(num)
        if not filepath:
            missing_files.append(num)
            continue

        result = audit_paper(filepath, num)
        results.append(result)

    # --- Summary Statistics ---
    total = len(results)
    tiers = {"GOLD": 0, "FULL": 0, "STANDARD": 0, "BRIEF": 0, "STUB": 0}
    gate_pass = {"G1": 0, "G2": 0, "G3": 0, "G4": 0, "G5": 0, "G6": 0}
    appendix_present = {"cosmogenesis": 0, "vds_dvp_bsh": 0}
    needs_upgrade = []

    for r in results:
        tiers[r["tier"]] += 1
        if r["gates"]["G1_header"]: gate_pass["G1"] += 1
        if r["gates"]["G2_abstract"]: gate_pass["G2"] += 1
        if r["gates"]["G3_equations"]: gate_pass["G3"] += 1
        if r["gates"]["G4_numerical"]: gate_pass["G4"] += 1
        if r["gates"]["G5_crossrefs"]: gate_pass["G5"] += 1
        if r["gates"]["G6_sm_anchors"]: gate_pass["G6"] += 1
        if r["appendices"]["cosmogenesis"]: appendix_present["cosmogenesis"] += 1
        if r["appendices"]["vds_dvp_bsh"]: appendix_present["vds_dvp_bsh"] += 1
        if r["critical"] or r["score"] < 70:
            needs_upgrade.append({
                "paper": r["paper"],
                "score": r["score"],
                "tier": r["tier"],
                "critical": r["critical"],
                "issues": r["issues"],
                "lines": r["lines"],
            })

    # Sort needs_upgrade by score ascending (worst first)
    needs_upgrade.sort(key=lambda x: x["score"])

    summary = {
        "audit_date": "2026-04-13",
        "range": f"PAPER_{RANGE_START:03d} to PAPER_{RANGE_END}",
        "total_scanned": total,
        "missing_files": len(missing_files),
        "missing_file_numbers": missing_files[:50],  # First 50
        "tiers": tiers,
        "gate_pass_rates": {k: f"{v}/{total} ({100*v/total:.0f}%)" for k, v in gate_pass.items()} if total > 0 else {},
        "appendix_rates": {k: f"{v}/{total} ({100*v/total:.0f}%)" for k, v in appendix_present.items()} if total > 0 else {},
        "avg_score": round(sum(r["score"] for r in results) / max(1, total), 1),
        "needs_upgrade_count": len(needs_upgrade),
        "needs_upgrade_papers": needs_upgrade,
    }

    output = {
        "summary": summary,
        "all_results": results,
    }

    out_path = os.path.join(os.path.dirname(__file__), "audit_gold_standard_results.json")
    with open(out_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2, ensure_ascii=False)

    # Print summary
    print(f"\n{'='*70}")
    print(f"  GOLD STANDARD AUDIT — PAPER_{RANGE_START:03d} to PAPER_{RANGE_END}")
    print(f"{'='*70}")
    print(f"  Scanned:        {total}")
    print(f"  Missing files:  {len(missing_files)}")
    print(f"  Average score:  {summary['avg_score']}/100")
    print(f"")
    print(f"  TIER DISTRIBUTION:")
    for tier, count in tiers.items():
        bar = '█' * (count // 5) if count > 0 else ''
        print(f"    {tier:10s} {count:5d}  {bar}")
    print(f"")
    print(f"  GATE PASS RATES:")
    for gate, rate in summary["gate_pass_rates"].items():
        print(f"    {gate}: {rate}")
    print(f"")
    print(f"  APPENDIX PRESENCE:")
    for app, rate in summary["appendix_rates"].items():
        print(f"    {app}: {rate}")
    print(f"")
    print(f"  NEEDS UPGRADE:  {len(needs_upgrade)} papers")
    if needs_upgrade:
        print(f"  Worst 20:")
        for p in needs_upgrade[:20]:
            crit = ', '.join(p['critical']) if p['critical'] else 'score<70'
            print(f"    {p['paper']:12s} score={p['score']:3d} tier={p['tier']:8s} lines={p['lines']:4d}  [{crit}]")
    print(f"\n  Results → audit_gold_standard_results.json")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
