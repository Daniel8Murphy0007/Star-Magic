"""whitepaper_audit.py — Comprehensive CVW G1-G6 + QS 1-5 audit across all whitepapers."""
import os, re, sys, json
from pathlib import Path

ROOT = Path(r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic")
WP_DIR = ROOT / "whitepapers"
PDF_DIR = ROOT / "pdf"

# ── Collect all PAPER_*.md files ─────────────────────────────────────────────
def collect_papers():
    papers = {}
    for d, label in [(WP_DIR, "whitepapers"), (ROOT, "root")]:
        if not d.exists():
            continue
        for f in d.iterdir():
            if f.is_file() and f.name.startswith("PAPER_") and f.name.endswith(".md"):
                m = re.match(r"PAPER_(\d+)", f.name)
                if m:
                    num = int(m.group(1))
                    if num not in papers:
                        papers[num] = []
                    papers[num].append({"file": f.name, "path": str(f), "loc": label})
    return papers

# ── Collect all PDFs ─────────────────────────────────────────────────────────
def collect_pdfs():
    pdfs = set()
    if PDF_DIR.exists():
        for f in PDF_DIR.iterdir():
            if f.is_file() and f.suffix == ".pdf":
                m = re.match(r"PAPER_(\d+)", f.name)
                if m:
                    pdfs.add(int(m.group(1)))
    return pdfs

# ── Gate checks on a single file ─────────────────────────────────────────────
def audit_file(filepath):
    try:
        text = Path(filepath).read_text(encoding="utf-8", errors="replace")
    except Exception as e:
        return {"error": str(e)}

    size_kb = len(text.encode("utf-8")) / 1024
    lines = text.split("\n")

    result = {"size_kb": round(size_kb, 1), "lines": len(lines)}

    # G1: Header — Author, Date, Session
    has_author = bool(re.search(r"(?i)\bauthors?\b", text))
    has_date = bool(re.search(r"(?i)\bdate\b.*\d{4}", text))
    has_session = bool(re.search(r"(?i)\bsession\b.*\d+", text))
    result["G1"] = has_author and has_date and has_session

    # G2: Abstract
    has_abstract = bool(re.search(r"(?i)##\s*(abstract|1\.\s*abstract)", text))
    abstract_text = ""
    if has_abstract:
        am = re.search(r"(?i)##\s*(?:abstract|1\.\s*abstract)\s*\n(.*?)(?=\n##|\Z)", text, re.DOTALL)
        if am:
            abstract_text = am.group(1).strip()
    result["G2"] = has_abstract and len(abstract_text.split()) >= 15

    # G3: Core Equation  (display math $$ or ``` code blocks with equations)
    has_display_math = bool(re.search(r"\$\$", text))
    has_code_eq = bool(re.search(r"```[\s\S]*?[=×·∑∫Σ][\s\S]*?```", text))
    has_inline_eq = len(re.findall(r"[=×·∑Σ∫]", text)) >= 3
    result["G3"] = has_display_math or has_code_eq or has_inline_eq

    # G4: Numerical Result  (scientific notation or decimal with units)
    sci_matches = re.findall(r"\d+\.?\d*\s*[×x]\s*10[\^⁻⁺⁰¹²³⁴⁵⁶⁷⁸⁹]+|\d+\.\d+e[+-]?\d+", text)
    result["G4"] = len(sci_matches) >= 1

    # G5: Anchor Cross-Reference (cites at least one PAPER_NNN)
    paper_refs = re.findall(r"PAPER_\d+", text)
    # Exclude self-reference
    own_num_match = re.match(r"PAPER_(\d+)", Path(filepath).name)
    own_num = int(own_num_match.group(1)) if own_num_match else -1
    external_refs = [r for r in paper_refs if int(re.search(r"\d+", r).group()) != own_num]
    result["G5"] = len(external_refs) >= 1
    result["paper_refs_count"] = len(external_refs)

    # G6: SM Anchor (SM comparison table or measured/alignment text)
    has_sm_section = bool(re.search(r"(?i)§?SM\s*Anchor|Standard\s*Model\s*Cross|SM\s*pred|alignment|measured.*PDG|arXiv.*\d{4}\.\d{4,5}", text))
    has_sm_table = bool(re.search(r"(?i)observable.*UQFF.*SM|observable.*prediction.*experiment", text))
    result["G6"] = has_sm_section or has_sm_table

    # QS dimensions
    # Q1: Novel physics claim
    result["Q1"] = bool(re.search(r"(?i)new\s*physics|novel|beyond.*standard|BSM|UQFF\s*(predict|explain|derive|unique|differ)", text))
    # Q2: >=2 display equations
    eq_count = len(re.findall(r"\$\$.*?\$\$", text, re.DOTALL)) + len(re.findall(r"```[^`]+[=×·].*?```", text, re.DOTALL))
    if eq_count < 2:
        eq_count += min(len(re.findall(r"^\s*[A-Za-z_]+\s*=\s*", text, re.MULTILINE)), 5)
    result["Q2"] = eq_count >= 2
    # Q3: Specific numerical result
    result["Q3"] = len(sci_matches) >= 1
    # Q4: SM/observational comparison
    result["Q4"] = result["G6"]
    # Q5: Testable prediction
    result["Q5"] = bool(re.search(r"(?i)testable|falsif|predict|measur.*future|observ.*test|experiment.*verify", text))

    return result

# ── MAIN ─────────────────────────────────────────────────────────────────────
def main():
    papers = collect_papers()
    pdfs = collect_pdfs()

    total = len(papers)
    dupes = {k: v for k, v in papers.items() if len(v) > 1}
    root_only = {k: v for k, v in papers.items() if all(e["loc"] == "root" for e in v)}

    print(f"Total unique paper numbers: {total}")
    print(f"Papers in whitepapers/: {sum(1 for v in papers.values() if any(e['loc']=='whitepapers' for e in v))}")
    print(f"Papers ONLY in root (unmigrated): {len(root_only)}")
    for k in sorted(root_only.keys()):
        print(f"  PAPER_{k:03d}: {root_only[k][0]['file']}")
    print(f"Papers with duplicate files: {len(dupes)}")
    for k in sorted(dupes.keys()):
        locs = [f"{e['loc']}/{e['file']}" for e in dupes[k]]
        print(f"  PAPER_{k:03d}: {', '.join(locs)}")
    print(f"Total PDFs in pdf/: {len(pdfs)}")

    # Audit each paper (prefer whitepapers/ copy)
    results = {}
    for num in sorted(papers.keys()):
        entries = papers[num]
        # Prefer whitepapers/ location
        entry = next((e for e in entries if e["loc"] == "whitepapers"), entries[0])
        audit = audit_file(entry["path"])
        results[num] = audit

    # Gate compliance summary
    gates = ["G1","G2","G3","G4","G5","G6"]
    qs = ["Q1","Q2","Q3","Q4","Q5"]
    
    gate_pass = {g: 0 for g in gates}
    gate_fail = {g: 0 for g in gates}
    qs_pass = {q: 0 for q in qs}
    qs_fail = {q: 0 for q in qs}
    
    for num, r in results.items():
        if "error" in r:
            continue
        for g in gates:
            if r.get(g):
                gate_pass[g] += 1
            else:
                gate_fail[g] += 1
        for q in qs:
            if r.get(q):
                qs_pass[q] += 1
            else:
                qs_fail[q] += 1

    audited = sum(1 for r in results.values() if "error" not in r)
    print(f"\n=== GATE COMPLIANCE ({audited} papers audited) ===")
    for g in gates:
        pct = gate_pass[g] / audited * 100 if audited else 0
        print(f"  {g}: {gate_pass[g]}/{audited} PASS ({pct:.1f}%) | {gate_fail[g]} FAIL")

    print(f"\n=== QS CONTENT DIMENSIONS ({audited} papers) ===")
    for q in qs:
        pct = qs_pass[q] / audited * 100 if audited else 0
        print(f"  {q}: {qs_pass[q]}/{audited} PASS ({pct:.1f}%) | {qs_fail[q]} FAIL")

    # Size distribution
    sizes = [r["size_kb"] for r in results.values() if "error" not in r]
    stubs = sum(1 for s in sizes if s < 5)
    standard = sum(1 for s in sizes if 5 <= s < 10)
    full = sum(1 for s in sizes if 10 <= s < 20)
    flagship = sum(1 for s in sizes if s >= 20)
    print(f"\n=== SIZE DISTRIBUTION ===")
    print(f"  Stub (<5KB):      {stubs}")
    print(f"  Standard (5-10KB): {standard}")
    print(f"  Full (10-20KB):    {full}")
    print(f"  Flagship (>20KB):  {flagship}")

    # PDF coverage
    papers_without_pdf = sorted(set(papers.keys()) - pdfs)
    pdfs_without_paper = sorted(pdfs - set(papers.keys()))
    print(f"\n=== PDF COVERAGE ===")
    print(f"  Papers with PDF:    {len(set(papers.keys()) & pdfs)}")
    print(f"  Papers WITHOUT PDF: {len(papers_without_pdf)}")
    if papers_without_pdf:
        print(f"    Missing: {papers_without_pdf[:30]}{'...' if len(papers_without_pdf)>30 else ''}")
    print(f"  PDFs without paper: {len(pdfs_without_paper)}")

    # Worst failures — papers failing 3+ gates
    multi_fail = []
    for num, r in sorted(results.items()):
        if "error" in r:
            continue
        fails = [g for g in gates if not r.get(g)]
        if len(fails) >= 3:
            multi_fail.append((num, fails, r["size_kb"]))
    print(f"\n=== PAPERS FAILING 3+ GATES ({len(multi_fail)}) ===")
    for num, fails, sz in multi_fail[:50]:
        print(f"  PAPER_{num:03d} ({sz:.0f}KB): FAIL {', '.join(fails)}")
    if len(multi_fail) > 50:
        print(f"  ... and {len(multi_fail)-50} more")

    # G6 failures breakdown (most common gap)
    g6_fail_list = sorted([num for num, r in results.items() if "error" not in r and not r.get("G6")])
    print(f"\n=== G6 (SM ANCHOR) FAILURES: {len(g6_fail_list)} papers ===")
    if g6_fail_list:
        # Show ranges
        ranges = []
        start = g6_fail_list[0]
        end = start
        for n in g6_fail_list[1:]:
            if n == end + 1:
                end = n
            else:
                ranges.append(f"{start}-{end}" if start != end else str(start))
                start = end = n
        ranges.append(f"{start}-{end}" if start != end else str(start))
        print(f"  Ranges: {', '.join(ranges[:20])}{'...' if len(ranges)>20 else ''}")

    # Write full results JSON
    output = {
        "audit_date": "2026-04-07",
        "total_papers": total,
        "total_pdfs": len(pdfs),
        "papers_without_pdf": papers_without_pdf,
        "root_only_papers": sorted(root_only.keys()),
        "duplicate_papers": sorted(dupes.keys()),
        "gate_summary": {g: {"pass": gate_pass[g], "fail": gate_fail[g]} for g in gates},
        "qs_summary": {q: {"pass": qs_pass[q], "fail": qs_fail[q]} for q in qs},
        "size_distribution": {"stub": stubs, "standard": standard, "full": full, "flagship": flagship},
        "multi_fail_papers": [{"paper": num, "fails": fails, "size_kb": sz} for num, fails, sz in multi_fail],
        "g6_fail_count": len(g6_fail_list),
    }
    with open(str(ROOT / "whitepaper_audit_results.json"), "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nFull results written to whitepaper_audit_results.json")

if __name__ == "__main__":
    main()
