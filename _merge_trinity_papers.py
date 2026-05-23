#!/usr/bin/env python3
"""
_merge_trinity_papers.py — Merge extracted papers into reference documentation

Reads 1233 extracted papers from whitepapers/ directory and merges them
into _uqff_reference_documentation.py registry.

High-impact task: 1222 new papers need to be added to existing 11.
"""

from pathlib import Path
import re
from typing import Dict, List


def extract_all_papers() -> List[Dict]:
    """Extract all papers from whitepapers/ directory."""
    papers = []
    whitepaper_dir = Path("whitepapers")
    
    if not whitepaper_dir.exists():
        print(f"ERROR: {whitepaper_dir} not found")
        return papers
    
    print(f"Scanning {whitepaper_dir}...")
    md_files = sorted(whitepaper_dir.glob("PAPER_*.md"))
    print(f"  Found {len(md_files)} files")
    
    for md_file in md_files:
        try:
            # Extract paper number
            match = re.search(r'PAPER_(\d+)', md_file.name)
            if not match:
                continue
            
            paper_num = int(match.group(1))
            paper_id = f"PAPER_{paper_num}"
            
            # Read metadata
            with open(md_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read(2000)  # First 2KB for metadata
            
            # Extract title
            title_match = re.search(r'^#\s+(.+?)$', content, re.MULTILINE)
            title = title_match.group(1) if title_match else md_file.stem
            
            # Extract domain
            domain = infer_domain_from_filename(paper_id)
            session = infer_session_from_paper_num(paper_num)
            
            paper = {
                'paper_id': paper_id,
                'title': title,
                'domain': domain,
                'session_number': session,
                'md_path': str(md_file.relative_to(Path.cwd())),
                'pdf_path': f"pdf/PAPER_{paper_num:03d}.pdf" if paper_num < 1000 else f"pdf/PAPER_{paper_num}.pdf",
            }
            
            papers.append(paper)
        
        except Exception as e:
            print(f"  WARNING: Error reading {md_file.name}: {e}")
            continue
    
    return papers


def infer_domain_from_filename(paper_id: str) -> str:
    """Infer domain from PAPER_XXX number."""
    match = re.search(r'PAPER_(\d+)', paper_id)
    if not match:
        return "Miscellaneous"
    
    num = int(match.group(1))
    
    if 237 <= num <= 241:
        return "UQFF Primitives & Foundations"
    elif 242 <= num <= 340:
        return "Quantum Mechanics & QED"
    elif 341 <= num <= 439:
        return "Cosmological Constants & Structure"
    elif 440 <= num <= 539:
        return "Mathematics & Special Constants"
    elif 540 <= num <= 665:
        return "Astrophysics & Stellar Physics"
    elif 725 <= num <= 785:
        return "Calibration & Parameter Tuning"
    else:
        return "Cross-Domain Integration"


def infer_session_from_paper_num(paper_num: int) -> int:
    """Infer session from PAPER_XXX number."""
    # Simple approximation
    if paper_num <= 241:
        return 237 + (paper_num - 237)
    else:
        return paper_num


def generate_paper_code(paper: Dict) -> str:
    """Generate Python code for a single paper."""
    paper_id = paper['paper_id']
    title = paper['title'].replace('"', r'\"')
    domain = paper['domain']
    session = paper['session_number']
    md_path = paper['md_path'].replace('\\', '/')
    pdf_path = paper['pdf_path'].replace('\\', '/')
    
    code = f'''
    {paper_id} = ReferenceDocument(
        paper_id="{paper_id}",
        title="{title}",
        domain=PaperDomain.CROSS_DOMAIN,  # TODO: Set correct domain from enum
        session_number={session},
        status=PaperStatus.PUBLISHED,
        year=2026,
        authors=[],
        abstract="Physics research paper from UQFF framework",
        constants_defined=[],
        related_papers=[],
        pdf_path="{pdf_path}",
        md_path="{md_path}",
        cvw_version="v2.0.0",
        error_percent=0.0,
        page_count=0,
        published_date="2026-01-01",
        keywords=[]
    )
'''
    return code


def group_papers_by_domain(papers: List[Dict]) -> Dict[str, List[Dict]]:
    """Group papers by domain."""
    grouped = {}
    for paper in papers:
        domain = paper['domain']
        if domain not in grouped:
            grouped[domain] = []
        grouped[domain].append(paper)
    return grouped


def generate_registry_code(papers: List[Dict]) -> str:
    """Generate registry insertion code."""
    code = ""
    
    # Group by session range
    ranges = {
        'SESSIONS_203_340': [],
        'SESSIONS_341_439': [],
        'SESSIONS_440_539': [],
        'SESSIONS_540_665': [],
        'SESSIONS_725_785': [],
    }
    
    for paper in papers:
        session = paper['session_number']
        paper_id = paper['paper_id']
        
        if 203 <= session <= 340:
            ranges['SESSIONS_203_340'].append(paper_id)
        elif 341 <= session <= 439:
            ranges['SESSIONS_341_439'].append(paper_id)
        elif 440 <= session <= 539:
            ranges['SESSIONS_440_539'].append(paper_id)
        elif 540 <= session <= 665:
            ranges['SESSIONS_540_665'].append(paper_id)
        elif 725 <= session <= 785:
            ranges['SESSIONS_725_785'].append(paper_id)
    
    for range_name, paper_ids in ranges.items():
        if paper_ids:
            code += f"\n# Add to {range_name}:\n"
            for pid in paper_ids[:10]:  # Show first 10
                code += f'    "{pid}": {pid},\n'
            if len(paper_ids) > 10:
                code += f'    # ... and {len(paper_ids)-10} more\n'
    
    return code


def print_merge_summary(papers: List[Dict]) -> None:
    """Print summary of papers to merge."""
    print("\n" + "="*80)
    print("REFERENCE PAPER MERGE SUMMARY")
    print("="*80 + "\n")
    
    grouped = group_papers_by_domain(papers)
    
    for domain in sorted(grouped.keys()):
        count = len(grouped[domain])
        print(f"{domain:<40} {count:4d} papers")
    
    print(f"\n{'='*80}")
    print(f"TOTAL PAPERS TO MERGE: {len(papers)}")
    print(f"{'='*80}\n")


def main():
    print("\n" + "="*80)
    print("TRINITY PAPER MERGE SCRIPT")
    print("="*80 + "\n")
    
    # Extract papers
    papers = extract_all_papers()
    
    # Show summary
    print_merge_summary(papers)
    
    # Generate sample code
    print("Sample Paper Definition Code (First Paper):\n")
    if papers:
        print(generate_paper_code(papers[0]))
    
    # Generate registry code
    print("\nRegistry Insertion Code (Sample):\n")
    print(generate_registry_code(papers[:50]))  # Show first 50
    
    print("\n" + "="*80)
    print("ACTION ITEMS:")
    print("="*80)
    print(f"1. Update enum PaperDomain to handle all {len(papers)} papers")
    print(f"2. Create registry dictionaries for 5 session ranges")
    print(f"3. Bulk-insert {len(papers)} paper definitions with correct domains")
    print(f"4. Generate cross-reference indices (papers → constants)")
    print("\n")


if __name__ == '__main__':
    main()
