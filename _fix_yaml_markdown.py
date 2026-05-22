#!/usr/bin/env python3
"""
Fix YAML frontmatter in markdown files with LaTeX expressions that break parsing.

Issue: Files have $ signs and underscores in YAML title fields, which breaks
pandoc YAML parsing. Solution: Remove LaTeX math from YAML title, move it to body.

Usage:
    python _fix_yaml_markdown.py
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WHITEPAPERS_DIR = ROOT / "whitepapers"

def strip_latex_from_title(title: str) -> str:
    """Remove LaTeX math expressions from title string."""
    # Remove $...$ (inline math)
    title = re.sub(r'\$[^$]*\$', '', title)
    # Remove \\(...\\) (LaTeX inline)
    title = re.sub(r'\\\([^)]*\\\)', '', title)
    # Clean up extra spaces
    title = re.sub(r'\s+', ' ', title).strip()
    return title

def fix_markdown_file(filepath: Path) -> tuple[bool, str]:
    """Fix a single markdown file. Returns (changed, message)."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except Exception as e:
        return False, f"Read error: {e}"
    
    # Check if file starts with YAML frontmatter
    if not lines or not lines[0].strip().startswith('---'):
        return False, "No YAML frontmatter"
    
    # Find end of frontmatter
    fm_end = -1
    for i in range(1, len(lines)):
        if lines[i].strip().startswith('---'):
            fm_end = i
            break
    
    if fm_end == -1:
        return False, "Malformed frontmatter"
    
    # Extract and fix frontmatter
    fm_lines = lines[:fm_end + 1]
    body_lines = lines[fm_end + 1:]
    
    changed = False
    for i in range(len(fm_lines)):
        line = fm_lines[i]
        if line.startswith('title:'):
            # Extract title value and clean it
            match = re.match(r'(title:\s*"?)(.+?)("?\s*)$', line)
            if match:
                prefix, title, suffix = match.groups()
                clean_title = strip_latex_from_title(title)
                if clean_title != title:
                    fm_lines[i] = f'{prefix}{clean_title}{suffix}\n'
                    changed = True
    
    if not changed:
        return False, "No LaTeX found in title"
    
    # Write back
    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.writelines(fm_lines + body_lines)
        return True, "Fixed"
    except Exception as e:
        return False, f"Write error: {e}"

def main():
    """Process all PAPER_*.md files."""
    papers = sorted(WHITEPAPERS_DIR.glob("PAPER_*.md"))
    
    if not papers:
        print("ERROR: No PAPER_*.md files found")
        return
    
    print(f"[Fix YAML] Found {len(papers)} markdown files")
    print(f"[Fix YAML] Scanning for LaTeX in YAML titles...\n")
    
    fixed = 0
    skipped = 0
    failed = 0
    
    for paper in papers:
        changed, msg = fix_markdown_file(paper)
        if changed:
            fixed += 1
            print(f"  FIXED: {paper.name}")
        elif "error" in msg.lower():
            failed += 1
            print(f"  ERROR: {paper.name} — {msg}")
        else:
            skipped += 1
    
    print(f"\n[Fix YAML] Summary")
    print(f"  Fixed: {fixed}")
    print(f"  Skipped (no LaTeX): {skipped}")
    print(f"  Failed: {failed}")
    print(f"  Total: {len(papers)}")

if __name__ == "__main__":
    main()
