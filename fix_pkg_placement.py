#!/usr/bin/env python3
"""
fix_pkg_placement.py — Strip incorrectly placed PKG blocks from file tops.

The integrate_body_physics.py script had a bug where PKG blocks were inserted
at position 0 (before YAML frontmatter) instead of between body and appendix.
This script strips those blocks so the fixed script can re-insert them correctly.
"""

import os
import re
import glob

REPO = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(REPO, "whitepapers")

def strip_bad_pkg(filepath):
    """Remove Session 225 PKG blocks that were incorrectly placed at file top."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
    
    # Only process files that have PKG markers
    if 'PKG-' not in content or '<!-- PKG-' not in content:
        return False
    
    # Pattern: file starts with optional whitespace/newlines, then ---,
    # then ## Session 225: Late-Corpus Physics Integration,
    # then PKG blocks, then the real YAML frontmatter (---\npaper_id:)
    #
    # We need to find where the real YAML frontmatter starts.
    # The real frontmatter has: ---\npaper_id: PAPER_xxx
    
    # Check if the PKG blocks are at the top (before YAML frontmatter)
    # Real YAML frontmatter pattern: ---\npaper_id: or ---\ntitle:
    yaml_start = re.search(r'^---\s*\npaper_id:', content, re.MULTILINE)
    if not yaml_start:
        yaml_start = re.search(r'^---\s*\ntitle:', content, re.MULTILINE)
    
    if not yaml_start:
        return False  # Can't find YAML frontmatter
    
    yaml_pos = yaml_start.start()
    
    if yaml_pos == 0:
        # YAML is at position 0 — PKG is NOT at top, it's somewhere else
        # Check if PKG is correctly placed (not at top)
        return False
    
    # PKG blocks are before position yaml_pos — strip them
    # The correct content starts at yaml_pos
    prefix = content[:yaml_pos]
    
    # Verify the prefix actually has Session 225 PKG content
    if '## Session 225: Late-Corpus Physics Integration' not in prefix:
        return False
    
    # Strip the bad prefix
    new_content = content[yaml_pos:]
    
    with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
        f.write(new_content)
    
    return True


def main():
    # Collect all paper files
    paper_files = sorted(glob.glob(os.path.join(WP_DIR, "PAPER_*.md")))
    root_papers = sorted(glob.glob(os.path.join(REPO, "PAPER_*.md")))
    all_files = paper_files + root_papers
    
    print(f"Scanning {len(all_files)} paper files for bad PKG placement...")
    
    stripped = 0
    for fp in all_files:
        fname = os.path.basename(fp)
        if strip_bad_pkg(fp):
            stripped += 1
            print(f"  [FIXED] {fname}")
    
    print(f"\nStripped bad PKG blocks from {stripped} files.")
    print("Now re-run integrate_body_physics.py with the fixed insertion logic.")


if __name__ == "__main__":
    main()
