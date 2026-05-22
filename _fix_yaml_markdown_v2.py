#!/usr/bin/env python3
"""
Fix YAML frontmatter in markdown files - more aggressive approach.

Remove all $ signs, backslashes in math contexts, and problematic
underscores in YAML title fields.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WHITEPAPERS_DIR = ROOT / "whitepapers"

def sanitize_title(title: str) -> str:
    """Aggressively remove LaTeX from title."""
    # Remove \kappa, \rho, \mu, etc. (LaTeX commands)
    title = re.sub(r'\\[a-zA-Z_]+', '', title)
    # Remove $ signs (math delimiters)
    title = re.sub(r'\$', '', title)
    # Remove _ followed by numbers/letters in math context
    title = re.sub(r'_([a-zA-Z0-9])', r'\1', title)
    # Remove \\( and \\)
    title = re.sub(r'\\\([^)]*\\\)', '', title)
    # Clean up extra spaces and dashes
    title = re.sub(r'\s+', ' ', title).strip()
    title = re.sub(r'\s*—\s*', ' - ', title)
    return title

def fix_markdown_file(filepath: Path) -> tuple[bool, str]:
    """Fix a single markdown file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find YAML frontmatter
    if not content.startswith('---'):
        return False, "No YAML"
    
    parts = content.split('---', 2)
    if len(parts) < 3:
        return False, "Malformed frontmatter"
    
    yaml_block = parts[1]
    body = '---'.join(parts[2:])
    
    # Find title line in YAML
    lines = yaml_block.split('\n')
    changed = False
    
    for i, line in enumerate(lines):
        if line.startswith('title:'):
            original = line
            # Extract the title value
            match = re.match(r'^(title:\s*)(.*?)$', line)
            if match:
                prefix, value = match.groups()
                # Remove quotes if present
                if value.startswith('"') and value.endswith('"'):
                    inner = value[1:-1]
                    clean = sanitize_title(inner)
                    lines[i] = f'{prefix}"{clean}"'
                elif value.startswith("'") and value.endswith("'"):
                    inner = value[1:-1]
                    clean = sanitize_title(inner)
                    lines[i] = f"{prefix}'{clean}'"
                else:
                    clean = sanitize_title(value)
                    lines[i] = f'{prefix}{clean}'
                
                if lines[i] != original:
                    changed = True
    
    if not changed:
        return False, "No changes"
    
    # Reconstruct and write
    new_yaml = '\n'.join(lines)
    new_content = f'---{new_yaml}---{body}'
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    return True, "Fixed"

def main():
    papers = sorted(WHITEPAPERS_DIR.glob("PAPER_*.md"))
    
    print(f"[Fix YAML v2] Processing {len(papers)} files...\n")
    
    fixed = 0
    failed = 0
    
    for paper in papers:
        changed, msg = fix_markdown_file(paper)
        if changed:
            fixed += 1
            print(f"  {paper.name}")
        elif "error" in msg.lower():
            failed += 1
    
    print(f"\n[Fix YAML v2] Fixed: {fixed}/{len(papers)}")

if __name__ == "__main__":
    main()
