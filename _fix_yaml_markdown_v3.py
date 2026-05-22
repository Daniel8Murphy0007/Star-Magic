#!/usr/bin/env python3
"""
Fix YAML frontmatter - v3: Rebuild YAML completely to ensure single-line format.
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
WHITEPAPERS_DIR = ROOT / "whitepapers"

def sanitize_title_aggressive(title: str) -> str:
    """Remove ALL special LaTeX characters."""
    # Remove \commands
    title = re.sub(r'\\[a-zA-Z]+\{[^}]*\}', '', title)
    title = re.sub(r'\\[a-zA-Z_]+', '', title)
    # Remove $ (math markers)
    title = re.sub(r'\$', '', title)
    # Remove { }
    title = re.sub(r'[{}]', '', title)
    # Fix underscores - keep them but surrounded by text
    title = re.sub(r'_+', '_', title)
    # Remove multiple spaces
    title = re.sub(r'\s{2,}', ' ', title)
    title = title.strip()
    return title

def fix_markdown_file_v3(filepath: Path) -> tuple[bool, str]:
    """Fix by rebuilding entire YAML structure."""
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Find frontmatter boundaries
    if not lines or not lines[0].startswith('---'):
        return False, "No FM"
    
    end_fm = -1
    for i in range(1, len(lines)):
        if lines[i].startswith('---'):
            end_fm = i
            break
    
    if end_fm == -1:
        return False, "Bad FM"
    
    # Parse YAML lines
    yaml_lines = lines[1:end_fm]
    body_lines = lines[end_fm + 1:]
    
    # Rebuild YAML with sanitized title
    new_yaml = []
    title_found = False
    title_value = None
    
    for line in yaml_lines:
        if line.startswith('title:'):
            title_found = True
            # Extract title value (may be quoted or unquoted)
            match = re.match(r'^title:\s*["\']?(.+?)["\']?\s*$', line)
            if match:
                raw_title = match.group(1).strip('"\'')
                title_value = sanitize_title_aggressive(raw_title)
                # Write title on single line with quotes
                new_yaml.append(f'title: "{title_value}"\n')
            else:
                new_yaml.append(line)
        else:
            # Keep other YAML lines as-is
            if line.strip():  # Non-empty
                new_yaml.append(line)
    
    if not title_found:
        return False, "No title"
    
    # Reconstruct file
    new_content = '---\n' + ''.join(new_yaml) + '---\n' + ''.join(body_lines)
    
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    return True, f"Fixed: {title_value[:50]}"

def main():
    papers = sorted(WHITEPAPERS_DIR.glob("PAPER_*.md"))
    
    print(f"[Fix YAML v3] Processing {len(papers)} files...\n")
    
    fixed = 0
    failed = 0
    
    for idx, paper in enumerate(papers):
        changed, msg = fix_markdown_file_v3(paper)
        if changed:
            fixed += 1
            if (idx + 1) % 100 == 0:
                print(f"  [{idx + 1}/1216] Fixed {fixed} so far...")
        elif "error" in msg.lower():
            failed += 1
    
    print(f"\n[Fix YAML v3] Summary")
    print(f"  Fixed: {fixed}/{len(papers)}")
    print(f"  Failed: {failed}")

if __name__ == "__main__":
    main()
