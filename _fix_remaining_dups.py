#!/usr/bin/env python3
"""Fix remaining 7 papers with duplicate author/header blocks."""
import glob, os, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    lines = c.split('\n')
    auth_lines = [i for i, l in enumerate(lines) if l.startswith('**Author:**')]
    
    if len(auth_lines) <= 1:
        continue
    
    # Keep the first author block, remove subsequent duplicate blocks
    # A "block" starts from the **Title:** or heading line before **Author:**
    # and ends at the next --- or ## or blank section
    
    first_auth = auth_lines[0]
    blocks_to_remove = []
    
    for auth_idx in auth_lines[1:]:
        # Find start of this duplicate block (scan backwards)
        block_start = auth_idx
        for j in range(auth_idx - 1, first_auth, -1):
            stripped = lines[j].strip()
            if stripped.startswith('**Title:') or stripped.startswith('**Session:') or stripped.startswith('## ') or stripped.startswith('# '):
                block_start = j
                break
            elif stripped == '' or stripped == '---':
                block_start = j
                continue
            elif stripped.startswith('**'):
                continue
            else:
                # Non-header text — this might be a title continuation
                if j > 0 and any(lines[k].strip().startswith('**Title:') or lines[k].strip().startswith('## ') for k in range(max(0, j-3), j)):
                    continue
                block_start = j + 1
                break
        
        # Find end of this duplicate block
        block_end = auth_idx
        for j in range(auth_idx + 1, min(len(lines), auth_idx + 20)):
            stripped = lines[j].strip()
            if stripped.startswith('## ') or stripped == '---':
                block_end = j - 1
                break
            elif stripped.startswith('**') or stripped == '':
                block_end = j
            else:
                block_end = j - 1
                break
        
        blocks_to_remove.append((block_start, block_end))
    
    # Remove blocks in reverse order to maintain indices
    for start, end in reversed(blocks_to_remove):
        del lines[start:end + 1]
    
    # Clean up multiple blank lines
    result = '\n'.join(lines)
    result = re.sub(r'\n{3,}', '\n\n', result)
    
    with open(f, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write(result)
    fixed += 1
    print(f'Fixed: {os.path.basename(f)} (removed {len(blocks_to_remove)} duplicate blocks)')

print(f'\nTotal fixed: {fixed}')
