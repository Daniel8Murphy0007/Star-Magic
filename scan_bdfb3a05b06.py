"""Scan grok_share_bdfb3a05b06.txt for structure and unique physics blocks."""
import re

FNAME = 'grok_share_bdfb3a05b06.txt'

with open(FNAME, encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

print(f'Total lines: {len(lines)}')

# Find all fence markers
all_fences = []
for i, line in enumerate(lines):
    s = line.strip()
    if s.startswith('```'):
        all_fences.append((i+1, s))  # 1-based line number

print(f'\nTotal fence markers: {len(all_fences)}')
print('\nAll fences:')
for ln, s in all_fences:
    print(f'  L{ln}: {repr(s)}')

# Find class definitions
print('\n\nClass definitions:')
for i, line in enumerate(lines):
    if re.search(r'^class\s+\w+', line.strip()):
        print(f'  L{i+1}: {line.strip()[:120]}')

# Find Thoughts markers
print('\n\nThoughts markers:')
for i, line in enumerate(lines):
    if 'Thoughts' in line and len(line.strip()) < 30:
        print(f'  L{i+1}: {line.strip()[:120]}')

# Find document attachment markers (File label lines)
print('\n\nAttachment/File markers:')
for i, line in enumerate(lines):
    s = line.strip()
    if s.endswith('.docx\nFile') or (s.endswith('.docx') and i+1 < len(lines) and lines[i+1].strip() == 'File'):
        print(f'  L{i+1}: {s[:120]}')
    elif s == 'File' and i > 0 and lines[i-1].strip().endswith('.docx'):
        print(f'  L{i+1}: FILE marker after {lines[i-1].strip()[:80]}')

# Sample first 100 non-empty lines for context
print('\n\nFirst 80 lines (non-chrome):')
for i, line in enumerate(lines[:80]):
    print(f'L{i+1}: {line.rstrip()[:120]}')
