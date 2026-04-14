#!/usr/bin/env python3
"""Fix PAPER_1013-1018: prose paragraphs wrongly wrapped in $$ LaTeX."""
import glob, os

for num in [1013, 1014, 1015, 1016, 1017, 1018]:
    pattern = f'whitepapers/PAPER_{num}_*.md'
    matches = glob.glob(pattern)
    for f in matches:
        c = open(f, encoding='utf-8').read()
        lines = c.split('\n')
        new_lines = []
        fixed = False
        for line in lines:
            if line.startswith('$$') and len(line) > 20:
                inner = line[2:]
                if inner.endswith('$$'):
                    inner = inner[:-2]
                prose_words = ['the ', 'and ', 'for ', 'with ', 'from ',
                               'across ', 'apply ', 'that ', 'this ',
                               'which ', 'where ', 'using ', 'between ']
                word_count = sum(1 for w in prose_words if w.lower() in inner.lower())
                if word_count >= 2:
                    new_lines.append(inner)
                    fixed = True
                    continue
            new_lines.append(line)
        if fixed:
            with open(f, 'w', encoding='utf-8', newline='\n') as fh:
                fh.write('\n'.join(new_lines))
            print(f'Fixed prose: {os.path.basename(f)}')

print('Done')
