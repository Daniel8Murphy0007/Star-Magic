#!/usr/bin/env python3
"""Fix broken SM Anchors LaTeX in 68 papers."""
import glob, os, re

CORRECT_TABLE = """| Observable | UQFF Prediction | SM / Experiment | Source | Alignment |
|------------|-----------------|-----------------|--------|-----------|
| $\\sin^2\\theta_W$ | Embedded in $U_{g2}$ charge coupling | $0.2312$ | PDG 2024 | 99.6% |
| Fine structure $\\alpha$ | UQFF reproduces via $U_{g1}$ dipole | $1/137.036$ | PDG 2024 | 99.9% |
| $m_Z$ | SCm phonon predicts $Z$ mass | $91.1876$ GeV | PDG 2024 | 99.8% |"""

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0
for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    if '\\$\\\\sin' in c or '\\$\\\\alpha' in c:
        # Find the broken table section
        lines = c.split('\n')
        new_lines = []
        skip = False
        i = 0
        while i < len(lines):
            line = lines[i]
            # Detect start of broken table
            if '| Observable | UQFF Prediction' in line:
                # Check if next lines have the broken pattern
                ahead = '\n'.join(lines[i:i+6])
                if '\\$\\\\sin' in ahead or '\\$\\\\alpha' in ahead:
                    # Replace with correct table
                    new_lines.append(CORRECT_TABLE)
                    # Skip the broken rows (header + separator + 3 data rows)
                    j = i + 1
                    while j < len(lines) and (lines[j].startswith('|') or lines[j].strip() == ''):
                        if not lines[j].startswith('|'):
                            break
                        j += 1
                    i = j
                    fixed += 1
                    continue
            new_lines.append(line)
            i += 1

        result = '\n'.join(new_lines)
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(result)

print(f'Fixed SM Anchors in {fixed} papers')
