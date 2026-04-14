#!/usr/bin/env python3
"""Fix merged LaTeX Greek letter commands (e.g., \Deltaomega -> \Delta\omega)."""
import glob, os, re

# Greek letters (lowercase) that can follow uppercase Greek
LOWER_GREEK = [
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'pi', 'rho', 'sigma',
    'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega',
    'varepsilon', 'vartheta', 'varpi', 'varrho', 'varsigma', 'varphi',
]

UPPER_GREEK = [
    'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta',
    'Iota', 'Kappa', 'Lambda', 'Mu', 'Nu', 'Xi', 'Pi', 'Rho', 'Sigma',
    'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi', 'Omega',
]

# Also common LaTeX commands that may have merged
LATEX_FUNCTIONS = [
    'sin', 'cos', 'tan', 'exp', 'log', 'ln', 'sqrt', 'frac',
    'partial', 'nabla', 'cdot', 'times', 'approx',
    'text', 'mathrm', 'mathcal', 'mathbf', 'left', 'right',
    'bigl', 'bigr', 'sum', 'int', 'prod',
]

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Fix: \UpperGreek + lowercase_greek -> \UpperGreek\lowercase_greek
    for upper in UPPER_GREEK:
        for lower in LOWER_GREEK:
            merged = f'\\{upper}{lower}'
            fixed_str = f'\\{upper}\\{lower}'
            if merged in c:
                c = c.replace(merged, fixed_str)

        # Also fix \Uppercase + single lowercase letter (could be variable)
        # e.g., \Deltat -> \Delta t, \DeltaE -> \Delta E
        pattern = f'\\\\{upper}([a-zA-Z])(?![a-z])'
        def repl(m):
            return f'\\{upper} {m.group(1)}'
        c = re.sub(pattern, repl, c)

    # Fix merged latex commands: \fracpartial -> \frac\partial etc.
    for cmd1 in LATEX_FUNCTIONS:
        for cmd2 in LATEX_FUNCTIONS:
            if cmd1 != cmd2:
                merged = f'\\{cmd1}{cmd2}'
                if merged in c and f'\\{cmd1}\\{cmd2}' not in c:
                    c = c.replace(merged, f'\\{cmd1}\\{cmd2}')

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed += 1

print(f'Fixed merged Greek in {fixed} papers')
