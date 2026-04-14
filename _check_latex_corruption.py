#!/usr/bin/env python3
"""Check if LaTeX math was corrupted by backslash-to-backtick conversion."""
import glob, re

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
problems = 0
for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    # Find $$ blocks with backticks inside
    for m in re.finditer(r'\$\$(.*?)\$\$', c, re.DOTALL):
        block = m.group(1)
        if '`' in block:
            problems += 1
            if problems <= 5:
                print(f'PROBLEM in {f}: {block[:150]}')
            break

print(f'Total files with backtick-in-math: {problems}')

# Check for `sin`, `cos`, `theta` etc that were LaTeX commands
latex_cmds = ['sin', 'cos', 'tan', 'log', 'exp', 'ln', 'sum', 'int', 'frac',
              'theta', 'alpha', 'beta', 'gamma', 'delta', 'lambda', 'mu', 'nu',
              'sigma', 'omega', 'phi', 'psi', 'rho', 'epsilon', 'kappa',
              'text', 'mathrm', 'mathcal', 'mathbf', 'left', 'right',
              'bigl', 'bigr', 'cdot', 'times', 'approx', 'equiv', 'neq',
              'partial', 'nabla', 'infty', 'sqrt', 'vec', 'hat', 'bar',
              'boxed', 'begin', 'end']
corrupted = 0
for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    for cmd in latex_cmds:
        if f'`{cmd}`' in c:
            corrupted += 1
            if corrupted <= 5:
                print(f'CORRUPTED LaTeX cmd `{cmd}` in {f}')
            break

print(f'Files with corrupted LaTeX commands: {corrupted}')
