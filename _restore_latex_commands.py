#!/usr/bin/env python3
"""
Reverse the damage from overly broad backslash-to-backtick conversion.
Restores `word` back to \word for known LaTeX commands.
"""
import glob, os, re

# All LaTeX commands that could have been corrupted
LATEX_COMMANDS = [
    # Greek letters
    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta',
    'iota', 'kappa', 'lambda', 'mu', 'nu', 'xi', 'pi', 'rho', 'sigma',
    'tau', 'upsilon', 'phi', 'chi', 'psi', 'omega',
    'Alpha', 'Beta', 'Gamma', 'Delta', 'Epsilon', 'Zeta', 'Eta', 'Theta',
    'Iota', 'Kappa', 'Lambda', 'Mu', 'Nu', 'Xi', 'Pi', 'Rho', 'Sigma',
    'Tau', 'Upsilon', 'Phi', 'Chi', 'Psi', 'Omega',
    'varepsilon', 'vartheta', 'varpi', 'varrho', 'varsigma', 'varphi',
    # Math operators
    'sin', 'cos', 'tan', 'cot', 'sec', 'csc',
    'arcsin', 'arccos', 'arctan',
    'sinh', 'cosh', 'tanh', 'coth',
    'log', 'ln', 'exp', 'det', 'dim', 'ker', 'hom',
    'deg', 'gcd', 'max', 'min', 'sup', 'inf', 'lim',
    'arg', 'Pr', 'mod',
    # Math symbols
    'sum', 'prod', 'int', 'oint', 'iint', 'iiint',
    'frac', 'dfrac', 'tfrac', 'binom', 'dbinom',
    'sqrt', 'root',
    'partial', 'nabla', 'infty', 'hbar',
    'cdot', 'cdots', 'ldots', 'vdots', 'ddots',
    'times', 'div', 'pm', 'mp', 'circ', 'bullet', 'star',
    'approx', 'equiv', 'neq', 'leq', 'geq', 'sim', 'simeq',
    'll', 'gg', 'subset', 'supset', 'subseteq', 'supseteq',
    'in', 'ni', 'notin', 'forall', 'exists', 'nexists',
    'wedge', 'vee', 'cap', 'cup',
    'oplus', 'otimes', 'ominus',
    'rightarrow', 'leftarrow', 'leftrightarrow',
    'Rightarrow', 'Leftarrow', 'Leftrightarrow',
    'implies', 'iff', 'to', 'mapsto',
    'perp', 'parallel', 'angle',
    # Formatting
    'text', 'textrm', 'textbf', 'textit', 'textsf', 'texttt',
    'mathrm', 'mathbf', 'mathit', 'mathsf', 'mathtt', 'mathcal', 'mathbb', 'mathfrak',
    'boldsymbol', 'bm',
    'hat', 'tilde', 'bar', 'vec', 'dot', 'ddot', 'check', 'breve', 'acute', 'grave',
    'overline', 'underline', 'overbrace', 'underbrace', 'widehat', 'widetilde',
    'overleftarrow', 'overrightarrow',
    # Delimiters
    'left', 'right', 'bigl', 'bigr', 'Bigl', 'Bigr', 'biggl', 'biggr', 'Biggl', 'Biggr',
    'langle', 'rangle', 'lfloor', 'rfloor', 'lceil', 'rceil', 'lvert', 'rvert',
    # Environments
    'begin', 'end',
    'boxed', 'cancel', 'bcancel', 'xcancel',
    # Spacing
    'quad', 'qquad', 'hspace', 'vspace', 'phantom', 'hphantom', 'vphantom',
    'kern', 'mkern', 'mskip',
    # Other common commands
    'label', 'ref', 'eqref', 'tag', 'notag', 'nonumber',
    'color', 'textcolor',
    'rule', 'strut',
    'limits', 'nolimits', 'displaystyle', 'textstyle', 'scriptstyle',
    'operatorname',
]

# Build a set for fast lookup
latex_set = set(LATEX_COMMANDS)

files = sorted(glob.glob('whitepapers/PAPER_*.md'))
fixed_count = 0

for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    original = c

    # Replace `word` with \word only for known LaTeX commands
    # But only when the backtick-wrapped word appears in math context
    # (near $, \, or other LaTeX)
    for cmd in LATEX_COMMANDS:
        backtick_cmd = f'`{cmd}`'
        if backtick_cmd in c:
            # Replace it back to \cmd
            c = c.replace(backtick_cmd, f'\\{cmd}')

    if c != original:
        with open(f, 'w', encoding='utf-8', newline='\n') as fh:
            fh.write(c)
        fixed_count += 1

print(f'Restored LaTeX commands in {fixed_count} papers')
