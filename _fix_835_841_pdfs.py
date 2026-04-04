#!/usr/bin/env python3
"""
_fix_835_841_pdfs.py
Fix the 4 failing whitepaper .mds (836, 837, 839, 841) that produce
LaTeX "Missing $ inserted" errors. Root cause: newunicodechar in
pdf_header.tex makes Unicode sub/superscripts active inside verbatim.

Fix strategy:
  1. In code blocks: replace Unicode sub/super chars with ASCII
  2. Outside code blocks: wrap bare N^M patterns in $...$
"""
import os
import re

WHITEPAPER_DIR = "whitepapers"
FENCE = '```'

# Pattern: number^number (scientific notation) outside of code blocks
_CARET_NUM = re.compile(r'(\d+)\^(\{[^}]+\}|[\-]?\d+)')

# Unicode -> ASCII replacements for code blocks
UNICODE_TO_ASCII = {
    '²': '2', '³': '3', '¹': '1',
    '⁰': '0', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
    '⁻': '-', '⁺': '+',
    '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
    '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    'ᵢ': 'i', 'ₙ': 'n',
    '∫': 'integral', '×': '*', '÷': '/', '±': '+-',
    '≈': '~', '≥': '>=', '≤': '<=', '≠': '!=',
    'θ': 'theta', 'ρ': 'rho', 'π': 'pi', 'σ': 'sigma',
    'ω': 'omega', 'α': 'alpha', 'β': 'beta', 'γ': 'gamma',
    'δ': 'delta', 'ε': 'epsilon', 'λ': 'lambda', 'μ': 'mu',
    'ν': 'nu', 'Ω': 'Omega', 'Λ': 'Lambda', 'Σ': 'Sigma',
    'Δ': 'Delta', 'Γ': 'Gamma', 'η': 'eta',
    '∞': 'inf', '∂': 'd',
    '\u2013': '-', '\u2014': '--',  # en/em dash
    '\u2018': "'", '\u2019': "'",   # smart quotes
    '\u201c': '"', '\u201d': '"',
}


def sanitize_for_code(line):
    """Replace Unicode chars that newunicodechar makes active."""
    for uc, asc in UNICODE_TO_ASCII.items():
        line = line.replace(uc, asc)
    return line


def fix_markdown_for_latex(text):
    """
    1. Inside code blocks: replace Unicode sub/super/Greek with ASCII.
    2. Outside code blocks: wrap N^M patterns in $...$ for LaTeX math.
    """
    lines = text.split('\n')
    result = []
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(FENCE):
            in_block = not in_block
            result.append(line)  # keep fences as-is
            continue
        if in_block:
            result.append(sanitize_for_code(line))
        else:
            # Fix caret notation outside code: 10^208 -> $10^{208}$
            def _fix_caret(m):
                base = m.group(1)
                exp = m.group(2)
                if exp.startswith('{'):
                    return f'${base}^{exp}$'
                else:
                    return f'${base}^{{{exp}}}$'
            line = _CARET_NUM.sub(_fix_caret, line)
            result.append(line)
    return '\n'.join(result)


TARGETS = [
    'PAPER_836_Chandra_35System_UQFF_Survey_NegativeBuoyancy.md',
    'PAPER_837_Fquark_Fneutrino_FALP_Fdark_ArXiv_Bridge_UQFF.md',
    'PAPER_839_ADD_Large_Extra_Dimensions_FLED_UQFF.md',
    'PAPER_841_UQFF_Millennium_Prize_Applications.md',
]

for fname in TARGETS:
    fpath = os.path.join(WHITEPAPER_DIR, fname)
    with open(fpath, encoding='utf-8', errors='replace') as f:
        original = f.read()
    fixed = fix_markdown_for_latex(original)
    if fixed != original:
        with open(fpath, 'w', encoding='utf-8') as f:
            f.write(fixed)
        count = sum(1 for a, b in zip(original.split('\n'), fixed.split('\n')) if a != b)
        print(f"Fixed: {fname} ({count} lines changed)")
    else:
        print(f"No change: {fname}")

print("\nDone. Re-run: python generate_pdfs.py 835 841")
