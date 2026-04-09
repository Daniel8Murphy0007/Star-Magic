#!/usr/bin/env python3
"""
fix_pdf_failures.py - Session 204
Fix the 35 whitepapers that fail PDF generation due to LaTeX/Unicode issues.

Issues fixed:
  1. Unicode ½ in text/code → 1/2
  2. Unicode math symbols inside code blocks (verbatim) that xelatex can't render
  3. Double superscripts in LaTeX equations (A_\mu^a^ → A_\mu^{a})
  4. Other problematic Unicode in non-math context

Strategy: Operate on the raw markdown. Two passes:
  A) Inside code blocks (``` fences): replace Unicode with ASCII equivalents
  B) Outside code blocks: replace isolated math Unicode with inline $...$-wrapped LaTeX
"""

import os, re, glob

WHITEPAPER_DIR = "whitepapers"

# The 35 failing papers
FAILING_PAPERS = [
    11, 14, 15, 19, 101, 102, 103, 104, 121, 172,
    183, 189, 190, 205, 246, 348, 376, 386, 476, 624,
    625, 649, 655, 741, 742, 743, 744, 745, 747, 748,
    750, 762, 763, 765, 841,
]

# Unicode → ASCII replacements INSIDE code blocks (verbatim environments)
CODE_BLOCK_REPLACEMENTS = {
    '½': '1/2',
    '⅓': '1/3',
    '⅔': '2/3',
    '¼': '1/4',
    '¾': '3/4',
    '∫': 'integral',
    '∇': 'nabla',
    '∂': 'd',
    '∝': '~',
    '∞': 'inf',
    '≈': '~=',
    '≤': '<=',
    '≥': '>=',
    '≠': '!=',
    '±': '+/-',
    '×': 'x',
    '·': '*',
    '÷': '/',
    '→': '->',
    '←': '<-',
    '↔': '<->',
    '⇒': '=>',
    '⇐': '<=',
    'ℏ': 'hbar',
    'ħ': 'hbar',
    '−': '-',
    '—': '--',
    '–': '-',
    'α': 'alpha',
    'β': 'beta',
    'γ': 'gamma',
    'δ': 'delta',
    'ε': 'epsilon',
    'ζ': 'zeta',
    'η': 'eta',
    'θ': 'theta',
    'ι': 'iota',
    'κ': 'kappa',
    'λ': 'lambda',
    'μ': 'mu',
    'ν': 'nu',
    'ξ': 'xi',
    'π': 'pi',
    'ρ': 'rho',
    'σ': 'sigma',
    'τ': 'tau',
    'υ': 'upsilon',
    'φ': 'phi',
    'χ': 'chi',
    'ψ': 'psi',
    'ω': 'omega',
    'Γ': 'Gamma',
    'Δ': 'Delta',
    'Θ': 'Theta',
    'Λ': 'Lambda',
    'Ξ': 'Xi',
    'Π': 'Pi',
    'Σ': 'Sigma',
    'Φ': 'Phi',
    'Ψ': 'Psi',
    'Ω': 'Omega',
    '²': '^2',
    '³': '^3',
    '⁴': '^4',
    '⁵': '^5',
    '⁶': '^6',
    '⁷': '^7',
    '⁸': '^8',
    '⁹': '^9',
    '⁰': '^0',
    '¹': '^1',
    '⁻': '^-',
    '₀': '_0',
    '₁': '_1',
    '₂': '_2',
    '₃': '_3',
    '₄': '_4',
    '₅': '_5',
    '₆': '_6',
    '₇': '_7',
    '₈': '_8',
    '₉': '_9',
    'ℓ': 'l',
    'χ̇': "chi'",
    'ẋ': "x'",
}

# Unicode → safe-text replacements OUTSIDE code blocks (in regular markdown text/tables)
# These get replaced with ASCII equivalents that pandoc/xelatex can handle
TEXT_REPLACEMENTS = {
    '½': '1/2',
    '⅓': '1/3',
    '⅔': '2/3',
    '¼': '1/4',
    '¾': '3/4',
    '∫': 'integral ',
    '∝': '~ ',
    'ℓ': 'l',
    'ℏ': 'hbar',
    'ħ': 'hbar',
    '∞': 'inf',
    '∇': 'nabla',
    '∂': 'd',
    '²': '^2',
    '³': '^3',
    '⁴': '^4',
    '⁻': '^-',
    '⁰': '^0',
    '¹': '^1',
    '⁵': '^5',
    '⁶': '^6',
    '⁷': '^7',
    '⁸': '^8',
    '⁹': '^9',
    '₀': '_0',
    '₁': '_1',
    '₂': '_2',
    '₃': '_3',
    '₄': '_4',
    '₅': '_5',
    '₆': '_6',
    '₇': '_7',
    '₈': '_8',
    '₉': '_9',
    '−': '-',
    '×': 'x',
    '·': '*',
    '→': '->',
    '←': '<-',
    '⇒': '=>',
    '—': '--',
    '–': '-',
}


def find_file_for_paper(paper_num):
    """Find the whitepaper .md file for a given paper number."""
    pattern = os.path.join(WHITEPAPER_DIR, f"PAPER_{paper_num:03d}_*.md")
    matches = glob.glob(pattern)
    if not matches:
        # Try 4-digit
        pattern = os.path.join(WHITEPAPER_DIR, f"PAPER_{paper_num:04d}_*.md")
        matches = glob.glob(pattern)
    return matches[0] if matches else None


def fix_code_blocks(text):
    """Replace Unicode chars inside markdown code blocks with ASCII equivalents."""
    lines = text.split('\n')
    result = []
    in_code_block = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('```'):
            in_code_block = not in_code_block
            result.append(line)
            continue

        if in_code_block:
            for uc, ascii_rep in CODE_BLOCK_REPLACEMENTS.items():
                if uc in line:
                    line = line.replace(uc, ascii_rep)
            result.append(line)
        else:
            result.append(line)

    return '\n'.join(result)


def fix_text_unicode(text):
    """Fix Unicode fractions and symbols in regular text (outside code blocks and $...$)."""
    lines = text.split('\n')
    result = []
    in_code_block = False

    for line in lines:
        stripped = line.strip()
        if stripped.startswith('```'):
            in_code_block = not in_code_block
            result.append(line)
            continue

        if not in_code_block:
            for uc, rep in TEXT_REPLACEMENTS.items():
                if uc in line:
                    line = line.replace(uc, rep)
            result.append(line)
        else:
            result.append(line)

    return '\n'.join(result)


def fix_double_superscript(text):
    """Fix LaTeX double superscript errors like A_\\mu^a^ → A_\\mu^{a}."""
    # Fix patterns like ^a^ (double superscript)
    # Match ^X^ where X is a single char (not { or })
    text = re.sub(r'\^([^{}\s\^])\^', r'^{\1}', text)
    return text


def fix_paper(filepath):
    """Apply all fixes to a single paper. Returns (changed, description)."""
    with open(filepath, 'r', encoding='utf-8') as f:
        original = f.read()

    text = original

    # Pass 1: Fix code blocks
    text = fix_code_blocks(text)

    # Pass 2: Fix text Unicode
    text = fix_text_unicode(text)

    # Pass 3: Fix double superscripts
    text = fix_double_superscript(text)

    if text != original:
        with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)
        return True
    return False


def main():
    fixed = 0
    skipped = 0
    errors = []

    print(f"Fixing {len(FAILING_PAPERS)} failing whitepapers...")
    print("=" * 60)

    for paper_num in FAILING_PAPERS:
        filepath = find_file_for_paper(paper_num)
        if not filepath:
            errors.append(f"PAPER_{paper_num:03d}: file not found")
            continue

        fname = os.path.basename(filepath)
        try:
            changed = fix_paper(filepath)
            if changed:
                fixed += 1
                print(f"  FIXED  PAPER_{paper_num:03d}  {fname}")
            else:
                skipped += 1
                print(f"  SKIP   PAPER_{paper_num:03d}  (no changes needed)")
        except Exception as e:
            errors.append(f"PAPER_{paper_num:03d}: {e}")
            print(f"  ERROR  PAPER_{paper_num:03d}  {e}")

    print("=" * 60)
    print(f"Fixed: {fixed}  |  Skipped: {skipped}  |  Errors: {len(errors)}")
    if errors:
        print("\nErrors:")
        for e in errors:
            print(f"  {e}")


if __name__ == '__main__':
    main()
