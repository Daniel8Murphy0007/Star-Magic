"""
upgrade_latex.py — Bulk arXiv LaTeX upgrade for Star-Magic whitepapers.

Converts Unicode math symbols to proper LaTeX commands in all PAPER_*.md files.
Handles three contexts:
  1. YAML frontmatter        → unchanged
  2. Code blocks / inline code → unchanged
  3. Existing $...$ / $$...$$ math → replace Unicode with LaTeX commands (no extra wrapping)
  4. Plain text               → wrap Unicode math in $...$

Run:  python upgrade_latex.py [glob_pattern]
      python upgrade_latex.py "whitepapers/PAPER_*.md"   (default)
"""

import re
import os
import sys
import glob

# ---------------------------------------------------------------------------
# Replacement tables
# ---------------------------------------------------------------------------

# Math chars: Unicode → LaTeX command (without $)
MATH_CHARS = {
    # Greek lowercase
    'α': r'\alpha',    'β': r'\beta',     'γ': r'\gamma',    'δ': r'\delta',
    'ε': r'\varepsilon','ζ': r'\zeta',    'η': r'\eta',      'θ': r'\theta',
    'ι': r'\iota',     'κ': r'\kappa',    'λ': r'\lambda',   'μ': r'\mu',
    'ν': r'\nu',       'ξ': r'\xi',       'π': r'\pi',       'ρ': r'\rho',
    'σ': r'\sigma',    'τ': r'\tau',      'υ': r'\upsilon',  'φ': r'\phi',
    'χ': r'\chi',      'ψ': r'\psi',      'ω': r'\omega',
    # Greek uppercase (those with dedicated LaTeX commands)
    'Γ': r'\Gamma',    'Δ': r'\Delta',    '∆': r'\Delta',    'Θ': r'\Theta',
    'Λ': r'\Lambda',   'Ξ': r'\Xi',       'Π': r'\Pi',       'Σ': r'\Sigma',
    'Υ': r'\Upsilon',  'Φ': r'\Phi',      'Ψ': r'\Psi',      'Ω': r'\Omega',
    # Math operators
    '×': r'\times',    '÷': r'\div',      '·': r'\cdot',
    '→': r'\to',       '←': r'\leftarrow','↔': r'\leftrightarrow',
    '↑': r'\uparrow',  '↓': r'\downarrow',
    '∑': r'\sum',      '∏': r'\prod',     '∫': r'\int',      '∮': r'\oint',
    '∂': r'\partial',  '∇': r'\nabla',
    '∞': r'\infty',    '±': r'\pm',       '∓': r'\mp',
    '≈': r'\approx',   '≠': r'\neq',      '≤': r'\leq',      '≥': r'\geq',
    '≡': r'\equiv',    '≃': r'\simeq',    '≅': r'\cong',     '∝': r'\propto',
    '∈': r'\in',       '∉': r'\notin',    '⊂': r'\subset',   '⊃': r'\supset',
    '⊆': r'\subseteq', '⊇': r'\supseteq',
    '∪': r'\cup',      '∩': r'\cap',
    '∀': r'\forall',   '∃': r'\exists',   '∄': r'\nexists',
    '∼': r'\sim',      '⊕': r'\oplus',    '⊗': r'\otimes',
    '∘': r'\circ',
    # Blackboard bold
    'ℚ': r'\mathbb{Q}','ℝ': r'\mathbb{R}','ℂ': r'\mathbb{C}',
    'ℤ': r'\mathbb{Z}','ℕ': r'\mathbb{N}',
    # Special
    '∓': r'\mp',
}

MATH_SET = set(MATH_CHARS.keys())

# Replaced unconditionally (no $ wrapping needed)
UNCONDITIONAL = {
    '−': '-',   # U+2212 minus sign → ASCII hyphen
}

# Subscript/superscript digits → plain digits (strip glyph)
SUBSCRIPT_DIGITS = {
    '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4',
    '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    '⁺': '+', '⁻': '-', '⁽': '(', '⁾': ')',
}

SUBSCRIPT_LETTERS = {
    'ᵢ': 'i', 'ₐ': 'a', 'ₑ': 'e', 'ₒ': 'o',
}

# Superscript digits → ^N (stripped to plain in text, useful in math context)
SUPERSCRIPT_DIGITS = {
    '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4',
    '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
}

ALL_REPLACEABLE = (set(MATH_CHARS) | set(UNCONDITIONAL) |
                   set(SUBSCRIPT_DIGITS) | set(SUBSCRIPT_LETTERS) |
                   set(SUPERSCRIPT_DIGITS))

# ---------------------------------------------------------------------------
# Processing helpers
# ---------------------------------------------------------------------------

def replace_in_math(text):
    """Replace Unicode inside an already-existing math block."""
    for u, l in MATH_CHARS.items():
        text = text.replace(u, l)
    for u, l in UNCONDITIONAL.items():
        text = text.replace(u, l)
    for u, l in SUBSCRIPT_DIGITS.items():
        text = text.replace(u, l)
    for u, l in SUBSCRIPT_LETTERS.items():
        text = text.replace(u, l)
    for u, l in SUPERSCRIPT_DIGITS.items():
        text = text.replace(u, l)
    return text


def replace_in_text(text):
    """Replace Unicode in plain text (outside math/code blocks)."""
    result = []
    i = 0
    while i < len(text):
        c = text[i]
        if c in UNCONDITIONAL:
            result.append(UNCONDITIONAL[c])
        elif c in SUBSCRIPT_DIGITS:
            result.append(SUBSCRIPT_DIGITS[c])
        elif c in SUBSCRIPT_LETTERS:
            # subscript letter like ᵢ: try to join with previous alphanum as _i
            result.append(SUBSCRIPT_LETTERS[c])
        elif c in SUPERSCRIPT_DIGITS:
            # superscript digits: strip to plain
            result.append(SUPERSCRIPT_DIGITS[c])
        elif c in MATH_SET:
            result.append('$' + MATH_CHARS[c] + '$')
        elif c == '√':
            # Handle √N or √(expr) — consume following char(s) into \sqrt{...}
            j = i + 1
            if j < len(text) and text[j] == '(':
                # √(...) → $\sqrt{...}$
                k = text.find(')', j)
                if k != -1:
                    result.append('$\\sqrt{' + text[j+1:k] + '}$')
                    i = k + 1
                    continue
            elif j < len(text) and (text[j].isdigit() or text[j].isalpha() or text[j] == '{'):
                result.append('$\\sqrt{' + text[j] + '}$')
                i = j + 1
                continue
            # fallback
            result.append('$\\sqrt{}$')
        else:
            result.append(c)
        i += 1
    return ''.join(result)


# ---------------------------------------------------------------------------
# Main file processor
# ---------------------------------------------------------------------------

def process_content(content):
    """Process a full markdown file content."""
    if not any(c in content for c in ALL_REPLACEABLE):
        return content   # fast path: nothing to do

    result = []
    i = 0
    n = len(content)

    # Skip YAML frontmatter (--- ... ---)
    if content.startswith('---\n') or content.startswith('---\r\n'):
        end = content.find('\n---', 4)
        if end != -1:
            # include up to and including the closing ---\n
            close_end = end + 1  # skip \n
            while close_end < n and content[close_end] == '-':
                close_end += 1
            # skip trailing newline
            if close_end < n and content[close_end] in ('\n', '\r'):
                close_end += 1
            result.append(content[:close_end])
            i = close_end

    while i < n:
        # -------------------------------------------------------------------
        # Fenced code block  ```...```  or  ~~~...~~~
        # -------------------------------------------------------------------
        if content[i:i+3] in ('```', '~~~'):
            fence = content[i:i+3]
            # find closing fence on its own line
            close_pat = '\n' + fence
            end = content.find(close_pat, i + 3)
            if end == -1:
                result.append(content[i:])
                break
            # include closing fence + any trailing chars until newline
            end2 = end + len(close_pat)
            while end2 < n and content[end2] not in ('\n', '\r'):
                end2 += 1
            if end2 < n:
                end2 += 1  # include the newline
            result.append(content[i:end2])
            i = end2
            continue

        # -------------------------------------------------------------------
        # Inline code  `...`
        # -------------------------------------------------------------------
        if content[i] == '`':
            j = i + 1
            while j < n and content[j] == '`':
                j += 1
            tick_count = j - i
            close = content.find('`' * tick_count, j)
            if close == -1:
                result.append(content[i])
                i += 1
                continue
            result.append(content[i:close + tick_count])
            i = close + tick_count
            continue

        # -------------------------------------------------------------------
        # Display math  $$...$$
        # -------------------------------------------------------------------
        if content[i:i+2] == '$$':
            end = content.find('$$', i + 2)
            if end == -1:
                result.append(content[i])
                i += 1
                continue
            math_content = replace_in_math(content[i+2:end])
            result.append('$$' + math_content + '$$')
            i = end + 2
            continue

        # -------------------------------------------------------------------
        # Inline math  $...$
        # (skip if preceded by $  i.e. part of $$)
        # -------------------------------------------------------------------
        if content[i] == '$' and content[i:i+2] != '$$' and (i == 0 or content[i-1] != '$'):
            # scan for closing $, not crossing paragraph breaks
            j = i + 1
            found_close = -1
            while j < n:
                if content[j] == '$' and content[j:j+2] != '$$':
                    found_close = j
                    break
                if content[j] == '\n' and j + 1 < n and content[j+1] == '\n':
                    break  # paragraph break: not a math span
                j += 1
            if found_close > i:
                math_content = replace_in_math(content[i+1:found_close])
                result.append('$' + math_content + '$')
                i = found_close + 1
                continue

        # -------------------------------------------------------------------
        # Plain text character
        # -------------------------------------------------------------------
        c = content[i]
        if c in ALL_REPLACEABLE:
            # process a run of plain text up to the next special boundary
            j = i
            while j < n:
                nc = content[j]
                if nc in ('`', '$') or content[j:j+3] in ('```', '~~~'):
                    break
                j += 1
            plain = content[i:j]
            result.append(replace_in_text(plain))
            i = j
            continue

        result.append(c)
        i += 1

    return ''.join(result)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    pattern = sys.argv[1] if len(sys.argv) > 1 else 'whitepapers/PAPER_*.md'
    files = sorted(glob.glob(pattern))
    if not files:
        print(f'No files matched: {pattern}')
        sys.exit(1)

    changed_files = []
    for filepath in files:
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                original = f.read()
        except Exception as e:
            print(f'ERR read {filepath}: {e}')
            continue

        new_content = process_content(original)
        if new_content != original:
            try:
                with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
                    f.write(new_content)
                changed_files.append(filepath)
                print(f'UPDATED: {filepath}')
            except Exception as e:
                print(f'ERR write {filepath}: {e}')

    print(f'\nDone. {len(changed_files)}/{len(files)} files updated.')
    # Write list of changed files for PDF regeneration
    if changed_files:
        with open('_changed_files.txt', 'w', encoding='utf-8') as f:
            f.write('\n'.join(changed_files))
        print(f'Changed list written to _changed_files.txt')
