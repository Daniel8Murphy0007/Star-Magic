#!/usr/bin/env python3
"""diagnose_and_fix.py - Fix LaTeX math errors in 4 specific failing papers."""
import re, os, subprocess, tempfile

PAPERS = [
    'whitepapers/PAPER_072_Red_Dwarf_Reactor_Physics_UQFF.md',
    'whitepapers/PAPER_114_EP07_ParkerProbe_Heliosheath_Proof.md',
    'whitepapers/PAPER_129_UQFF_Triadic_3C273_Jet_NegativeTime_N13.md',
    'whitepapers/PAPER_137_UQFF_26QuantumLevels_EnergyLadder_E0to10n_Higgs_GalacticVacuum.md',
]

def read_src(path):
    for enc in ('utf-8', 'cp1252', 'latin-1'):
        try:
            return open(path, encoding=enc).read()
        except UnicodeDecodeError:
            continue
    return open(path, encoding='utf-8', errors='replace').read()

def fix_paper(text):
    """Apply targeted surgical fixes to math content in the markdown source."""

    # ── Fix 1: sanitize problematic content inside \text{...} ──────────────────
    # Handles: % (comment char), — en/em dashes (math-active via unicode-math),
    # and ^{n} superscripts (math-only syntax invalid in text mode).
    # The pattern allows ONE level of nested braces so ^{10} is captured fully.
    _TEXT_PAT = re.compile(r'\\text\{((?:[^{}]|\{[^{}]*\})*)\}')

    def sanitize_text_cmd(m):
        inner = m.group(1)
        inner = re.sub(r'(?<!\\)%', r'\\%', inner)              # % → \%
        inner = inner.replace('\u2014', '---')                    # — em dash
        inner = inner.replace('\u2013', '--')                     # – en dash
        inner = re.sub(r'(?<!\\)_', r'\\_', inner)              # _ → \_ (math-active in text)
        inner = re.sub(r'\^\{([^}]*)\}', r'\\textsuperscript{\1}', inner)  # ^{n}
        return r'\text{' + inner + '}'

    def fix_math_block(m):
        return _TEXT_PAT.sub(sanitize_text_cmd, m.group(0))

    text = re.sub(r'\$\$.*?\$\$', fix_math_block, text, flags=re.DOTALL)
    text = re.sub(r'\\\[.*?\\\]',  fix_math_block, text, flags=re.DOTALL)

    # ── Fix 2: replace ÷ with \ensuremath{\div} in math blocks ─────────────────
    # \ensuremath works in both text and math mode; \div alone requires math.
    def replace_div(m):
        return m.group(0).replace('\u00f7', r'\ensuremath{\div}')

    text = re.sub(r'\$\$.*?\$\$', replace_div, text, flags=re.DOTALL)
    text = re.sub(r'\\\[.*?\\\]',  replace_div, text, flags=re.DOTALL)

    # ── Fix 3: remove [ ] brackets wrapping \text{...} in math blocks ──────────
    # Affects PAPER_129: [\text{wrong — use t_n/3 form}] and
    #                    [\text{zero-crossings of cos}(\pi t_n)]
    text = re.sub(r'\[\\text\{([^}]*)\}([^\]]*)\]',
                  lambda m: r'\text{' + m.group(1) + '}' + m.group(2),
                  text)

    # ── Fix 4: remove extra } after \right)^2 (fixed PAPER_239 previously) ────
    text = re.sub(r'(\\right\)\^2)\}(\s*\$)', r'\1\2', text)

    return text


def try_generate_pdf(path, pdf_path):
    """Try to generate PDF using pandoc with fixed temp file."""
    text = read_src(path)
    text = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ufffd\u202f]', '', text)
    text = fix_paper(text)

    tmp = tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', suffix='.md', delete=False)
    tmp.write(text)
    tmp.close()

    BASE_CMD = [
        'pandoc', '--pdf-engine=xelatex',
        '-V', 'geometry:margin=1in', '-V', 'fontsize=11pt',
        '-V', 'documentclass=article',
        '-H', 'pdf_header.tex',
        '--pdf-engine-opt=-interaction=nonstopmode',
        '--from=markdown-yaml_metadata_block-raw_tex+smart',
        '--standalone', '--wrap=none',
    ]
    r = subprocess.run(BASE_CMD + [tmp.name, '-o', pdf_path], capture_output=True, timeout=120)
    os.unlink(tmp.name)
    
    if r.returncode == 0 and os.path.exists(pdf_path):
        print(f'  OK: {os.path.getsize(pdf_path)//1024}KB')
        return True
    
    stderr = r.stderr.decode('utf-8', errors='replace')
    # Get first error line
    for line in stderr.split('\n'):
        if '!' in line or 'Error' in line:
            print(f'  FAIL: {line.strip()[:100]}')
            break
    return False

if __name__ == '__main__':
    os.makedirs('pdf', exist_ok=True)
    for path in PAPERS:
        fname = os.path.basename(path)
        pdf_path = f"pdf/{fname.replace('.md', '.pdf')}"
        print(f'\n{fname}:')
        try_generate_pdf(path, pdf_path)
