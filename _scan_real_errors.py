"""Find real formatting errors: encoding artifacts and genuine LaTeX issues
in PAPER_1109-1125 merged files."""
import re
from pathlib import Path

ARTIFACTS = {
    'â€"': '—',      # em dash
    'â€™': "'",      # right single quote
    'â€œ': '"',      # left double quote
    'â€': '"',       # right double quote
    'â€¦': '…',      # ellipsis
    'Ã©': 'é',
    'Ã¨': 'è',
    'Ã ': 'à',
    'Ã¼': 'ü',
    'Â ': ' ',       # non-breaking space artifact
    'â€˜': "'",      # left single quote
    'Ã¢': 'â',
    'Ã¯': 'ï',
}

# Genuine LaTeX errors (not false positives from scanner)
LATEX_ERRORS = [
    # \text{} or \mathrm{} NOT inside any math delimiters on the same line
    # i.e. in pure prose paragraphs
    (r'^(?![^$]*\$)[^$]*\\text\{', 'TEXT_IN_PROSE'),
    # Undefined units: \mus, \mum, \mbox used wrong
    (r'\\mus\b', 'UNDEF_MUS'),
    (r'\\mum\b', 'UNDEF_MUM'),
    # Mixed inline math delimiters \(...\) which pandoc sometimes mangles
    (r'\\\(.*\\\)', 'PAREN_INLINE_MATH'),
]

total_fixes = 0
issues_by_file = {}

for n in range(1109, 1126):
    candidates = sorted(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if not candidates:
        continue
    p = candidates[0]
    text = p.read_text(encoding='utf-8', errors='replace')
    lines = text.split('\n')

    file_issues = []

    # --- Encoding artifacts ---
    for i, line in enumerate(lines, 1):
        found = [(artifact, replacement)
                 for artifact, replacement in ARTIFACTS.items()
                 if artifact in line]
        for artifact, replacement in found:
            count = line.count(artifact)
            file_issues.append((i, 'ENCODING', artifact, replacement, line.strip()[:100]))

    # --- Genuine LaTeX errors ---
    for i, line in enumerate(lines, 1):
        # Skip lines that are inside math blocks (start with $$)
        stripped = line.strip()
        if stripped.startswith('$$') and stripped.endswith('$$') and len(stripped) > 4:
            continue  # display math line, \text{} is valid

        for pattern, label in LATEX_ERRORS:
            if re.search(pattern, line):
                file_issues.append((i, label, '', '', line.strip()[:100]))
                break

    if file_issues:
        issues_by_file[p.name] = file_issues

# Print summary
total = sum(len(v) for v in issues_by_file.values())
encoding_total = sum(1 for errs in issues_by_file.values()
                     for _, kind, *_ in errs if kind == 'ENCODING')
latex_total = sum(1 for errs in issues_by_file.values()
                  for _, kind, *_ in errs if kind != 'ENCODING')
print(f'Total issues: {total}  (encoding: {encoding_total}, latex: {latex_total})\n')

for fname, errs in issues_by_file.items():
    enc = [e for e in errs if e[1] == 'ENCODING']
    lat = [e for e in errs if e[1] != 'ENCODING']
    print(f'{fname}: {len(enc)} encoding, {len(lat)} latex')
    for e in enc[:5]:
        print(f'  L{e[0]} [{e[2]}] -> [{e[3]}]  | {e[4]}')
    for e in lat[:5]:
        print(f'  L{e[0]} {e[1]}: {e[4]}')
    print()
