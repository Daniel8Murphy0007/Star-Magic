"""Scan all 17 merged PAPER_1109-1125 files for formatting issues."""
import re
from pathlib import Path

issues = {}
for n in range(1109, 1126):
    candidates = sorted(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if not candidates:
        continue
    p = candidates[0]
    text = p.read_text(encoding='utf-8', errors='replace')
    lines = text.split('\n')
    paper_issues = []
    for i, line in enumerate(lines, 1):
        # Encoding artifacts (mojibake from latin-1/windows-1252 misread as utf-8)
        artifacts = ['â€"', 'â€™', 'â€œ', 'â€', 'Ã©', 'Ã¨', 'Ã ', 'â€¦', 'Ã¼', 'Â ', 'â€˜']
        if any(a in line for a in artifacts):
            paper_issues.append(f'  L{i} ENCODING: {line[:120]}')
            continue

        # Unicode Greek/math characters outside of $...$
        greek = 'κβραδΔλμωΩπσΣφΦθΘεηγΓξζ'
        for ch in greek:
            if ch in line:
                # Skip if inside math delimiters on same line
                in_math = bool(re.search(r'\$[^$]*' + re.escape(ch) + r'[^$]*\$', line))
                if not in_math:
                    # Skip table rows and headers - they often have Greek legitimately
                    stripped = line.strip()
                    if not stripped.startswith('|'):
                        paper_issues.append(f'  L{i} UNICODE_GREEK({ch}): {line[:120]}')
                        break

        # \mus, \mum or other undefined units
        if re.search(r'\\mu[^_\s{]', line) and not re.search(r'\$.*\\mu[^_\s{].*\$', line):
            paper_issues.append(f'  L{i} UNDEFINED_UNIT: {line[:120]}')

        # \text{} or \mathrm{} outside math
        if re.search(r'[^$]\\(text|mathrm)\{', line):
            before = line.split('\\text')[0].split('\\mathrm')[0]
            if before.count('$') % 2 == 0:
                paper_issues.append(f'  L{i} CMD_OUTSIDE_MATH: {line[:120]}')

        # \\( ... \\) inline math (pandoc usually handles but can break)
        if r'\(' in line and r'\)' in line:
            paper_issues.append(f'  L{i} PAREN_MATH: {line[:120]}')

    if paper_issues:
        issues[p.name] = paper_issues

total = sum(len(v) for v in issues.values())
print(f'Found {total} issues in {len(issues)} files:\n')
for fname, errs in issues.items():
    print(f'{fname}:')
    for e in errs[:20]:
        print(e)
    if len(errs) > 20:
        print(f'  ... and {len(errs)-20} more')
    print()
