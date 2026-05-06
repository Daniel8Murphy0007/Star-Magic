"""
Relocate the ## Supplementary Derivations block in PAPER_1109-1125
from its current position (before the short ## References section)
to the end of each document (after ## B. VDS/DVP/BSH, before
### Key References if it exists, otherwise at EOF).

This fixes the formatting error where ## References appears mid-document
between the Supplementary Derivations and ## Session 225 content.
"""
import re
from pathlib import Path
import sys
sys.path.insert(0, '.')
from generate_pdfs import generate_pdf

SUPP_START_PAT = re.compile(r'\n---\n\n## Supplementary Derivations \(Polylogarithmic / VDS\)')
SUPP_END_PAT   = re.compile(r'\n---\n## References\b')         # trailing --- + short refs (no blank line)
KEY_REFS_PAT   = re.compile(r'\n\n### Key References with arXiv/DOI Identifiers\b')

results = []

for n in range(1109, 1126):
    candidates = sorted(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if not candidates:
        print(f'PAPER_{n}: file not found')
        continue
    p = candidates[0]
    text = p.read_text(encoding='utf-8')

    # --- Locate the supplement block ---
    m_start = SUPP_START_PAT.search(text)
    if not m_start:
        print(f'{p.name}: WARNING — no supplement block found, skipping')
        continue

    supp_block_start = m_start.start()   # position of the leading \n before ---

    # End of supplement block = the \n before the trailing ---\n## References
    m_end = SUPP_END_PAT.search(text, supp_block_start + 1)
    if not m_end:
        print(f'{p.name}: WARNING — no end marker (---\\n## References) found, skipping')
        continue

    supp_block_end = m_end.start()       # position of \n before trailing ---

    supplement_text = text[supp_block_start:supp_block_end]
    # supplement_text looks like: \n---\n\n## Supplementary Derivations\n...\n---\n

    # --- Remove supplement from current position ---
    text_no_supp = text[:supp_block_start] + text[supp_block_end:]

    # --- Find insertion point ---
    m_key = KEY_REFS_PAT.search(text_no_supp)
    if m_key:
        # Insert before ### Key References (at the double-newline before it)
        insert_pos = m_key.start()       # position of first \n of \n\n### Key References
        new_text = text_no_supp[:insert_pos] + supplement_text + '\n' + text_no_supp[insert_pos:]
    else:
        # Append at end of file
        # Make sure there's a clean separator
        stripped_end = text_no_supp.rstrip('\n')
        new_text = stripped_end + supplement_text + '\n'

    p.write_text(new_text, encoding='utf-8')
    lines = new_text.count('\n')
    print(f'{p.name}: relocated supplement ({supp_block_end - supp_block_start}B) -> {"before Key References" if m_key else "EOF"} | {lines} lines')

    # --- Regenerate PDF ---
    result = generate_pdf(p)
    ok = 'OK' if result[2] else 'FAIL'
    print(f'  PDF: {ok} {result[3] or result[4] or ""}')
    results.append((p.name, ok))

print(f'\n=== SUMMARY: {sum(1 for _,s in results if s=="OK")}/{len(results)} OK ===')
for name, status in results:
    print(f'  {status}  {name}')
