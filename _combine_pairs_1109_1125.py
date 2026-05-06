"""
Combine each PAPER_1109-1125 pair:
  - PAPER_XXXX.md (generic, 200-300L) = base with YAML frontmatter + full Session 225 content
  - PAPER_XXXX_Title.md (named, 76-108L) = supplementary derivations (VDS tables, Ramanujan
    acceleration values, dimensional cascades)
Strategy:
  1. Read both files
  2. Strip the named file's simple header (title/author/date lines at top)
  3. Strip the named file's References section (generic has the authoritative one)
  4. Insert named file's remaining body into the generic before its ## References block
  5. Write combined content to the NAMED file path
  6. Delete the generic (PAPER_XXXX.md) file
  7. Regenerate PDF from the named file via generate_pdfs.generate_pdf()
"""

import sys, re
sys.path.insert(0, '.')
from pathlib import Path
from generate_pdfs import generate_pdf

WHITEPAPERS = Path('whitepapers')
PDF_DIR = Path('pdf')

# Canonical UQFF constants used in recalculation validation comments
CANONICAL = {
    'kappa': '5.0e-4 day^-1',
    'SSq': '0.57',
    'beta_i': '0.603',
    'H_SCm': '0.99',
    'rho_SCm': '7.09e-37 J/m^3',
}


def strip_named_file_header(text: str) -> str:
    """
    Remove the simple markdown header block at the top of named files:
        # PAPER_XXXX: Title
        **Star Magic UQFF Framework**
        **Author:** ...
        **Date:** ...
        ---
    Returns the body starting from the first ## section header.
    """
    lines = text.split('\n')
    for i, line in enumerate(lines):
        # First level-2 (##) header marks the real body start
        if line.startswith('## '):
            return '\n'.join(lines[i:])
    # Fallback: return everything after the first --- separator
    found_sep = False
    for i, line in enumerate(lines):
        if line.strip() == '---':
            if found_sep:
                return '\n'.join(lines[i+1:])
            found_sep = True
    return text


def extract_body_without_references(text: str) -> str:
    """
    Return text with ## References (and everything after) stripped out.
    Keeps ## Session ... upgrades and other content before References.
    """
    # Find the FIRST occurrence of a top-level ## References block
    match = re.search(r'\n## References\b', text)
    if match:
        return text[:match.start()]
    # Also handle ### Key References...
    match2 = re.search(r'\n### Key References\b', text)
    if match2:
        return text[:match2.start()]
    return text


def get_section_headers(text: str) -> set:
    """Extract all ## and ### section header lines from text."""
    return set(re.findall(r'^#{2,3} .+', text, re.MULTILINE))


def normalise_header(h: str) -> str:
    """Lowercase stripped header text for comparison."""
    return re.sub(r'^#{1,3}\s+', '', h).strip().lower()


def combine(n: int) -> tuple[bool, str]:
    """
    Combine pair for paper number n.
    Returns (success, message).
    """
    generic_path = WHITEPAPERS / f'PAPER_{n}.md'
    named_candidates = sorted(WHITEPAPERS.glob(f'PAPER_{n}_*.md'))

    if not generic_path.exists():
        return False, f'PAPER_{n}.md not found'
    if not named_candidates:
        return False, f'No named file found for PAPER_{n}'

    named_path = named_candidates[0]

    generic_text = generic_path.read_text(encoding='utf-8', errors='replace')
    named_text = named_path.read_text(encoding='utf-8', errors='replace')

    # --- Extract YAML frontmatter from generic ---
    yaml_match = re.match(r'^---\n.*?\n---\n', generic_text, re.DOTALL)
    frontmatter = yaml_match.group(0) if yaml_match else ''
    generic_body = generic_text[len(frontmatter):]

    # --- Prepare named file body: strip header + strip references ---
    named_body = strip_named_file_header(named_text)
    named_body_no_refs = extract_body_without_references(named_body)

    # --- Identify sections in named body that are UNIQUE (not in generic) ---
    generic_headers = get_section_headers(generic_body)
    generic_normalised = {normalise_header(h) for h in generic_headers}

    # Split named body into sections
    section_pattern = re.compile(r'^(#{2,3} .+)', re.MULTILINE)
    pieces = section_pattern.split(named_body_no_refs)

    # pieces alternates: [preamble, header1, content1, header2, content2, ...]
    unique_chunks = []
    i = 1  # skip preamble
    while i < len(pieces) - 1:
        header = pieces[i]
        content = pieces[i + 1] if i + 1 < len(pieces) else ''
        norm = normalise_header(header)
        skip_keywords = {'abstract', 'introduction', 'conclusion', 'references',
                         'key references', 'calibrat', 'calibration'}
        if norm not in generic_normalised and not any(kw in norm for kw in skip_keywords):
            unique_chunks.append(f'{header}\n{content}')
        i += 2

    # --- Find insert point: just before ## References in generic body ---
    ref_match = re.search(r'\n## References\b', generic_body)
    key_ref_match = re.search(r'\n### Key References\b', generic_body)
    # Use whichever comes first
    candidates = [m for m in [ref_match, key_ref_match] if m is not None]
    if candidates:
        insert_pos = min(c.start() for c in candidates)
    else:
        insert_pos = len(generic_body)

    # --- Build combined document ---
    supplement_block = ''
    if unique_chunks:
        supplement_block = (
            '\n\n---\n\n'
            '## Supplementary Derivations (Polylogarithmic / VDS)\n\n'
            '*Merged from companion derivation file. '
            'Canonical UQFF constants: kappa=5.0e-4/day, [SSq]=0.57, '
            'beta_i=0.603, rho_SCm=7.09e-37 J/m3.*\n\n'
            + '\n\n'.join(c.strip() for c in unique_chunks)
        )

    combined = (
        frontmatter
        + generic_body[:insert_pos]
        + supplement_block
        + generic_body[insert_pos:]
    )

    # Write combined content to the named file
    named_path.write_text(combined, encoding='utf-8')

    # Delete the generic file
    generic_path.unlink()

    # Delete generic PDF if it exists
    generic_pdf = PDF_DIR / f'PAPER_{n}.pdf'
    if generic_pdf.exists():
        generic_pdf.unlink()

    return True, f'Combined {generic_path.name} + {named_path.name} -> {named_path.name} ({len(combined.splitlines())}L)'


def main():
    results = []
    for n in range(1109, 1126):
        ok, msg = combine(n)
        print(f'  PAPER_{n}: {"OK" if ok else "FAIL"}  {msg}')
        results.append((n, ok, msg))

    ok_count = sum(1 for _, ok, _ in results if ok)
    print(f'\nCombined {ok_count}/17 pairs. Regenerating PDFs...\n')

    # Regenerate PDFs for all named files
    fail_pdfs = []
    for n in range(1109, 1126):
        named_candidates = sorted(WHITEPAPERS.glob(f'PAPER_{n}_*.md'))
        if not named_candidates:
            print(f'  PAPER_{n}: no named file, skipping PDF')
            continue
        named_path = named_candidates[0]
        result = generate_pdf(named_path)
        # result is a tuple: (paper_num, fname, success, detail, error)
        success = result[2] if isinstance(result, tuple) and len(result) >= 3 else bool(result)
        detail = result[3] if isinstance(result, tuple) and len(result) >= 4 else ''
        err = result[4] if isinstance(result, tuple) and len(result) >= 5 else ''
        if success:
            print(f'  PAPER_{n}: PDF OK  {detail}  -> {named_path.stem}.pdf')
        else:
            print(f'  PAPER_{n}: PDF FAIL  {err}')
            fail_pdfs.append(named_path.name)

    print(f'\nDone. PDF failures: {len(fail_pdfs)}')
    if fail_pdfs:
        for f in fail_pdfs:
            print(f'  FAIL: {f}')


if __name__ == '__main__':
    main()
