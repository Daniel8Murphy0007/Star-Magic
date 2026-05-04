# -*- coding: utf-8 -*-
"""Regenerate PDFs for the 716 whitepapers that had kg/m3 -> J/m3 unit fix."""
from pathlib import Path
import subprocess, sys, os

PYTHON = sys.executable
PANDOC = Path(os.environ.get('LOCALAPPDATA', '')) / 'Pandoc' / 'pandoc.exe'
XELATEX = Path(r'C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe')

def generate_pdf(md_path: Path) -> bool:
    pdf_path = Path('pdf') / (md_path.stem + '.pdf')
    cmd = [
        str(PANDOC), str(md_path),
        '--pdf-engine', str(XELATEX),
        '-V', 'geometry:margin=1in',
        '-V', 'fontsize=11pt',
        '-V', 'colorlinks=true',
        '--highlight-style=tango',
        '-o', str(pdf_path)
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, timeout=180)
        return r.returncode == 0
    except Exception:
        return False

if __name__ == '__main__':
    # Find whitepapers that now have J/m3 (just regenerate all that were fixed)
    wp_dir = Path('whitepapers')
    pdf_dir = Path('pdf')
    pdf_dir.mkdir(exist_ok=True)

    # The 716 fixed papers - find those with J/m^{3} in content (freshly fixed)
    import re
    targets = []
    for md in sorted(wp_dir.glob('*.md')):
        txt = md.read_text('utf-8', errors='replace')
        if r'\text{J/m}^3' in txt:
            targets.append(md)

    print(f'Regenerating {len(targets)} PDFs...')
    ok = fail = 0
    for i, md in enumerate(targets, 1):
        success = generate_pdf(md)
        if success:
            ok += 1
        else:
            fail += 1
        if i % 50 == 0:
            print(f'  [{i}/{len(targets)}] ok={ok} fail={fail}')

    print(f'\nDone: {ok} OK, {fail} failed out of {len(targets)}')
