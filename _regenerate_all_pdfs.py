#!/usr/bin/env python3
"""
Regenerate ALL PDFs from upgraded whitepapers using pandoc + xelatex.
Processes in batches of 20 for progress reporting.
"""
import os
import glob
import subprocess
import sys
from datetime import datetime

WHITEPAPERS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'whitepapers')
PDF_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pdf')
PANDOC = r'C:\Users\tmsjd\AppData\Local\Pandoc\pandoc.exe'

os.makedirs(PDF_DIR, exist_ok=True)

files = sorted(glob.glob(os.path.join(WHITEPAPERS_DIR, 'PAPER_*.md')))
total = len(files)
success = 0
failed = []

print(f'Regenerating {total} PDFs...')
print(f'Started: {datetime.now().strftime("%H:%M:%S")}')
print()

for i, md_path in enumerate(files):
    basename = os.path.splitext(os.path.basename(md_path))[0]
    pdf_path = os.path.join(PDF_DIR, basename + '.pdf')

    try:
        result = subprocess.run(
            [
                PANDOC,
                md_path,
                '-o', pdf_path,
                '--pdf-engine=xelatex',
                '-V', 'geometry:margin=1in',
                '-V', 'fontsize=11pt',
                '-V', 'mainfont=Segoe UI',
                '-V', 'monofont=Consolas',
                '--wrap=auto',
                '--columns=80',
            ],
            capture_output=True,
            text=True,
            timeout=120
        )

        if result.returncode == 0:
            success += 1
        else:
            # Try fallback without xelatex font options
            result2 = subprocess.run(
                [
                    PANDOC,
                    md_path,
                    '-o', pdf_path,
                    '--pdf-engine=xelatex',
                    '-V', 'geometry:margin=1in',
                    '-V', 'fontsize=11pt',
                    '--wrap=auto',
                    '--columns=80',
                ],
                capture_output=True,
                text=True,
                timeout=120
            )
            if result2.returncode == 0:
                success += 1
            else:
                # Final fallback: strip problematic LaTeX
                result3 = subprocess.run(
                    [
                        PANDOC,
                        md_path,
                        '-o', pdf_path,
                        '--pdf-engine=xelatex',
                        '-V', 'geometry:margin=1in',
                    ],
                    capture_output=True,
                    text=True,
                    timeout=120
                )
                if result3.returncode == 0:
                    success += 1
                else:
                    failed.append((basename, result3.stderr[-200:] if result3.stderr else 'Unknown error'))
                    print(f'  FAILED: {basename}')

    except subprocess.TimeoutExpired:
        failed.append((basename, 'TIMEOUT'))
        print(f'  TIMEOUT: {basename}')
    except Exception as e:
        failed.append((basename, str(e)))
        print(f'  ERROR: {basename}: {e}')

    if (i + 1) % 50 == 0 or i == total - 1:
        print(f'  [{i+1:4d}/{total}] {success} OK, {len(failed)} failed')

print()
print(f'Finished: {datetime.now().strftime("%H:%M:%S")}')
print(f'Success: {success}/{total}')
print(f'Failed: {len(failed)}')

if failed:
    print('\nFailed papers:')
    for name, err in failed:
        print(f'  {name}: {err[:100]}')

# Save results
import json
with open('pdf_regen_stats.json', 'w') as f:
    json.dump({
        'total': total,
        'success': success,
        'failed': len(failed),
        'failed_list': [name for name, _ in failed]
    }, f, indent=2)
