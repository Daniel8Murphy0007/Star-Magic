"""Batch regenerate ALL PDFs from upgraded whitepapers.
Uses pandoc + xelatex with Segoe UI font and 1-inch margins."""
import os
import subprocess
import time
import json
import sys

WHITEPAPERS_DIR = 'whitepapers'
PDF_DIR = 'pdf'
PANDOC = r'C:\Users\tmsjd\AppData\Local\Pandoc\pandoc.exe'

# Base pandoc arguments
BASE_ARGS = [
    PANDOC,
    '--pdf-engine=xelatex',
    '-V', 'geometry:margin=1in',
    '-V', 'mainfont:Segoe UI',
    '-V', 'monofont:Consolas',
]

os.makedirs(PDF_DIR, exist_ok=True)

# Get all whitepaper files
files = sorted([f for f in os.listdir(WHITEPAPERS_DIR) if f.endswith('.md')])
total = len(files)

successes = 0
skipped = 0
failures = []
start_time = time.time()

# Cutoff: skip PDFs regenerated after this script was first run (within last 2 hours)  
cutoff = time.time() - 7200  # 2 hours ago

for i, fn in enumerate(files, 1):
    input_path = os.path.join(WHITEPAPERS_DIR, fn)
    pdf_name = fn.replace('.md', '.pdf')
    output_path = os.path.join(PDF_DIR, pdf_name)
    
    # Skip if PDF already exists and is recent (regenerated this session)
    if os.path.exists(output_path) and os.path.getmtime(output_path) > cutoff and os.path.getsize(output_path) > 1000:
        successes += 1
        skipped += 1
        continue
    
    try:
        result = subprocess.run(
            BASE_ARGS + [input_path, '-o', output_path],
            capture_output=True, text=True, timeout=120,
            encoding='utf-8', errors='replace'
        )
        
        if result.returncode == 0 and os.path.exists(output_path):
            successes += 1
        elif os.path.exists(output_path) and os.path.getsize(output_path) > 1000:
            # PDF created despite non-zero exit (font warnings cause this)
            successes += 1
        else:
            # Actual failure
            err_snippet = result.stderr[-200:] if result.stderr else 'no stderr'
            failures.append({'file': fn, 'error': err_snippet})
    except subprocess.TimeoutExpired:
        failures.append({'file': fn, 'error': 'TIMEOUT (120s)'})
    except Exception as e:
        failures.append({'file': fn, 'error': str(e)})
    
    # Progress every 50 papers
    if i % 50 == 0 or i == total:
        elapsed = time.time() - start_time
        rate = i / elapsed if elapsed > 0 else 0
        eta = (total - i) / rate if rate > 0 else 0
        print(f"[{i}/{total}] {successes} OK, {len(failures)} failed | "
              f"{elapsed:.0f}s elapsed, ~{eta:.0f}s remaining", flush=True)
    elif i % 10 == 0:
        print(f"  [{i}]...", end='', flush=True)

elapsed = time.time() - start_time
print(f"\n=== COMPLETE ===")
print(f"Total: {total}")
print(f"Skipped (already done): {skipped}")
print(f"Success: {successes}")
print(f"Failed: {len(failures)}")
print(f"Time: {elapsed:.1f}s ({elapsed/60:.1f}m)")

if failures:
    print(f"\n--- FAILURES ({len(failures)}) ---")
    for f in failures[:20]:
        print(f"  {f['file']}: {f['error'][:100]}")
    if len(failures) > 20:
        print(f"  ... and {len(failures) - 20} more")
    
    # Save failure list
    with open('pdf_failures.json', 'w') as fh:
        json.dump(failures, fh, indent=2)
    print(f"\nFull failure list saved to pdf_failures.json")
