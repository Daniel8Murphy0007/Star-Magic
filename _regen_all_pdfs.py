"""Regenerate PDFs for all 1,119 modified whitepapers using parallel workers."""
import sys
sys.path.insert(0, '.')
from generate_pdfs import generate_pdf
from pathlib import Path
import concurrent.futures, time

papers = sorted(Path('whitepapers').glob('PAPER_*.md'))

ok_count = 0
fail_count = 0
fails = []

start = time.time()

# Use 2 workers (matching generate_pdfs.py WORKERS constant)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
    futures = {ex.submit(generate_pdf, p): p for p in papers}
    for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
        p = futures[fut]
        try:
            ok = fut.result()
        except Exception as e:
            ok = False
        if ok:
            ok_count += 1
        else:
            fail_count += 1
            fails.append(p.name)
        if i % 100 == 0:
            elapsed = time.time() - start
            print(f'  [{i}/{len(papers)}] OK={ok_count} FAIL={fail_count}  elapsed={elapsed:.0f}s')

elapsed = time.time() - start
print(f'\nDone in {elapsed:.0f}s')
print(f'  OK:   {ok_count}')
print(f'  FAIL: {fail_count}')
if fails:
    print('Failed papers:')
    for f in fails[:30]:
        print('  ', f)
