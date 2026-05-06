import subprocess, re

result = subprocess.run(['git', 'log', '--oneline', '-25', '--name-status'], capture_output=True, text=True)
lines = result.stdout.splitlines()

modified_md_nums = set()
pdf_touched_nums = set()

for line in lines:
    if line.startswith(('A\t','M\t')):
        path = line.split('\t')[1]
        m = re.search(r'PAPER_0*(\d+)', path)
        if not m:
            continue
        n = int(m.group(1))
        if 'whitepapers/' in path and path.endswith('.md'):
            modified_md_nums.add(n)
        elif 'pdf/' in path and path.endswith('.pdf'):
            pdf_touched_nums.add(n)

stale = sorted(modified_md_nums - pdf_touched_nums)
print(f'MD files changed in last 25 commits:    {len(modified_md_nums)}')
print(f'PDF files touched in last 25 commits:   {len(pdf_touched_nums)}')
print(f'Whitepapers with stale PDF (no regen):  {len(stale)}')
print()

# Group into ranges for readability
if stale:
    ranges = []
    start = stale[0]
    end = stale[0]
    for n in stale[1:]:
        if n == end + 1:
            end = n
        else:
            ranges.append((start, end))
            start = n
            end = n
    ranges.append((start, end))
    for s, e in ranges:
        if s == e:
            print(f'  PAPER_{s:04d}')
        else:
            print(f'  PAPER_{s:04d} – PAPER_{e:04d}  ({e-s+1} papers)')
