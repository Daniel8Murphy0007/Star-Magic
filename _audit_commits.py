import subprocess, re

result = subprocess.run(['git', 'log', '--oneline', '-25', '--name-status'], capture_output=True, text=True)
lines = result.stdout.splitlines()

commits = []
current = None
for line in lines:
    m = re.match(r'^([0-9a-f]{7,40}) (.+)', line)
    if m:
        current = {'hash': m.group(1), 'msg': m.group(2), 'md_mod': [], 'pdf_add': []}
        commits.append(current)
    elif current and line.startswith('M\t') and 'whitepapers/PAPER_' in line:
        current['md_mod'].append(line.split('\t')[1])
    elif current and (line.startswith('A\t') or line.startswith('M\t')) and 'pdf/PAPER_' in line:
        current['pdf_add'].append(line.split('\t')[1])

print('Commits with MD changes but NO PDF updates:')
for c in commits:
    if c['md_mod'] and not c['pdf_add']:
        print(f"  {c['hash']} -- {len(c['md_mod'])} md changed -- {c['msg'][:70]}")

print()
print('Summary of all 25 commits:')
for c in commits:
    print(f"  {c['hash']}  md_mod={len(c['md_mod'])}  pdf_add={len(c['pdf_add'])}  {c['msg'][:55]}")
