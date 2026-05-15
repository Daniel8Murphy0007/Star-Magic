"""Walk full git history of PAPER_020 to find when ? characters appeared."""
import subprocess
r = subprocess.run(['git','log','--reverse','--format=%H',
                    '--', 'whitepapers/PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md'],
                   capture_output=True, text=True)
shas = r.stdout.strip().split('\n')
print(f'Total commits: {len(shas)}')
prev_q = None
for sha in shas:
    r2 = subprocess.run(['git','show', f'{sha}:whitepapers/PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md'],
                        capture_output=True)
    txt = r2.stdout.decode('utf-8', errors='replace')
    qc = txt.count('?')
    flag = ''
    if prev_q is not None and qc != prev_q:
        flag = f'  <-- delta {qc-prev_q:+d}'
    # Also get short message
    r3 = subprocess.run(['git','log','-1','--format=%s', sha], capture_output=True, text=True)
    msg = r3.stdout.strip()[:80]
    print(f'  {sha[:10]}  size={len(txt):6}  ?={qc:3}{flag}  | {msg}')
    prev_q = qc
