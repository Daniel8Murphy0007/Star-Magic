"""For each affected file, walk full history and print ?-count per commit so we
can identify the EXACT corruption commit (max positive delta)."""
import subprocess
TARGETS = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]
for stem in TARGETS:
    path = f'whitepapers/{stem}.md'
    r = subprocess.run(['git','log','--reverse','--format=%H %ci %s','--', path],
                       capture_output=True, text=True, encoding='utf-8')
    lines = r.stdout.strip().split('\n')
    print(f'\n=== {stem} ===')
    prev_q = 0
    max_delta = 0
    max_info = None
    for line in lines:
        sha = line.split()[0]
        rest = line[len(sha)+1:]
        r2 = subprocess.run(['git','show', f'{sha}:{path}'], capture_output=True)
        txt = r2.stdout.decode('utf-8', errors='replace')
        qc = txt.count('?')
        delta = qc - prev_q
        if abs(delta) >= 3:
            print(f'  {sha[:10]} ?={qc:3} d={delta:+3}  {rest[:80]}')
        if delta > max_delta:
            max_delta = delta
            max_info = (sha, qc, prev_q, rest)
        prev_q = qc
    print(f'  >>> MAX JUMP: {max_info[0][:10]} prev?={max_info[2]} new?={max_info[1]} delta=+{max_delta}')
    print(f'      msg: {max_info[3][:100]}')
