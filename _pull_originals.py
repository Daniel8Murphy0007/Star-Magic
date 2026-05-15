import subprocess, os, re
files = {
 'PAPER_009b_Aether_String_TRZ_Damping_GW': 'ca2c552d',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF': 'ca2c552d',
 'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations': 'f705451e',
 'PAPER_016b_White_Dwarf_Foreground_UQFF': 'ca2c552d',
 'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA': '2201be5b',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime': '98b86da4',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density': 'eae91218',
}
os.makedirs('_orig_md', exist_ok=True)
for s, sha in files.items():
    r = subprocess.run(['git','show', f'{sha}:whitepapers/{s}.md'], capture_output=True)
    open(f'_orig_md/{s}.md','wb').write(r.stdout)
    cur = open(f'whitepapers/{s}.md','rb').read()
    # decode for ? counts
    orig_q = r.stdout.decode('utf-8', errors='replace').count('?')
    cur_q = cur.decode('utf-8', errors='replace').count('?')
    print(f'{s:65}  orig_size={len(r.stdout):6} cur_size={len(cur):6}  orig_?={orig_q:3} cur_?={cur_q:3}')
