import re, os
stems = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]
out = []
for s in stems:
    t = open(f'whitepapers/{s}.md', encoding='utf-8').read()
    out.append(f'\n========== {s}  ?-count={t.count("?")} ==========')
    for m in re.finditer(r'.{0,55}\?.{0,55}', t):
        ctx = m.group(0).replace('\n', ' ')
        out.append(f'  | {ctx}')
open('_q_contexts.txt','w',encoding='utf-8').write('\n'.join(out))
print('wrote _q_contexts.txt', len(out), 'lines')
