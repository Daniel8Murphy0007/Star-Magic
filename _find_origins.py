stems = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]
import subprocess
for s in stems:
    p = f'whitepapers/{s}.md'
    # find first commit that added the file
    r = subprocess.run(['git','log','--oneline','--diff-filter=A','--all','--',p], capture_output=True, text=True)
    print(f'{s}: {r.stdout.strip()}')
