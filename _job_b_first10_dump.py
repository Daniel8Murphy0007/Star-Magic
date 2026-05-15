import os, sys
sys.stdout.reconfigure(encoding='utf-8')

papers = [
    'PAPER_001_GW170817_UQFF_Damping_Analysis.md',
    'PAPER_002_GW190425_Mass_Gap_Interpretation.md',
    'PAPER_003_GW150914_UQFF_vs_LIGO_Strain.md',
    'PAPER_004_GW170817_BNS_Chirp_Phase_Evolution.md',
    'PAPER_005_BH_Merger_Energy_Retention_UQFF.md',
    'PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md',
    'PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF.md',
    'PAPER_008_UQFF_Waveform_Phase_Evolution_Template_Mismatch.md',
    'PAPER_008b_Full_Inspiral_Waveform_UQFF.md',
    'PAPER_009_Damping_Mechanism_Decomposition_UQFF.md',
]

for fn in papers:
    path = os.path.join('whitepapers', fn)
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()
    print('='*78)
    print(fn)
    # First 1200 chars (abstract/intro)
    print('--- HEAD ---')
    print(text[:1200])
    # Last 1500 chars (conclusion, refs, session-225 block)
    print('--- TAIL ---')
    print(text[-1500:])
