import os, re, sys
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

v578_signals = {
    '27-decade ledger':       [r'27[- ]decade', r'PAPER_1170', r'vacuum[- ]energy ledger'],
    'xi=13/3':                [r'\bxi\s*=\s*13/3\b', r'13/3\s*(lock|saturation|closure)', r'PAPER_1171', r'PAPER_1172'],
    'G1-G8 closed Lagrangian':[r'G1[-\s]?G8', r'PAPER_116[0-7]', r'closed Lagrangian'],
    'three-anchor SI':        [r'three-?anchor', r'PAPER_59[0-3]'],
    'P-suite P6-P14':         [r'\bP[6-9]\b', r'\bP1[0-4]\b', r'PAPER_117[4-9]', r'PAPER_1180', r'DESI Y5', r'CMB-?S4'],
    'KK regulator':           [r'KK regulator', r'sub-?mm Yukawa', r'PAPER_1173'],
    'Lambda closure':         [r'PAPER_1156', r'Lambda closure', r'Lambda\s*=\s*0\.002%'],
    'Session 225 YM/BCS':     [r'Session 225', r'BCS phonon', r'Yang-Mills.{0,40}phonon'],
}

# Heuristic topic tags for relevance
topic_kw = {
    'GW propagation/cosmology': [r'cosmolog', r'redshift\s*propagation', r'GW propagation', r'tensor mode'],
    'Lambda / vacuum / DE':     [r'dark energy', r'cosmological constant', r'\bLambda\b.{0,15}(CDM|vacuum)', r'27-decade'],
    'Lagrangian derivation':    [r'Lagrangian', r'master equation'],
    'Falsifiable prediction':   [r'falsif', r'sub-?mm', r'DESI', r'Euclid', r'CMB-?S4'],
    'KK / extra-dim':           [r'\bKK\b', r'extra dimension', r'Kaluza', r'compactification'],
}

for fn in papers:
    path = os.path.join('whitepapers', fn)
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()
    size = len(text)
    print('='*78)
    print(f'{fn}  [{size} bytes]')
    # Title line
    title = ''
    for line in text.split('\n')[:30]:
        if line.startswith('# '):
            title = line[2:].strip(); break
    print(f'  Title: {title[:90]}')
    # v5.78 signals already present?
    present = []
    for sig, pats in v578_signals.items():
        for p in pats:
            if re.search(p, text, re.IGNORECASE):
                present.append(sig); break
    print(f'  v5.78 signals already in body: {present if present else "(none)"}')
    # Topic tags
    tags = []
    for tag, pats in topic_kw.items():
        for p in pats:
            if re.search(p, text, re.IGNORECASE):
                tags.append(tag); break
    print(f'  Topic tags: {tags if tags else "(none)"}')
