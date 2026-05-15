import os, re, sys
sys.stdout.reconfigure(encoding='utf-8')

papers = [
    'PAPER_009b_Aether_String_TRZ_Damping_GW.md',
    'PAPER_015_Cosmological_Implications_UQFF_Modified_GW_Propagation.md',
    'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF.md',
    'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations.md',
    'PAPER_016b_White_Dwarf_Foreground_UQFF.md',
    'PAPER_017_Redshift_Corrections_z1_in_UQFF_GW_Propagation.md',
    'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA.md',
    'PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF.md',
    'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md',
    'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density.md',
]

SCAN = [
    ('beta_i', r'beta[_ ]?i|β[_ ]?i|\\beta_i'),
    ('F_TRZ', r'F[_ ]?TRZ|f_TRZ'),
    ('rho_SCm', r'rho[_ ]?SCm|ρ[_ ]?SCm|7\.09e-?37|7\.09\s*[x×]\s*10[\^\u207b\u00ad-]*37'),
    ('rho_UA', r'rho[_ ]?UA|ρ[_ ]?UA|7\.09e-?36|7\.09\s*[x×]\s*10[\^\u207b\u00ad-]*36'),
    ('SSq', r'\[SSq\]|SSq|0\.57'),
    ('kappa', r'\bkappa\b|κ'),
    ('Lambda_cosmo', r'\bLambda\b|\\Lambda|cosmological\s+constant|dark\s+energy|σ_?8|sigma_?8|w\(z\)'),
    ('ledger_27dec', r'27[- ]?decade|ledger|vacuum[- ]energy'),
    ('KK_xi', r'\bKK\b|Kaluza|extra dimension|13/3|sub-?mm'),
    ('GW_obs', r'LIGO|LISA|Euclid|DESI|CMB-?S4|pulsar timing|PTA|ringdown'),
]

print(f'{"FILE":58s} {"LINES":6s} | hooks')
print('-'*180)
for fn in papers:
    p = os.path.join('whitepapers', fn)
    if not os.path.exists(p):
        print(f'{fn}: MISSING'); continue
    with open(p, 'r', encoding='utf-8') as f: txt = f.read()
    hits = []
    for n, pat in SCAN:
        c = len(re.findall(pat, txt, re.IGNORECASE))
        if c: hits.append(f'{n}={c}')
    has = 'v5.78 Closure' in txt
    print(f'{fn:58s} {txt.count(chr(10)):6d} | block={"YES" if has else "no"} | {", ".join(hits)}')
    # Find unusual Unicode that might break pdflatex
    bad = sorted({hex(ord(c)) for c in txt if ord(c) > 0x2200 and ord(c) < 0x2300 and c not in '−–—…•·×÷'})
    if bad:
        # Look for the canonical problem chars
        problem = [c for c in '∝∈∉∀∃∇∂∞≈≃≅≠≤≥⊂⊃⊆⊇∪∩' if c in txt]
        if problem:
            print(f'  WARN unicode math chars: {problem}')
