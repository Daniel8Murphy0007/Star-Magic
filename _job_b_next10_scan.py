"""Dump first ~120 lines + scan key-constant + cross-ref summary for next 10 papers."""
import os, re, sys
sys.stdout.reconfigure(encoding='utf-8')

papers = [
    'PAPER_010_Post_Merger_Oscillations_Remnant_Mass_UQFF.md',
    'PAPER_010b_Time_Domain_Chirp_23Hz_UQFF.md',
    'PAPER_011_Stochastic_GW_Background_UQFF_Implications.md',
    'PAPER_011b_Amplitude_Reduction_Factor_UQFF.md',
    'PAPER_012_Eccentric_Binary_Circularization_UQFF.md',
    'PAPER_012b_GW150914_Waveform_Validation.md',
    'PAPER_013_Magnetar_Spin_Down_UQFF_Framework.md',
    'PAPER_013b_LISA_SMBH_Merger_Rate_UQFF.md',
    'PAPER_014_Primordial_Black_Holes_UQFF_Formation.md',
    'PAPER_014b_EMRI_Aether_Damping_UQFF.md',
]

# Tokens that indicate consumption of v5.78-now-derived constants
SCAN = [
    ('beta_i',          r'\bbeta[_ ]?i\b|β[_ ]?i|β_i|\\beta_i'),
    ('F_TRZ',           r'\bF[_ ]?TRZ\b|F_\{TRZ\}|f_TRZ'),
    ('rho_SCm',         r'rho[_ ]?SCm|ρ[_ ]?SCm|ρ_\{SCm\}|7\.09\s*[x×]\s*10\^?[-−]?37|7\.09e-?37'),
    ('rho_UA',          r'rho[_ ]?UA|ρ[_ ]?UA|7\.09\s*[x×]\s*10\^?[-−]?36|7\.09e-?36'),
    ('SSq_0.57',        r'\[SSq\]|0\.57'),
    ('kappa_5e-4',      r'\bkappa\b|κ|5\.?0?\s*[x×]\s*10\^?[-−]?4'),
    ('R26',             r'\bR26\b|26[- ]?layer|26[- ]?decade'),
    ('KK_tower',        r'\bKK\b|Kaluza|extra dimension|13/3'),
    ('damping_FU',      r'F_\{?U\}?|F_U_Bi|compute_FU|compute_compressed_MUGE'),
    ('Mexican_hat',     r'Mexican|broken[- ]symmetry|spontaneous'),
    ('P_suite',         r'\bP1[1-4]\b|P11|P12|P13|P14|LIGO O5|Euclid|DESI|CMB-?S4'),
]

print(f'{"FILE":58s} {"LINES":6} {"BYTES":7} | hooks-detected')
print('-'*180)
for fn in papers:
    p = os.path.join('whitepapers', fn)
    if not os.path.exists(p):
        print(f'{fn:58s}  MISSING')
        continue
    with open(p, 'r', encoding='utf-8') as f:
        txt = f.read()
    nlines = txt.count('\n')
    hits = []
    for name, pat in SCAN:
        c = len(re.findall(pat, txt, re.IGNORECASE))
        if c:
            hits.append(f'{name}={c}')
    has_block = 'v5.78 Closure' in txt
    print(f'{fn:58s} {nlines:6d} {len(txt):7d} | block={"YES" if has_block else "no"} | {", ".join(hits)}')

print()
print('=== TOPIC + REFERENCED PAPERS (Cross-link / See also lines) per file ===')
for fn in papers:
    p = os.path.join('whitepapers', fn)
    if not os.path.exists(p): continue
    with open(p, 'r', encoding='utf-8') as f: txt = f.read()
    # Grab title line and any "Cross-link" or "See also" lines
    title = re.search(r'^#\s+(.+)$', txt, re.MULTILINE)
    print(f'\n--- {fn} ---')
    print(f'  Title: {title.group(1) if title else "?"}')
    # Look for explicit paper cross-refs
    refs = sorted(set(re.findall(r'PAPER_(\d{3,4})[a-z]?', txt)))
    self_id = re.match(r'PAPER_(\d+)', fn).group(1)
    refs = [r for r in refs if r != self_id]
    print(f'  Cross-refs (3-digit): {[r for r in refs if len(r)==3]}')
    # Search for context strings that hint topic
    for kw in ['observation', 'predicts', 'damping', 'frequency', 'mass', 'energy']:
        m = re.search(rf'[^.\n]*\b{kw}\b[^.\n]*\.', txt, re.IGNORECASE)
        if m:
            snippet = m.group(0).strip()[:140]
            print(f'  [{kw}] {snippet}')
            break
