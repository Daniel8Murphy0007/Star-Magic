import os
from collections import Counter
stems = [
 'PAPER_009b_Aether_String_TRZ_Damping_GW',
 'PAPER_015_Cosmological_Implications_UQFF_Modified_GW_Propagation',
 'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
 'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations',
 'PAPER_016b_White_Dwarf_Foreground_UQFF',
 'PAPER_017_Redshift_Corrections_z1_in_UQFF_GW_Propagation',
 'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA',
 'PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF',
 'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
 'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]
agg = Counter()
for s in stems:
    p = f'whitepapers/{s}.md'
    if not os.path.exists(p):
        print('MISSING', p); continue
    t = open(p, encoding='utf-8').read()
    c = Counter(ch for ch in t if ord(ch) > 127)
    print(f'== {s} MD non-ASCII unique={len(c)} total={sum(c.values())}')
    agg.update(c)
print('\n== AGGREGATE MD ACROSS 10 BATCH-3 SOURCES ==')
for k, v in agg.most_common(120):
    try:
        import unicodedata
        name = unicodedata.name(k, '?')
    except Exception:
        name = '?'
    print(f'  U+{ord(k):04X}  count={v:5}  {name}')
