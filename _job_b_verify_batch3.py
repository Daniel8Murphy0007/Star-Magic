"""Verify batch 3 papers: block present, diff additions-only, refs found."""
import os, re, subprocess, sys
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

print(f'{"FILE":62s} {"BLOCK":6} {"DIFF":10} {"PDF":8} REFS')
print('-'*130)
cross = {}
for fn in papers:
    p = os.path.join('whitepapers', fn)
    pdf = os.path.join('pdf', fn.replace('.md','.pdf'))
    with open(p,'r',encoding='utf-8') as f: txt=f.read()
    has = 'v5.78 Closure' in txt
    # git diff stat against parent of 67e6a25a (which was d2665820, but 009b parent is c43804fd...)
    # Use 67e6a25a~1 ... 67e6a25a
    r = subprocess.run(['git','log','--numstat','--format=','67e6a25a~1..67e6a25a','--', p],
                       capture_output=True, text=True, encoding='utf-8', errors='replace')
    add=rem='?'
    for line in r.stdout.splitlines():
        m=re.match(r'^(\d+)\t(\d+)\t',line)
        if m: add=m.group(1); rem=m.group(2); break
    sz = os.path.getsize(pdf)//1024 if os.path.exists(pdf) else 0
    refs = sorted(set(re.findall(r'PAPER_(\d{3,4})[a-z]?', txt)))
    self_id = re.match(r'PAPER_(\d+)', fn).group(1)
    refs = [x for x in refs if x != self_id]
    cross[fn] = refs
    print(f'{fn:62s} {"YES" if has else "NO":6} +{add}/-{rem:6s} {sz}KB    refs[3-digit]={[r for r in refs if len(r)==3]}')

# Aggregate unique cross-refs not in done set
DONE3 = {str(i).zfill(3) for i in range(1,22)} | {'008','009'}
# Actually done IDs (without variant): 1-9 + 9b + 10-14 + 15-21 = 1..21
done_ids = {str(i).zfill(3) for i in range(1,22)}
all_refs = set()
for fn,refs in cross.items():
    for r in refs: all_refs.add(r)
print()
print('Unique 3-digit refs NOT yet updated:')
not_done = sorted(int(r) for r in all_refs if len(r)==3 and r not in done_ids)
print(not_done)
print()
print('Unique 4-digit refs (Session 204-225 extensions / v5.78 papers):')
ext = sorted(int(r) for r in all_refs if len(r)==4)
print(ext)
