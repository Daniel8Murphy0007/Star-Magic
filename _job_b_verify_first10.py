"""Verify the 10 updated papers: confirm closure block present, no content lost."""
import os, re, subprocess, sys
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

print(f'{"FILE":60s} {"BLOCK":6} {"DIFF_LINES":10} {"REF_PAPERS"}')
print('-'*120)

cross_refs = {}

for fn in papers:
    path = os.path.join('whitepapers', fn)
    pdf = os.path.join('pdf', fn.replace('.md','.pdf'))
    with open(path, 'r', encoding='utf-8') as f:
        text = f.read()
    has_block = '§v5.78 Closure' in text or 'v5.78 Closure' in text
    # git diff to confirm only additions (insertions count) and the previous content is preserved
    r = subprocess.run(['git','log','--all','--oneline','--follow','--numstat','HEAD~2..HEAD','--', path],
                       capture_output=True, text=True, encoding='utf-8', errors='replace')
    # numstat lines look like "85\t1\twhitepapers/..."
    add=rem=None
    for line in r.stdout.splitlines():
        m = re.match(r'^(\d+)\t(\d+)\t', line)
        if m:
            add = int(m.group(1)); rem = int(m.group(2)); break
    pdf_size = os.path.getsize(pdf) if os.path.exists(pdf) else 0
    # Find Cross-link references PAPER_XXX inside this paper (full body)
    refs = sorted(set(re.findall(r'PAPER_(\d{3,4})[a-z]?', text)))
    # Drop self
    self_id = re.match(r'PAPER_(\d+)', fn).group(1).lstrip('0') or '0'
    refs = [r for r in refs if r.lstrip('0') != self_id]
    cross_refs[fn] = refs
    print(f'{fn:60s} {"YES" if has_block else "NO":6} +{add}/-{rem if rem is not None else "?"}     pdf={pdf_size//1024}KB  refs={refs[:8]}{"..." if len(refs)>8 else ""}')

print()
print('='*80)
print('UNIQUE PAPERS CROSS-REFERENCED FROM THE 10 UPDATED PAPERS')
print('='*80)
all_refs = set()
for fn, refs in cross_refs.items():
    for r in refs:
        all_refs.add(r.lstrip('0') or '0')
# Filter out the 10 we already updated
already = {'1','2','3','4','5','6','7','8','9'}
referenced_not_yet_updated = sorted(int(r) for r in all_refs if r not in already and r.isdigit())
print(f'Total unique refs (excluding self): {len(all_refs)}')
print(f'References NOT in our updated batch: {referenced_not_yet_updated}')
