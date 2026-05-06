import re, os

targets = [
    'whitepapers/PAPER_001_GW170817_UQFF_Damping_Analysis.md',
    'whitepapers/PAPER_002_GW190425_Mass_Gap_Interpretation.md',
    'whitepapers/PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md',
    'whitepapers/PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF.md',
    'whitepapers/PAPER_013_Magnetar_Spin_Down_UQFF_Framework.md',
    'whitepapers/PAPER_014_Primordial_Black_Holes_UQFF_Formation.md',
    'whitepapers/PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF.md',
    'whitepapers/PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density.md',
    'whitepapers/PAPER_025_Dark_Matter_Direct_Detection_UQFF.md',
    'whitepapers/PAPER_034_Higgs_Kappa_t_Coupling_UQFF.md',
]
for f in targets:
    c = open(f, encoding='utf-8', errors='replace').read()
    title_m = re.search(r'title: "(.+?)"', c)
    title = title_m.group(1) if title_m else f
    has_refs = '## References' in c
    idx = c.rfind('## References')
    refs_snippet = c[idx:idx+500] if idx >= 0 else 'NO REFS SECTION'
    print(f'=== {f.split("/")[-1]} ===')
    print(f'Title: {title}')
    print(f'Has refs: {has_refs}')
    # first 600 chars of body after frontmatter
    body_start = c.find('\n# PAPER')
    if body_start >= 0:
        body = c[body_start:body_start+700]
        print('Body snippet:', body[:300])
    print('REFS:', refs_snippet[:250])
    print()
