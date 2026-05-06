files = [
    ('whitepapers/PAPER_002_GW190425_Mass_Gap_Interpretation.md', 800, 2000),
    ('whitepapers/PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md', 800, 2000),
    ('whitepapers/PAPER_034_Higgs_Kappa_t_Coupling_UQFF.md', 600, 1800),
]
for f, s, e in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    name = f.split('/')[-1]
    print('=== ' + name + ' ===')
    print(c[s:e])
    print()
