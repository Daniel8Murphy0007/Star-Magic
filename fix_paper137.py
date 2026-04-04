#!/usr/bin/env python3
"""Fix PAPER_137 L108: move 10^{10} outside \\text{}."""
with open('whitepapers/PAPER_137_UQFF_26QuantumLevels_EnergyLadder_E0to10n_Higgs_GalacticVacuum.md', 'r', encoding='utf-8') as f:
    content = f.read()

# The replacement char in the file is U+FFFD
old = '$$E_{18} = 10^{-2} \\text{ J} = 62.4 \\text{ MeV} \\quad \\text{(within factor 2 of Higgs at 125.1 GeV \ufffd 10^{10} ? normalized)}$$'
new = '$$E_{18} = 10^{-2} \\text{ J} = 62.4 \\text{ MeV} \\quad \\text{(within factor 2 of Higgs:} \\sim 10^{10} \\text{ normalized)}$$'

print('OLD in content:', old in content)

if old in content:
    content = content.replace(old, new)
    with open('whitepapers/PAPER_137_UQFF_26QuantumLevels_EnergyLadder_E0to10n_Higgs_GalacticVacuum.md', 'w', encoding='utf-8') as f:
        f.write(content)
    print('Fixed PAPER_137')
else:
    # Try to find what's actually there
    idx = content.find('62.4 \\text{ MeV}')
    if idx != -1:
        ctx = content[idx-50:idx+200]
        print('Context:', repr(ctx))
