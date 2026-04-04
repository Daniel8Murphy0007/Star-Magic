#!/usr/bin/env python3
"""Fix missing H1 headings in PAPER_072, 114, 129, 137.
Also removes the corrupted duplicate block from PAPER_072."""

fixes = {
    'whitepapers/PAPER_072_Red_Dwarf_Reactor_Physics_UQFF.md': {
        'heading': '# PAPER_072: Red Dwarf Reactor (RDR) Physics — UQFF TRZ Factor Validation, COP > 1, and Plasma Temperature Agreement',
        'strip_lines': 9,   # Remove first 9 corrupt/duplicate lines
    },
    'whitepapers/PAPER_114_EP07_ParkerProbe_Heliosheath_Proof.md': {
        'heading': '# PAPER_114: Empirical Proof EP-07 — Parker Solar Probe Heliosheath: UQFF Ug2 Charge-Reactivity Field Validated',
        'strip_lines': 1,   # Remove leading blank line only
    },
    'whitepapers/PAPER_129_UQFF_Triadic_3C273_Jet_NegativeTime_N13.md': {
        'heading': '# PAPER_129: UQFF Triadic Mode Negative Time Discovery — 3C273 Asymmetric Quasar Jet: t_n < 0, R=130, N=13 Zero-Crossings',
        'strip_lines': 1,
    },
    'whitepapers/PAPER_137_UQFF_26QuantumLevels_EnergyLadder_E0to10n_Higgs_GalacticVacuum.md': {
        'heading': '# PAPER_137: UQFF Compressed Mode Universal Energy Ladder — 26 Quantum Levels (E_n = E_0 × 10ⁿ): Atomic to Galactic Vacuum',
        'strip_lines': 1,
    },
}

for path, cfg in fixes.items():
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Safety check: confirm the heading is not already there
    if lines[0].startswith('# PAPER_'):
        print(f'SKIP {path} - already has heading')
        continue

    # Remove leading blank/corrupt lines per config
    n = cfg['strip_lines']
    remaining = lines[n:]

    # Prepend the H1 heading with a blank line separator
    new_content = cfg['heading'] + '\n\n' + ''.join(remaining)

    with open(path, 'w', encoding='utf-8') as f:
        f.write(new_content)

    print(f'Fixed: {path}')

print('Done.')
