#!/usr/bin/env python3
"""Diagnose the 5 failing papers by generating .tex and showing error lines."""
import subprocess, tempfile, os, re

_BAD_CHARS = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f\ufffd\u202f]')

papers = {
    72: 'whitepapers/PAPER_072_Red_Dwarf_Reactor_Physics_UQFF.md',
    114: 'whitepapers/PAPER_114_EP07_ParkerProbe_Heliosheath_Proof.md',
    129: 'whitepapers/PAPER_129_UQFF_Triadic_3C273_Jet_NegativeTime_N13.md',
    137: 'whitepapers/PAPER_137_UQFF_26QuantumLevels_EnergyLadder_E0to10n_Higgs_GalacticVacuum.md',
    239: 'whitepapers/PAPER_239_UQFF_THz_Conduit_Shock_StarFormation_Forces.md',
}

# Error lines reported by xelatex from last run
error_lines = {72: 293, 114: 255, 129: 220, 137: 279, 239: 262}

BASE_CMD = [
    'pandoc','--pdf-engine=xelatex',
    '-V','geometry:a4paper,top=0.75in,bottom=0.75in,left=0.75in,right=0.75in',
    '-V','fontsize=11pt','-V','documentclass=article',
    '-H','pdf_header.tex',
    '--pdf-engine-opt=-interaction=nonstopmode',
    '--from=markdown-yaml_metadata_block-raw_tex+smart',
    '--standalone','--wrap=none','-t','latex',
]

for num, wp in papers.items():
    el = error_lines[num]
    with open(wp, encoding='utf-8') as f:
        text = f.read()
    text = _BAD_CHARS.sub('', text)
    tmp = tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', suffix='.md', delete=False)
    tmp.write(text)
    tmp.close()
    tex_tmp = tempfile.mktemp(suffix='.tex')
    subprocess.run(BASE_CMD + [tmp.name, '-o', tex_tmp], capture_output=True)
    if os.path.exists(tex_tmp):
        with open(tex_tmp, 'r', encoding='utf-8', errors='replace') as f:
            lines = f.readlines()
        print(f'=== PAPER_{num} around line {el} (total {len(lines)} tex lines) ===')
        start = max(0, el - 8)
        end = min(len(lines), el + 6)
        for i, line in enumerate(lines[start:end], start=start+1):
            marker = ' <-- ERROR' if i == el else ''
            safe = line.rstrip().encode('ascii', 'replace').decode()
            print(f'{i:4d}{marker}: {safe[:120]}')
        print()
        os.unlink(tex_tmp)
    os.unlink(tmp.name)
