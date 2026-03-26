"""build_papers_541_545.py -- Generate PDFs for PAPER_541-545 (Session 145)
Source: grok_share_22e7a1abb.txt — DPM-Proplyd Bidirectional Encompassment,
        UQFF Off-Diagonal Orion 4-Telescope Fit, Navier-Stokes Discrete Hypergraph
        Regularity, Yang-Mills DPM Gauge Field Mass Gap, Simultaneous Multi-Method
        Equivalence Merger Hub
"""
import sys, pathlib, importlib.util

spec = importlib.util.spec_from_file_location("genpdf", pathlib.Path(__file__).parent / "generate_pdfs.py")
genpdf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genpdf)

repo    = pathlib.Path(__file__).parent
out_dir = repo / 'pdf'
out_dir.mkdir(parents=True, exist_ok=True)
styles  = genpdf.make_styles()

targets = [
    # Session 145: grok_share_22e7a1abb.txt — CP4 #136-#140
    'PAPER_541_DPM_Proplyd_Bidirectional_Encompassment_Framework.md',
    'PAPER_542_UQFF_OffDiag_Proplyd_Orion_Four_Telescope_Fit.md',
    'PAPER_543_Navier_Stokes_Discrete_Hypergraph_Regularity_Proof.md',
    'PAPER_544_YangMills_DPM_Gauge_Field_Mass_Gap_Proof.md',
    'PAPER_545_Simultaneous_Multi_Method_Equivalence_Merger_Hub.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_541-545 (Session 145)')
print('  Source: grok_share_22e7a1abb.txt')
print('  Physics: DPM-Proplyd Bidirectional Encompassment,')
print('           UQFF Off-Diagonal Orion 4-Telescope Fit,')
print('           Navier-Stokes Discrete Hypergraph Regularity,')
print('           Yang-Mills DPM Gauge Field Mass Gap,')
print('           Simultaneous Multi-Method Equivalence Merger Hub')
print(f'  Output: {out_dir}')
print('=' * 70)

errors = []
for fname in targets:
    src = repo / 'whitepapers' / fname
    if not src.exists():
        print(f'  [SKIP]  {fname}  (not found)')
        errors.append(f'NOT FOUND: {fname}')
        continue
    try:
        out  = genpdf.md_file_to_pdf(src, out_dir, styles)
        size = out.stat().st_size / 1024
        print(f'  [OK]    {out.name:<65s}  ({size:6.1f} KB)')
    except Exception as e:
        print(f'  [FAIL]  {fname}: {e}')
        errors.append(f'FAILED: {fname}: {e}')

print('=' * 70)
if errors:
    print(f'  {len(errors)} error(s):')
    for e in errors:
        print(f'    {e}')
    sys.exit(1)
else:
    print(f'  All {len(targets)} PDFs generated successfully.')
    print(f'  Session 145 complete: PAPER_541-545; CP4 v5.05 #136-#140')
    sys.exit(0)
