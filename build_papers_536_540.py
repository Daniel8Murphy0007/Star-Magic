"""build_papers_536_540.py -- Generate PDFs for PAPER_536-540 (Session 144)
Source: grok_share_dbd886661cd.txt — DPM Split-Monopole MHD Proplyd Topology,
        Solar Body Proplyd Legacy 10-Body Table, UQFF Orion Triple-Telescope Fit,
        Extended Centripetal NS Residual, Yang-Mills DPM Quantization Millennium Hub
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
    # Session 144: grok_share_dbd886661cd.txt — CP4 #131-#135
    'PAPER_536_DPM_Split_Monopole_MHD_Proplyd_Topology.md',
    'PAPER_537_Solar_Body_Proplyd_Legacy_10_Body_Table.md',
    'PAPER_538_UQFF_Orion_Triple_Telescope_Encompassment_Fit.md',
    'PAPER_539_Extended_10_Body_Centripetal_Table_NS_Residual.md',
    'PAPER_540_YangMills_DPM_Quantization_Millennium_Hub.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_536-540 (Session 144)')
print('  Source: grok_share_dbd886661cd.txt')
print('  Physics: DPM Split-Monopole MHD, Solar Body Proplyd Legacy,')
print('           UQFF Orion Triple-Telescope, Extended Centripetal NS,')
print('           Yang-Mills DPM Quantization Millennium Hub')
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
    sys.exit(0)
