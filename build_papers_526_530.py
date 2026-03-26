"""build_papers_526_530.py -- Generate PDFs for PAPER_526-530 (Session 142)
Source: grok_share_2515709ed.txt — BigBangHypergraphTheory Millennium proof set
5 papers: 3D-IPO Helical Overlay, Pymander Sphere, UQFF_comp Eigenvalue,
          Navier-Stokes UQFF Encompassment, Millennium Hub YM+Riemann+PvsNP
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
    # Session 142: grok_share_2515709ed.txt — Millennium proof set
    'PAPER_526_3D_IPO_Non_Linear_Three_Helix_Progression_Overlay.md',
    'PAPER_527_Pymander_Sphere_Six_Pyramid_Prob_order_Geometry.md',
    'PAPER_528_UQFF_comp_Spectral_Compression_Eigenvalue_Stability.md',
    'PAPER_529_Navier_Stokes_UQFF_Quasar_Jet_Regularity.md',
    'PAPER_530_Session142_Hub_Millennium_YangMills_Riemann_PvsNP_UQFF.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_526-530 (Session 142)')
print('  Source: grok_share_2515709ed.txt')
print('  Physics: 3D-IPO Helical Overlay, Pymander Sphere Prob_order,')
print('           UQFF_comp Eigenvalue, NS-UQFF Encompassment, Millennium Hub')
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
