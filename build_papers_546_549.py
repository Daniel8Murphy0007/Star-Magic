"""build_papers_546_549.py -- Generate PDFs for PAPER_546-549 (Session 146)
Source: grok_share_366dc393a37.txt — UgUb Boundary Overlap Simultaneous Displacement,
        Ug4 BH Tidal Time-reversal Stability, FUBi Collapse Prevention Eigenproof,
        Galaxy Merger UQFF vs Newton/Einstein Three-Method Hub
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
    # Session 146: grok_share_366dc393a37.txt — CP4 #141-#144
    'PAPER_546_UgUb_Boundary_Overlap_Simultaneous_Displacement.md',
    'PAPER_547_Ug4_BH_Tidal_Timereversal_Stability.md',
    'PAPER_548_FUBi_Universal_Buoyancy_Collapse_Prevention_Eigenproof.md',
    'PAPER_549_Galaxy_Merger_UQFF_vs_Newton_Einstein_ThreeMethod_Hub.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_546-549 (Session 146)')
print('  Source: grok_share_366dc393a37.txt')
print('  Physics: UgUb Boundary Overlap Simultaneous Displacement,')
print('           Ug4 BH Tidal Time-reversal Stability,')
print('           FUBi Collapse Prevention Eigenproof,')
print('           Galaxy Merger UQFF vs Newton/Einstein Three-Method Hub')
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
    print(f'  Session 146 complete: PAPER_546-549; CP4 v5.06 #141-#144')
    sys.exit(0)
