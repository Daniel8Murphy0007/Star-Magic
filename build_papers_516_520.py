"""build_papers_516_520.py -- Generate PDFs for PAPER_516-520 (Session 140)
Source: grok_share_0f5d4c91f2c.txt — BigBangHypergraphTheory recalculation
5 papers: DPM Shell Radiance, Negative Time Proof, DPM Forces, Shell Prototype, Hub
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
    # Session 140: grok_share_0f5d4c91f2c.txt — DPM Shell Radiance & Negative Time
    'PAPER_516_DPM_Layered_Shell_Energy_Radiance_Phase_Cascade.md',
    'PAPER_517_Negative_Time_Dilation_Proof_Spooky_Distance_Dual_Existence.md',
    'PAPER_518_DPM_Unified_Inertia_Centripetal_Centrifugal_Forces.md',
    'PAPER_519_Shell_Radiance_Prototype_Full_26D_Layer_Formulation.md',
    'PAPER_520_Session140_Hub_DPM_Shell_Radiance_Negative_Time_Forces.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_516-520 (Session 140)')
print('  Source: grok_share_0f5d4c91f2c.txt')
print('  Physics: DPM Shell Radiance, Negative Time, DPM-Unified Forces')
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
    print('  Session 140 complete: PAPER_516-520')
