"""build_papers_559_562.py -- Generate PDFs for PAPER_559-562 (Session 149)
BSFG Open Questions Resolved:
    PAPER_559  BSFG Einstein Tensor & Self-Sourced Field Equations (#154)
    PAPER_560  BSFG Holonomy Group & Parallel Transport (#155)
    PAPER_561  BSFG Black Hole Horizon Solution (#156)
    PAPER_562  BSFG Bohr-Sommerfeld Aether Quantization (#157)
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
    'PAPER_559_BSFG_Einstein_Tensor_Field_Equations.md',
    'PAPER_560_BSFG_Holonomy_Group_Parallel_Transport.md',
    'PAPER_561_BSFG_BlackHole_Horizon_Solution.md',
    'PAPER_562_BSFG_BohrSommerfeld_Aether_Quantization.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_559-562 (Session 149)')
print('  Physics: BSFG Open Questions Resolved')
print('           Einstein tensor field eq (amp=1.8e4), holonomy SO+(3,1)xU(1)^22,')
print('           BH horizon r_h=0.23Rs T_H=3.4e-12K, Bohr-Sommerfeld r_cross=0.36AU')
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
    print(f'  Session 149 complete: PAPER_559-562; CP4 v5.09 #154-#157')
    print(f'  BSFG: amp=1.8e4, G_hol=SO+(3,1)xU(1)^22, r_h=1.62e8m, r_cross=0.36AU')
    sys.exit(0)
