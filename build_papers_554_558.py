"""build_papers_554_558.py -- Generate PDFs for PAPER_554-558 (Session 148)
Buoyancy-Stratified Factorial Geometry (BSFG) -- Complete Geometric System:
    PAPER_554  BSFG Riemann Curvature Aether Metric (#149)
    PAPER_555  BSFG Geodesic + Metric Compatibility (#150)
    PAPER_556  BSFG 26D Line Element + Factorial Compactification (#151)
    PAPER_557  BSFG Symmetry Group + Isometry Analysis (#152)
    PAPER_558  BSFG Unification Atlas Theorem Hub (#153)
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
    'PAPER_554_BSFG_Riemann_Curvature_Aether_Metric.md',
    'PAPER_555_BSFG_Geodesic_Metric_Compatibility.md',
    'PAPER_556_BSFG_26D_Line_Element_Factorial_Compactification.md',
    'PAPER_557_BSFG_Symmetry_Group_Isometry_Analysis.md',
    'PAPER_558_BSFG_Unification_Atlas_Theorem_Hub.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_554-558 (Session 148)')
print('  Physics: Buoyancy-Stratified Factorial Geometry (BSFG)')
print('           Riemann curvature, geodesic, 26D compactification,')
print('           SO(3)xU(1)^23 symmetry, VDS+DVP+BH26 atlas theorem')
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
    print(f'  Session 148 complete: PAPER_554-558; CP4 v5.08 #149-#153')
    print(f'  BSFG: R^r_0r0=1.56e-19 m^-2, G=SO(3)xU(1)^23, 26 generators=13+13')
    sys.exit(0)
