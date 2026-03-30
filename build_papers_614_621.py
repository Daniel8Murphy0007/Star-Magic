"""build_papers_614_621.py -- Generate PDFs for PAPER_614-621 (Session 160)
26th-Order Complete Incorporation into UQFF Core Equations:
    PAPER_614  UQFFFUComplete26DProjectionOperatorCalculator (#201)
    PAPER_615  UQFFUg26DPolynomialDefectExpansionCalculator (#202)
    PAPER_616  UQFFUmDPMTimeDerivative26thOrderCalculator (#203)
    PAPER_617  UQFFSCmLaurentSeries26DExpansionCalculator (#204)
    PAPER_618  UQFFUbDensityGradient26thDerivativeCalculator (#205)
    PAPER_619  UQFFCompTensorFull26D13DCrossCalculator (#206)
    PAPER_620  UQFF3DIPODegree26TensorOverlayCalculator (#207)
    PAPER_621  UQFFPymanderSphere26DPyramidThreadCalculator (#208)
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
    'PAPER_614_UQFF_FU_Complete_26D_Projection_Operator.md',
    'PAPER_615_UQFF_Ug_26D_Polynomial_Defect_Expansion.md',
    'PAPER_616_UQFF_Um_DPM_Time_Derivative_26th_Order.md',
    'PAPER_617_UQFF_SCm_Laurent_Series_26D_Expansion.md',
    'PAPER_618_UQFF_Ub_Density_Gradient_26th_Derivative.md',
    'PAPER_619_UQFF_Comp_Tensor_Full_26D_13D_Cross.md',
    'PAPER_620_UQFF_3DIPO_Degree26_Tensor_Overlay.md',
    'PAPER_621_UQFF_Pymander_Sphere_26D_Pyramid_Thread.md',
]

print('=' * 74)
print('  UQFF PDF Generator -- PAPER_614-621 (Session 160)')
print('  Physics: 26th-Order Complete Incorporation into UQFF Core Equations')
print('           F_U 26D Projection / Ug Polynomial Defect / Um Time-Deriv 26')
print('           SCm Laurent Series / Ub Density 26th / Comp Tensor 3x3')
print('           3D-IPO Degree-26 Tensor / Pymander Sphere Pyramid Thread')
print(f'  Output: {out_dir}')
print('=' * 74)

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

print('=' * 74)
if errors:
    print(f'  {len(errors)} error(s):')
    for e in errors:
        print(f'    {e}')
    sys.exit(1)
else:
    print(f'  All {len(targets)} PDFs generated successfully.')
    print('  Session 160 complete -- CP4 v5.17, classes #201-208, PAPER_614-621')
    sys.exit(0)
