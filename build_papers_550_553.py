"""build_papers_550_553.py -- Generate PDFs for PAPER_550-553 (Session 147)
Source: grok_share_b08cc4e3684.txt -- 26th-Order Polynomial Proofs:
        Um26D DPM Quantization Confinement (PAPER_550),
        Ug26D Factorial AntiCollapse Ug4 Split (PAPER_551),
        UQFFComp 26D Tensor OffDiag13 NS YM Hub (PAPER_552),
        FUBi 26th Gaussian Truncated Polynomial Bound (PAPER_553)
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
    # Session 147: grok_share_b08cc4e3684.txt -- CP4 #145-#148
    'PAPER_550_Um26D_Polynomial_DPM_Quantization_Confinement.md',
    'PAPER_551_Ug26D_Factorial_AntiCollapse_Ug4_Split.md',
    'PAPER_552_UQFFComp26D_Tensor_OffDiag13_NS_YM_Hub.md',
    'PAPER_553_FUBi26th_Gaussian_Polynomial_Bounded_Proof.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_550-553 (Session 147)')
print('  Source: grok_share_b08cc4e3684.txt')
print('  Physics: 26th-Order U_m DPM Quantization + Confinement (r_q=0.097AU),')
print('           26th-Order U_g Anti-Collapse + Ug4 13+13 Split,')
print('           UQFF_comp 26D Tensor Hub (T12=13!, YM gap=4.033e26),')
print('           F_U_Bi_i 26th Gaussian Polynomial Bounded Integral')
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
    print(f'  Session 147 complete: PAPER_550-553; CP4 v5.07 #145-#148')
    print(f'  26D manifold: r_q=0.097AU, rho_min=2.48e-30, YM gap=4.033e26')
    sys.exit(0)
