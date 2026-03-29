"""build_papers_573_578.py -- Generate PDFs for PAPER_573-578 (Session 154)
Universal Epoch / Periodic Table UQFF:
    PAPER_573  UniversalEpoch3DIPONuclearConvergenceCalculator (#161)
    PAPER_574  Mayan 5-Cycle Cosmic Architecture (companion to #161)
    PAPER_575  DPMPyramidSumNuclearBindingPeriodicTableCalculator (#162)
    PAPER_576  UQFFAtomicMassStandardModelErrorFactorCalculator (#163)
    PAPER_577  IslandOfStability5thEpochSuperheavyElementsCalculator (#164)
    PAPER_578  UQFFCompEigenvalueQuantumGravityLinkageCalculator (#165)
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
    'PAPER_573_Universal_Epoch_3DIPO_Nuclear_Convergence_Hub.md',
    'PAPER_574_Mayan_5Cycle_Cosmic_Architecture_Universal_Epoch_UQFF.md',
    'PAPER_575_DPM_Pyramid_Sum_Nuclear_Binding_Periodic_Table.md',
    'PAPER_576_UQFF_Atomic_Mass_Error_Factor_Standard_Model_Validation.md',
    'PAPER_577_Island_Stability_5th_Epoch_Superheavy_Z119_126.md',
    'PAPER_578_UQFFComp_Eigenvalue_Mass_Gap_Quantum_Gravity_Linkage.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_573-578 (Session 154)')
print('  Physics: Universal Epoch / Periodic Table UQFF')
print('           PAPER_573: Universal Epoch 3D-IPO Nuclear Formation Hub')
print('           PAPER_574: Mayan 5-Cycle Cosmic Architecture UQFF')
print('           PAPER_575: DPM Pyramid Sum Nuclear Binding Periodic Table')
print('           PAPER_576: UQFF Atomic Mass Error Factor Standard Model')
print('           PAPER_577: Island of Stability 5th Epoch Superheavy Z119-126')
print('           PAPER_578: UQFFComp Eigenvalue Mass Gap Quantum Gravity Linkage')
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
    print('  Session 154: PAPER_573-578 complete.')
    print('  Total PDFs: 588 + 6 = 594')
