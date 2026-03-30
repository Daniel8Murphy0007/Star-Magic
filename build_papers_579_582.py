"""build_papers_579_582.py -- Generate PDFs for PAPER_579-582 (Session 156)
UQFF All-Forms / GW Amplitude / LQG Triple / String Planar:
    PAPER_579  UQFFAllFormsEvolutionCatalogueCalculator (#166)
    PAPER_580  UQFFGWAmplitudeLambdaCDMEmergenceCalculator (#167)
    PAPER_581  UQFFLQGLambdaCDMTripleSystemComparisonCalculator (#168)
    PAPER_582  StringGWPlanarFrequencyReboundDiskFormationCalculator (#169)
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
    'PAPER_579_UQFF_All_Forms_Evolution_Catalogue_Triadic_Solution.md',
    'PAPER_580_UQFF_GW_Amplitude_Lambda_CDM_Emergence.md',
    'PAPER_581_UQFF_LQG_LambdaCDM_Triple_System_QG_Comparison.md',
    'PAPER_582_String_GW_Planar_Frequency_Rebound_Disk_Formation.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_579-582 (Session 156)')
print('  Physics: UQFF All-Forms / GW Amplitude / LQG Triple / String Planar')
print('           PAPER_579: UQFF All Four Forms Evolution Catalogue + Triadic Solution')
print('           PAPER_580: GW Amplitude h=26!kQ/f27r + Lambda_CDM Dynamical Emergence')
print('           PAPER_581: UQFF vs LQG vs LambdaCDM Simultaneous Three-System QG Comparison')
print('           PAPER_582: String GW Planar Model + Universal Frequency Rebound + Disk Formation')
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
