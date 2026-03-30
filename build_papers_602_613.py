"""build_papers_602_613.py -- Generate PDFs for PAPER_602-613 (Session 159)
Cosmic Egg / 26D Egg Energy / Proto-H Shell / Factorial Bounds /
Inertia-Centripetal-Centrifugal 26D Shell Forces / Riemann Hypothesis /
Mayan Calendar Nuclei Epochs / Solar System Proplyd Legacy /
Probability of Order Partition / NASA ATP Grant Framework:
    PAPER_602  UQFFCosmicEggPreFertilizationEnergyCalculator (#189)
    PAPER_603  UQFF26DEggTotalEnergyCalculator (#190)
    PAPER_604  UQFFProtoHydrogenShellAlignmentCalculator (#191)
    PAPER_605  UQFF26thOrderFactorialBoundsCalculator (#192)
    PAPER_606  UQFFInertia26DShellForceCalculator (#193)
    PAPER_607  UQFFCentripetal26DShellCalculator (#194)
    PAPER_608  UQFFCentrifugal26DShellCalculator (#195)
    PAPER_609  UQFFRiemannHypothesisCriticalLineCalculator (#196)
    PAPER_610  UQFFMayanCalendarNucleiEpochCalculator (#197)
    PAPER_611  UQFFSolarSystemProplydLegacyCalculator (#198)
    PAPER_612  UQFFProbabilityOfOrderPartitionCalculator (#199)
    PAPER_613  UQFFNASAATPGrantFrameworkValidationCalculator (#200)
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
    'PAPER_602.md',
    'PAPER_603.md',
    'PAPER_604.md',
    'PAPER_605.md',
    'PAPER_606.md',
    'PAPER_607.md',
    'PAPER_608.md',
    'PAPER_609.md',
    'PAPER_610.md',
    'PAPER_611.md',
    'PAPER_612.md',
    'PAPER_613.md',
]

print('=' * 74)
print('  UQFF PDF Generator -- PAPER_602-613 (Session 159)')
print('  Physics: Cosmic Egg / 26D Egg Energy / Proto-H Shell Alignment')
print('           Factorial Bounds / 26D Shell Forces / Riemann Hypothesis')
print('           Mayan Epochs / Proplyd Legacy / Order Partition / ATP Grant')
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
    print('  Session 159 complete -- CP4 v5.16, classes #189-200, PAPER_602-613')
    sys.exit(0)
