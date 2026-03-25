"""build_papers_509_515.py -- Generate PDFs for PAPER_509-515 (Session 138)
Source module: source179.cpp (PI Co-Resonance Field, namespace SOURCE179)
7 papers: PCR field theory + 6 astrophysical validation targets
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
    # Session 138: source179.cpp — PI Co-Resonance Field + 6 astro targets
    'PAPER_509_PI_Co_Resonance_Field_Equations.md',
    'PAPER_510_GW150914_LIGO_BinaryBH_UQFF_Validation.md',
    'PAPER_511_PSR_J0437_SacredQuantumOrbit.md',
    'PAPER_512_Eta_Carinae_BuoyantGravity_PCR.md',
    'PAPER_513_NGC1277_Hypergraph_SpacetimeDimension.md',
    'PAPER_514_TON618_SacredTimePhase_Integral.md',
    'PAPER_515_TXS0506_IceCube_PICoSum_SpectralIndex.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_509-515 (Session 138)')
print('  Source: source179.cpp | Namespace: SOURCE179')
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
        print(f'  [OK]    {out.name:<60s}  ({size:6.1f} KB)')
    except Exception as e:
        print(f'  [FAIL]  {fname}: {e}')
        errors.append(f'FAILED: {fname}: {e}')

print('=' * 70)
if errors:
    print(f'  {len(errors)} error(s):')
    for e in errors: print(f'    {e}')
    sys.exit(1)
else:
    print(f'  All {len(targets)} PDFs generated successfully.')
print('=' * 70)
