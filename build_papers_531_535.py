"""build_papers_531_535.py -- Generate PDFs for PAPER_531-535 (Session 143)
Source: grok_share_fd81483544d.txt — VDS SCm Expansion, Quantum Plasma Orb,
        Solar System Proplyd DVP, Centripetal UQFF Encompassment Proof,
        VDS-DVP-BH Number Systems Unified Catalogue Hub
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
    # Session 143: grok_share_fd81483544d.txt — CP4 #126-#130
    'PAPER_531_BB_Hypergraph_Origin_VDS_SCm_Expansion.md',
    'PAPER_532_Quantum_Plasma_Orb_USorb_BH_Harmonic_Spectrum.md',
    'PAPER_533_Solar_System_Proplyd_DVP_Orbital_Quantization.md',
    'PAPER_534_Centripetal_Centrifugal_UQFF_Encompassment_Proof.md',
    'PAPER_535_VDS_DVP_BH_Number_Systems_Unified_Catalogue_Hub.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_531-535 (Session 143)')
print('  Source: grok_share_fd81483544d.txt')
print('  Physics: BB Hypergraph Origin VDS SCm, Quantum Plasma Orb,')
print('           Solar System Proplyd DVP, Centripetal UQFF Proof,')
print('           VDS-DVP-BH Number Systems Hub')
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
