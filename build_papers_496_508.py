"""build_papers_496_508.py -- Generate PDFs for PAPER_496-508 (Sessions 136-137)"""
import sys, pathlib, importlib.util

spec = importlib.util.spec_from_file_location("genpdf", pathlib.Path(__file__).parent / "generate_pdfs.py")
genpdf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(genpdf)

repo    = pathlib.Path(__file__).parent
out_dir = repo / 'pdf'
out_dir.mkdir(parents=True, exist_ok=True)
styles  = genpdf.make_styles()

targets = [
    # Session 136: grok_share_1jkdsgv7
    'PAPER_496_DPM_DiPseudoMonopole_Full_Formulation.md',
    'PAPER_497_26D_Downward_Projection_Framework.md',
    'PAPER_498_3D_IPO_SCm_UA_Grinding_Sequence.md',
    'PAPER_499_Higgs_Inertial_Gradient_Shift_Marker.md',
    'PAPER_500_Proto_Hydrogen_26Shell_First_Atom.md',
    'PAPER_501_BBDT_Feynman_Globular_Clusters_1st_Epoch_BH.md',
    # Session 137: grok_share_84a767d3
    'PAPER_502_WSTP_Embedded_Kernel_Bridge.md',
    'PAPER_503_UQFF_Lagrangian_Wolfram_Export.md',
    'PAPER_504_WOLFRAM_TERM_AutoCollection_Framework.md',
    'PAPER_505_MSVC_Release_MaxCompress_Build_Profile.md',
    'PAPER_506_PI_Infinity_Decoder_Quantum_Mapping.md',
    'PAPER_507_Wolfram_Field_Unity_Engine_Hypergraph.md',
    'PAPER_508_Sacred_Time_Constants_Phase_Modulation.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_496-508 (Sessions 136-137)')
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
