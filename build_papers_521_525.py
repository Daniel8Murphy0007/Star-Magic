"""build_papers_521_525.py -- Generate PDFs for PAPER_521-525 (Session 141)
Source: grok_share_3b6f26809.txt — BigBangHypergraphTheory continuation
5 papers: US Spectral Divisions, DPM Frequency Drive, Quantum Egg Sim, Plasma Orb Emergence, Hub
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
    # Session 141: grok_share_3b6f26809.txt — US Spectral / DPM / Proplyds
    'PAPER_521_Universal_Spectrum_Spectral_Divisions_ReRinging_BigBang_Vacuum_Gradient.md',
    'PAPER_522_DPM_Frequency_Drive_Ug1_Spectra_UQFF_Spectral_Tensor.md',
    'PAPER_523_Quantum_Egg_Frequency_Numerical_Simulation_Orion_Nebula_Validation.md',
    'PAPER_524_Plasma_Orb_Emergence_Threshold_Orion_Proplyd_Calibration.md',
    'PAPER_525_Session141_Hub_Universal_Spectrum_DPM_Quantum_Egg_Plasma_Orb_Proplyds.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_521-525 (Session 141)')
print('  Source: grok_share_3b6f26809.txt')
print('  Physics: US Spectral Divisions, DPM Frequency Drive,')
print('           Quantum Egg Sim, Plasma Orb Emergence, Hub')
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
