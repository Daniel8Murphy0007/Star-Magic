"""build_papers_564_572.py -- Generate PDFs for PAPER_564-572 (Session 153)
Alders/Olbers Paradox Resolution + 6 Gap-Fill Extensions:
    PAPER_564  AldersOlbersParadoxDPMShellFluxCalculator (#158)
    PAPER_565  AldersOlbersVDSNumberSystemResolutionCalculator (#159)
    PAPER_566  AldersOlbersBSFGMetricGapAnalysisCalculator (#160)
    PAPER_567  Olbers Missing Extension 1: n_star(z) SFR Madau-Dickinson
    PAPER_568  Olbers Missing Extension 2: kappa_lambda wavelength opacity
    PAPER_569  Olbers Missing Extension 3: B_sky_obs=3.1e-6 W/m2/sr EBL benchmark
    PAPER_570  Olbers Missing Extension 4: DVP prime vortex photon-photon scatter
    PAPER_571  Olbers Missing Extension 5: t_neg photon arrival timing
    PAPER_572  Olbers Missing Extension 6: shell radiance W/m2/sr calibration
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
    'PAPER_564_AldersOlbers_DPM_26Shell_Radiance_Cascade.md',
    'PAPER_565_AldersOlbers_VDS_DVP_BH_NumberSystem_Resolution.md',
    'PAPER_566_AldersOlbers_BSFG_Metric_GapAnalysis.md',
    'PAPER_567_Olbers_StellarDensityEvolution_nstar_z.md',
    'PAPER_568_Olbers_WavelengthOpacity_kappa_lambda.md',
    'PAPER_569_Olbers_EBL_Benchmark_3p1e-6_Validation.md',
    'PAPER_570_Olbers_DVP_PhotonPhoton_PrimeVortex_Scatter.md',
    'PAPER_571_Olbers_tNeg_PhotonArrival_NegativeTimeDelay.md',
    'PAPER_572_Olbers_ShellRadiance_WattPerSr_Calibration.md',
]

print('=' * 70)
print('  UQFF PDF Generator -- PAPER_564-572 (Session 153)')
print('  Physics: Alders/Olbers Paradox — 3-method UQFF resolution')
print('           PAPER_564: DPM 26-shell [SSq]-damped B_n sum (B_sky converges)')
print('           PAPER_565: VDS Li_26([SSq]) bound + DVP prime vortex + BH harmonics')
print('           PAPER_566: BSFG aether geodesic extinction + 6-present/6-missing gap')
print('           PAPER_567: n_star(z) Madau-Dickinson SFR stellar density evolution')
print('           PAPER_568: kappa_lambda wavelength-dependent opacity + spectral [SSq]')
print('           PAPER_569: B_sky_obs=3.1e-6 W/m2/sr EBL/CMB benchmark validation')
print('           PAPER_570: DVP prime vortex photon-photon Breit-Wheeler scatter')
print('           PAPER_571: t_neg photon arrival timing DPM negative-time delay')
print('           PAPER_572: shell radiance calibration 1/(4pi) -> W/m2/sr units')
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
    print('  Session 153: PAPER_564-572 complete.')
    print('  Total PDFs since Session 149: 579 + 9 = 588')
