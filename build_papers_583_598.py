"""build_papers_583_598.py -- Generate PDFs for PAPER_583-598 (Session 157)
Six-Form UQFF / Collatz / Euler / Big Bang / Inflation / Maxwell26 / Dark Energy /
h-alpha-c-G derivations / BH Finite Bound / Sgr A* / QG Unification / t_neg / VDS-DVP-BH26:
    PAPER_583  UQFFSixFormSimultaneousSolverCalculator (#170)
    PAPER_584  UQFFCollatzConvergence26DCalculator (#171)
    PAPER_585  UQFFEulerEquationsInviscidProofCalculator (#172)
    PAPER_586  UQFFBigBangExpansionDynamicsCalculator (#173)
    PAPER_587  UQFFInflationaryEpochDetailsCalculator (#174)
    PAPER_588  UQFFMaxwellPowerLarge26thOrderCalculator (#175)
    PAPER_589  UQFFDarkEnergyVoidBuoyancyCalculator (#176)
    PAPER_590  UQFFPlanckConstantDerivedCalculator (#177)
    PAPER_591  UQFFFineStructureConstantDerivedCalculator (#178)
    PAPER_592  UQFFSpeedOfLightTriadEquilibriumCalculator (#179)
    PAPER_593  UQFFGravitationalConstantVoidCouplingCalculator (#180)
    PAPER_594  UQFFBlackHoleFiniteBoundCalculator (#181)
    PAPER_595  UQFFSgrAStarBoundApplicationCalculator (#182)
    PAPER_596  UQFFQuantumGravityUnificationCalculator (#183)
    PAPER_597  UQFFNegativeTimeDualExistenceCalculator (#184)
    PAPER_598  UQFFVDSDVPBH26IntegrationReferenceCalculator (#185)
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
    'PAPER_583_UQFF_Six_Form_Simultaneous_Solver.md',
    'PAPER_584_UQFF_Collatz_Convergence_26D.md',
    'PAPER_585_UQFF_Euler_Equations_Inviscid_Proof.md',
    'PAPER_586_UQFF_Big_Bang_Expansion_Dynamics.md',
    'PAPER_587_UQFF_Inflationary_Epoch_Details.md',
    'PAPER_588_UQFF_Maxwell_Power_Large_26th_Order.md',
    'PAPER_589_UQFF_Dark_Energy_Void_Buoyancy.md',
    'PAPER_590_UQFF_Planck_Constant_Derived.md',
    'PAPER_591_UQFF_Fine_Structure_Constant_Derived.md',
    'PAPER_592_UQFF_Speed_of_Light_Triad_Equilibrium.md',
    'PAPER_593_UQFF_Gravitational_Constant_Void_Coupling.md',
    'PAPER_594_UQFF_Black_Hole_Finite_Bound.md',
    'PAPER_595_UQFF_Sgr_A_Star_Bound_Application.md',
    'PAPER_596_UQFF_Quantum_Gravity_Unification.md',
    'PAPER_597_UQFF_Negative_Time_Dual_Existence.md',
    'PAPER_598_VDS_DVP_BH26_Integration_Reference.md',
]

print('=' * 74)
print('  UQFF PDF Generator -- PAPER_583-598 (Session 157)')
print('  Physics: 6-Form UQFF / Collatz / Euler / BB / Inflation / Maxwell26')
print('           Dark Energy / h,alpha,c,G derived / BH Bounds / Sgr A*')
print('           QG Unification / t_neg / VDS-DVP-BH26 Integration')
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
    sys.exit(0)
