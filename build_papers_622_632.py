"""build_papers_622_632.py -- Generate PDFs for PAPER_622-632 (Session 161)
BigBang Hypergraph Theory: Zero-Mass UA, 9D Wolfram, 26D Sculpting, Jets, Grants
    PAPER_622  UQFFZeroMassAetherVacuumGradientReformulationCalculator (#209)
    PAPER_623  UQFFNineDimensionalWolframForceTroadProjectionCalculator (#210)
    PAPER_624  UQFF26DSimultaneousGeometricInfinitySculptingCalculator (#211)
    PAPER_625  UQFFExoticPocketedShellQuantumFrequencyCalculator (#212)
    PAPER_626  UQFFM87JetNineDHypergraphPocketShellSimulationCalculator (#213)
    PAPER_627  UQFFCentaurusAKnottedJetVHEHypergraphCalculator (#214)
    PAPER_628  UQFFNGC6278DwarfGalaxyVoidPocketShellCalculator (#215)
    PAPER_629  UQFFMS073567421ClusterAGNJetVoidPocketCalculator (#216)
    PAPER_630  UQFFPerseusClusterIXPEXRayPolarizationJetCalculator (#217)
    PAPER_631  UQFFMultiSystemJetHypergraphComparisonCalculator (#218)
    PAPER_632  UQFFGrantProposalDatasetCompressionFrameworkCalculator (#219)
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
    'PAPER_622_UQFF_Zero_Mass_Aether_Vacuum_Gradient_Reformulation.md',
    'PAPER_623_UQFF_Nine_Dimensional_Wolfram_Force_Triad_Projection.md',
    'PAPER_624_UQFF_26D_Simultaneous_Geometric_Infinity_Sculpting.md',
    'PAPER_625_UQFF_Exotic_Pocketed_Shell_Quantum_Frequency_Events.md',
    'PAPER_626_UQFF_M87_Jet_9D_Hypergraph_Pocket_Shell_Simulation.md',
    'PAPER_627_UQFF_Centaurus_A_Knotted_Jet_VHE_Hypergraph.md',
    'PAPER_628_UQFF_NGC6278_Dwarf_Galaxy_Void_Pocket_Shell.md',
    'PAPER_629_UQFF_MS073567421_Cluster_AGN_Jet_Void_Pocket.md',
    'PAPER_630_UQFF_Perseus_Cluster_IXPE_XRay_Polarization_Jet_Solution.md',
    'PAPER_631_UQFF_Multi_System_Jet_Hypergraph_Comparison.md',
    'PAPER_632_UQFF_Grant_Proposal_Dataset_Compression_Framework.md',
]

# Whitepapers for this session are in repo root (not whitepapers/ subfolder)
search_dirs = [repo, repo / 'whitepapers']

print('=' * 78)
print('  UQFF PDF Generator -- PAPER_622-632 (Session 161)')
print('  Source: grok_share_6322ac199.txt (BigBang Hypergraph Theory, 27 Topics)')
print('  Physics: Zero-Mass UA / 9D Wolfram Force-Triad / 26D Simultaneous Sculpting')
print('           Exotic Pocket Shells / M87+CenA Jet Simulations')
print('           NGC6278/MS0735/Perseus Dataset Calculators / Multi-System / Grant')
print(f'  Output: {out_dir}')
print('=' * 78)

errors = []
for fname in targets:
    src = None
    for d in search_dirs:
        candidate = d / fname
        if candidate.exists():
            src = candidate
            break

    if src is None:
        print(f'  [SKIP]  {fname}  (not found in repo root or whitepapers/)')
        errors.append(f'NOT FOUND: {fname}')
        continue
    try:
        out  = genpdf.md_file_to_pdf(src, out_dir, styles)
        size = out.stat().st_size / 1024
        print(f'  [OK]    {out.name:<70s}  ({size:6.1f} KB)')
    except Exception as e:
        print(f'  [FAIL]  {fname}: {e}')
        errors.append(f'FAILED: {fname}: {e}')

print('=' * 78)
if errors:
    print(f'  {len(errors)} error(s):')
    for e in errors:
        print(f'    {e}')
    sys.exit(1)
else:
    print(f'  All {len(targets)} PDFs generated successfully.')
    print(f'  Total PDFs in pdf/: 637 prior + {len(targets)} = {637 + len(targets)}')
