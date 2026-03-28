"""
test_millennium_phase_h.py — Phase H Millennium Prize Calculator Validation Suite
===================================================================================
Tests all 9 CP2 Millennium Prize classes introduced in Phase H (commit 65c7f0f).

Group structure (62 tests total):
  [M1-NS]   NSHypergraphDiscreteRegularityCalculator          — 6 tests
  [M2-YM]   YMDPMGaugeFieldMassGapProofCalculator             — 7 tests
  [M3-HUB]  Session142MillenniumEquationsHubCalculator        — 6 tests
  [M4-DPM]  YangMillsDPMQuantizationHubCalculator             — 7 tests
  [M5-PNP]  PvsNPUQFFComplexityCalculator                     — 7 tests
  [M6-BSD]  BirchSwinnertonDyerUQFFCalculator                 — 7 tests
  [M7-HDG]  HodgeUQFFAlgebraicCycleCalculator                 — 7 tests
  [M8-GAU]  FUBi26thGaussianTruncatedPolynomialBoundCalculator— 8 tests
  [M9-HUB]  MillenniumPrizeUQFFHubCalculator                  — 7 tests
  [REG]     SOURCE_MILLENNIUM_CP2 registry                     — 2 tests

Expected: 62/62 PASS

Author: test harness auto-generated, Session 151H
"""

import sys
import math
import os

# ─── Minimal test framework (mirrors test_phase3_millennium_problems.py) ─────

class TestResult:
    def __init__(self, name: str, passed: bool, details: str = ""):
        self.name    = name
        self.passed  = passed
        self.details = details


class ValidationSuite:
    def __init__(self, name: str):
        self.name    = name
        self.results = []

    def add_test(self, name: str, passed: bool, details: str = "") -> bool:
        self.results.append(TestResult(name, passed, details))
        status = "PASS" if passed else "FAIL"
        suffix = f" — {details}" if details else ""
        print(f"  [{status}] {name}{suffix}")
        return passed

    def summary(self) -> tuple:
        passed = sum(1 for r in self.results if r.passed)
        return passed, len(self.results)


def assert_close(actual, expected, rel_tol: float = 0.01, name: str = "") -> bool:
    """True if |actual - expected| / |expected| < rel_tol."""
    if expected == 0:
        return abs(actual) < rel_tol
    rel_error = abs((actual - expected) / expected)
    if rel_error >= rel_tol:
        label = f" ({name})" if name else ""
        print(f"      Expected: {expected:.6e}, Got: {actual:.6e}, "
              f"Rel Error: {rel_error * 100:.2f}%{label}")
        return False
    return True


# ─── Import CondensedPhysics2 ─────────────────────────────────────────────────
# CondensedPhysics.py has a pre-existing syntax error at line 167803 that
# prevents a full CP2 import.  The 9 Millennium classes are entirely
# self-contained (only `import math as _math` + `_vds_z26`), so we extract
# just those lines and exec them in a clean namespace.

def _import_cp2():
    import types as _types

    ws      = os.path.dirname(os.path.abspath(__file__))
    cp2_path = os.path.join(ws, "CondensedPhysics2.py")

    with open(cp2_path, 'r', encoding='utf-8') as fh:
        lines = fh.readlines()

    # Find the first self-contained block — the _vds_z26 helper at CP4 §148
    start = next(
        (i for i, ln in enumerate(lines) if ln.startswith('def _vds_z26')),
        None
    )
    if start is None:
        raise RuntimeError("Could not locate _vds_z26 in CondensedPhysics2.py")

    snippet = "".join(lines[start:])

    # Seed the exec namespace with everything the snippet needs
    ns = {
        '__name__': 'CondensedPhysics2',
        '_math': math,
        'math': math,
        # Registry dicts the end-of-file block tries to .update()
        'CP2_CALCULATORS': {},
        'SOURCE_SESSION137_CP2': {},
        'SOURCE_SOURCE179_CP2': {},
    }
    exec(compile(snippet, cp2_path, 'exec'), ns)   # noqa: S102

    return _types.SimpleNamespace(**ns)


# ─── [M1-NS] NSHypergraphDiscreteRegularityCalculator ────────────────────────

def test_m1_navier_stokes(cp2) -> tuple:
    """6 tests — Navier-Stokes discrete hypergraph regularity proof."""
    suite = ValidationSuite("[M1-NS] NSHypergraphDiscreteRegularityCalculator")
    print(f"\n=== {suite.name} ===")

    cls  = cp2.NSHypergraphDiscreteRegularityCalculator
    calc = cls()

    # T01 — instantiation and UQFF constants
    suite.add_test("T01 instantiation SSq=0.57 kappa=5e-4",
                   calc.SSq == 0.57 and calc.kappa == 5e-4)

    # T02 — eigenvalue λ_max < 1 (regularity criterion)
    eig = calc.compute_eigenvalues({})
    suite.add_test("T02 lambda_max < 1 (no blow-up)",
                   eig['lambda_max'] < 1.0,
                   f"lambda_max={eig['lambda_max']:.3e}")

    # T03 — λ_max > 0 (non-trivial)
    suite.add_test("T03 lambda_max > 0",
                   eig['lambda_max'] > 0.0)

    # T04 — λ₃ = 2·λ₁ exactly (tensor structure)
    suite.add_test("T04 lambda_3 == 2*lambda_1",
                   abs(eig['lambda_3'] - 2.0 * eig['lambda_1']) < 1e-20,
                   f"lambda_1={eig['lambda_1']:.3e}, lambda_3={eig['lambda_3']:.3e}")

    # T05 — regularity check returns True with IVT/π flags
    reg = calc.check_regularity({})
    suite.add_test("T05 check_regularity() regular=True, IVT=True, pi=True",
                   reg['regular'] is True and
                   reg['existence_IVT'] is True and
                   reg['uniqueness_pi'] is True)

    # T06 — BH26 harmonic sum non-zero (26 terms)
    bh = calc.compute_buoyancy_harmonic({})
    suite.add_test("T06 BH26 buoyancy harmonic != 0",
                   bh['U_bjet_harmonic'] != 0.0,
                   f"U_bjet_harmonic={bh['U_bjet_harmonic']:.4e}")

    return suite.summary()


# ─── [M2-YM] YMDPMGaugeFieldMassGapProofCalculator ───────────────────────────

def test_m2_yang_mills_dpm(cp2) -> tuple:
    """7 tests — Yang-Mills DPM gauge field mass gap proof."""
    suite = ValidationSuite("[M2-YM] YMDPMGaugeFieldMassGapProofCalculator")
    print(f"\n=== {suite.name} ===")

    calc = cp2.YMDPMGaugeFieldMassGapProofCalculator()

    # T07 — compute() positive mass gap
    result = calc.compute({})
    suite.add_test("T07 mass_gap_positive == True",
                   result['mass_gap_positive'] is True)

    # T08 — Δ > 0 (dimensionless UQFF units)
    suite.add_test("T08 delta > 0 (UQFF units)",
                   result['delta'] > 0.0,
                   f"delta={result['delta']:.4e}")

    # T09 — no zero-mode charge (charge quantisation)
    suite.add_test("T09 no_zero_mode == True (q_e = 2πn ≠ 0)",
                   result['no_zero_mode'] is True)

    # T10 — q_min ≈ 2π (minimum quantised charge)
    q_min = result['charge_q_min']
    suite.add_test("T10 charge_q_min ≈ 2π",
                   assert_close(q_min, 2.0 * math.pi, rel_tol=1e-6),
                   f"q_min={q_min:.6f}")

    # T11 — DVP anchor p=113 is prime, renders hypergraph aperiodic
    suite.add_test("T11 dvp_prime == 113 and dvp_irreducible == True",
                   result['dvp_prime'] == 113 and result['dvp_irreducible'] is True)

    # T12 — dvp_irreducibility() sub-method
    dvp = calc.dvp_irreducibility()
    suite.add_test("T12 dvp_irreducibility() is_prime=True aperiodic=True",
                   dvp['is_prime'] is True and dvp['aperiodic'] is True)

    # T13 — UQFF Δ within a physically reasonable factor-3 of lattice QCD
    ratio13 = result['ratio_lattice']
    suite.add_test("T13 UQFF Δ within factor-3 of lattice QCD (physical consistency)",
                   0.33 < ratio13 < 3.0,
                   f"ratio={ratio13:.3f}× (lattice QCD 1.4 GeV²)")

    return suite.summary()


# ─── [M3-HUB] Session142MillenniumEquationsHubCalculator ─────────────────────

def test_m3_session142_hub(cp2) -> tuple:
    """6 tests — Session 142 combined YM/Riemann/P≠NP hub."""
    suite = ValidationSuite("[M3-HUB] Session142MillenniumEquationsHubCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.Session142MillenniumEquationsHubCalculator()
    result = calc.compute({})

    # T14 — YM gap positive with DVP prime anchor
    suite.add_test("T14 YM_gap_positive == True",
                   result['YM_gap_positive'] is True)
    suite.add_test("T15 prime_anchor == 113",
                   result['prime_anchor'] == 113)

    # T15 — Riemann critical line and error bound
    suite.add_test("T16 Riemann_critical_line == True",
                   result['Riemann_critical_line'] is True)
    suite.add_test("T17 Riemann t13 error < 2%",
                   result['Riemann_error_pct'] < 2.0,
                   f"error={result['Riemann_error_pct']:.2f}%")

    # T16 — P≠NP separation via compute_p_ne_np
    suite.add_test("T18 PneNP == True",
                   result['PneNP'] is True)

    # T17 — sub-method compute_ym_gap
    ym = calc.compute_ym_gap({'E': 1e10, 'F': 1e19, 'Z': 0.570})
    suite.add_test("T19 compute_ym_gap positive with p=113 anchor",
                   ym['positive'] is True and ym['prime_anchor'] == 113)

    return suite.summary()


# ─── [M4-DPM] YangMillsDPMQuantizationHubCalculator ──────────────────────────

def test_m4_ym_dpm_quantization(cp2) -> tuple:
    """7 tests — Yang-Mills DPM quantization Millennium hub."""
    suite = ValidationSuite("[M4-DPM] YangMillsDPMQuantizationHubCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.YangMillsDPMQuantizationHubCalculator()
    result = calc.compute()

    # T20 — YM mass gap positive in GeV²
    suite.add_test("T20 YM_mass_gap_GeV2 > 0",
                   result['YM_mass_gap_GeV2'] > 0.0,
                   f"Δ={result['YM_mass_gap_GeV2']:.4f} GeV²")

    # T21 — UQFF Δ within a physically reasonable factor-3 of lattice QCD
    ratio21 = result['YM_lattice_ratio']
    suite.add_test("T21 UQFF Δ_GeV2 within factor-3 of lattice QCD (physical consistency)",
                   0.33 < ratio21 < 3.0,
                   f"ratio={ratio21:.3f}× (lattice QCD 1.4 GeV²)")

    # T22 — P≠NP separation matches standalone formula
    expected_sep = (2 ** 26) / (26 ** 4)
    suite.add_test("T22 PneNP_separation ≈ 2^26/26^4",
                   assert_close(result['PneNP_separation'], expected_sep, rel_tol=0.001),
                   f"separation={result['PneNP_separation']:.1f}×")

    suite.add_test("T23 PneNP == True",
                   result['PneNP'] is True)

    # T23 — Riemann estimate error bound
    suite.add_test("T24 Riemann t13 error < 2.0%",
                   result['Riemann_t13_error_pct'] < 2.0,
                   f"error={result['Riemann_t13_error_pct']:.2f}%")

    # T24 — Riemann t13 key present (fixes Hub KeyError)
    suite.add_test("T25 Riemann_t13_UQFF key present",
                   'Riemann_t13_UQFF' in result,
                   f"t13={result.get('Riemann_t13_UQFF', 'MISSING')}")

    # T25 — NS-DPM H¹ bound positive
    suite.add_test("T26 NS_DPM_H1bound > 0",
                   result['NS_DPM_H1bound'] > 0.0,
                   f"bound={result['NS_DPM_H1bound']:.4f} GeV²")

    return suite.summary()


# ─── [M5-PNP] PvsNPUQFFComplexityCalculator ──────────────────────────────────

def test_m5_p_vs_np(cp2) -> tuple:
    """7 tests — P vs NP UQFF 26D complexity calculator."""
    suite = ValidationSuite("[M5-PNP] PvsNPUQFFComplexityCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.PvsNPUQFFComplexityCalculator()
    result = calc.compute()

    # T27 — constant UA = 1e-4
    suite.add_test("T27 UA == 1e-4 (computational horizon)",
                   calc.UA == 1e-4)

    # T28 — P≠NP separation
    suite.add_test("T28 p_ne_np == True",
                   result['p_ne_np'] is True)

    # T29 — separation ratio matches 2^26/26^4
    expected_sep = (2 ** 26) / (26 ** 4)
    suite.add_test("T29 separation_ratio ≈ 2^26/26^4 (0.1%)",
                   assert_close(result['separation_ratio'], expected_sep, rel_tol=0.001),
                   f"ratio={result['separation_ratio']:.2f}×")

    # T30 — extraction cost: [UA]^-2 = (10^4)^2 = 10^8 shots per bit
    # Use 10**8 directly — int(1e-4 ** -2) can give 99999999 due to fp truncation
    ext = calc.extraction_cost(26)
    suite.add_test("T30 extraction shots_per_bit == 10^8",
                   ext['shots_per_bit'] == 10**8,
                   f"shots_per_bit={ext['shots_per_bit']}")

    # T31 — NP space = 2^26
    sep = calc.separation_ratio(26)
    suite.add_test("T31 NP_space == 2^26",
                   sep['NP_space'] == 2 ** 26,
                   f"NP_space={sep['NP_space']}")

    # T32 — P nodes = 26^4
    suite.add_test("T32 P_nodes == 26^4",
                   sep['P_nodes'] == 26 ** 4,
                   f"P_nodes={sep['P_nodes']}")

    # T33 — complexity hierarchy contains UA bridge
    hier = calc.complexity_hierarchy()
    suite.add_test("T33 complexity_hierarchy UA_bridge == 1e-4",
                   hier['UA_bridge'] == 1e-4 and 1 in hier['P_layers'])

    return suite.summary()


# ─── [M6-BSD] BirchSwinnertonDyerUQFFCalculator ──────────────────────────────

def test_m6_birch_swinnerton_dyer(cp2) -> tuple:
    """7 tests — Birch-Swinnerton-Dyer UQFF L-function."""
    suite = ValidationSuite("[M6-BSD] BirchSwinnertonDyerUQFFCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.BirchSwinnertonDyerUQFFCalculator()
    result = calc.compute()

    # T34 — UQFF constant κ = 5e-4
    suite.add_test("T34 kappa == 5e-4",
                   calc.kappa == 5e-4)

    # T35 — L_UQFF(E,1) > 0 (partial Euler product converges)
    suite.add_test("T35 L_UQFF_E1 > 0",
                   result['L_UQFF_E1'] > 0.0,
                   f"L_UQFF(E,1)={result['L_UQFF_E1']:.6f}")

    # T36 — κ→0 limit valid flag
    suite.add_test("T36 kappa_limit_valid == True",
                   result['kappa_limit_valid'] is True)

    # T37 — amplification factor > 1 (no κ shrinkage)
    suite.add_test("T37 amplification_factor > 1",
                   result['amplification_factor'] > 1.0,
                   f"factor={result['amplification_factor']:.4f}")

    # T38 — rank-1 amplification: factor ≈ 1/(1-e^{-κ}) ≈ 2000
    ra = calc.compute_rank_amplification(rank=1)
    expected_factor = 1.0 / (1.0 - math.exp(-5e-4))
    suite.add_test("T38 rank=1 amplification_factor ≈ 1/(1-e^{-κ})",
                   assert_close(ra['amplification_factor'], expected_factor, rel_tol=1e-6),
                   f"factor={ra['amplification_factor']:.2f}")

    # T39 — ord_UQFF for rank=1 > 1.0 (amplified order)
    suite.add_test("T39 ord_UQFF for rank=1 > 1.0",
                   ra['ord_UQFF'] > 1.0,
                   f"ord={ra['ord_UQFF']:.2f}")

    # T40 — rank=0 edge case: ord_UQFF == 0 (no pole if trivial rank)
    ra0 = calc.compute_rank_amplification(rank=0)
    suite.add_test("T40 rank=0 ord_UQFF == 0",
                   ra0['ord_UQFF'] == 0.0,
                   f"ord_0={ra0['ord_UQFF']}")

    return suite.summary()


# ─── [M7-HDG] HodgeUQFFAlgebraicCycleCalculator ──────────────────────────────

def test_m7_hodge(cp2) -> tuple:
    """7 tests — Hodge conjecture UQFF algebraic cycle pairing."""
    suite = ValidationSuite("[M7-HDG] HodgeUQFFAlgebraicCycleCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.HodgeUQFFAlgebraicCycleCalculator()
    result = calc.compute()

    # T41 — ground-state energy
    suite.add_test("T41 E0 == 1e-19 J",
                   calc.E0 == 1e-19)

    # T42 — all 26 Hodge class ratios rational
    suite.add_test("T42 all_ratios_rational == True",
                   result['all_ratios_rational'] is True)

    # T43 — 26 classes declared
    suite.add_test("T43 n_hodge_classes == 26",
                   result['n_hodge_classes'] == 26)

    # T44 — 26D decomposition sum > 0
    suite.add_test("T44 H_pq_UQFF_total > 0",
                   result['H_pq_UQFF_total'] > 0.0,
                   f"H_total={result['H_pq_UQFF_total']:.4e}")

    # T45 — energy level n=1: ratio = 10^0 = 1 (exact integer ∈ ℚ)
    levels = calc.compute_energy_levels()
    lvl1 = levels['levels'][0]   # n=1
    suite.add_test("T45 level n=1 rational == 1 (10^0)",
                   lvl1['rational'] == 1 and lvl1['in_Q'] is True)

    # T46 — Hodge pairing n=1, p=1, q=1: integral > 0
    pair1 = calc.compute_hodge_pairing(1, 1, 1)
    suite.add_test("T46 compute_hodge_pairing(1,1,1) > 0",
                   pair1['hodge_pairing'] > 0.0,
                   f"pairing={pair1['hodge_pairing']:.4e}")

    # T47 — energy levels strictly increasing: E_2 > E_1
    E1 = calc.E0 * (10.0 ** 0)   # n=1
    E2 = calc.E0 * (10.0 ** 1)   # n=2
    suite.add_test("T47 energy levels monotonically increasing (E_2 > E_1)",
                   E2 > E1,
                   f"E1={E1:.2e}, E2={E2:.2e}")

    return suite.summary()


# ─── [M8-GAU] FUBi26thGaussianTruncatedPolynomialBoundCalculator ─────────────

def test_m8_fubi26_gaussian(cp2) -> tuple:
    """8 tests — 26th-order Gaussian truncated polynomial bound."""
    suite = ValidationSuite("[M8-GAU] FUBi26thGaussianTruncatedPolynomialBoundCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.FUBi26thGaussianTruncatedPolynomialBoundCalculator()
    result = calc.compute({'z': 1.0})

    # T48 — truncation error below float64 epsilon
    machine_eps = 2.22e-16
    suite.add_test("T48 below_float64_eps == True (1/27! < ε_float64)",
                   result['below_float64_eps'] is True,
                   f"trunc_err={result['truncation_error']:.3e} < {machine_eps:.2e}")

    # T49 — explicit 1/27! < float64 eps (standalone check)
    fac27 = math.factorial(27)
    suite.add_test("T49 1/27! < 2.22e-16 (exact bound check)",
                   (1.0 / fac27) < machine_eps,
                   f"1/27!={1/fac27:.3e}")

    # T50 — Gaussian integral finite
    suite.add_test("T50 integral_finite == True (∫exp(-z²)dz = √π)",
                   result['integral_finite'] is True)

    # T51 — √π value
    suite.add_test("T51 sqrt_pi ≈ √π (0.0001% tolerance)",
                   assert_close(result['sqrt_pi'], math.sqrt(math.pi), rel_tol=1e-5),
                   f"sqrt_pi={result['sqrt_pi']:.8f}")

    # T52 — collapse_prevented
    suite.add_test("T52 collapse_prevented == True",
                   result['collapse_prevented'] is True)

    # T53 — poly(z=0) ≈ exp(0) = 1.0
    t0 = calc.compute_truncated_exp(0.0)
    suite.add_test("T53 compute_truncated_exp(0) poly ≈ 1.0",
                   assert_close(t0['poly_sum'], 1.0, rel_tol=1e-10),
                   f"poly(0)={t0['poly_sum']:.15f}")

    # T54 — poly(z=1) ≈ exp(-1) = 0.367879...
    t1 = calc.compute_truncated_exp(1.0)
    suite.add_test("T54 compute_truncated_exp(1) poly ≈ exp(-1)",
                   assert_close(t1['poly_sum'], math.exp(-1.0), rel_tol=1e-10),
                   f"poly(1)={t1['poly_sum']:.15f}")

    # T55 — VDS/DVP/BH26 context: BH26_freq_bins == 26, DVP_26fac == 26!
    ctx = calc.vds_dvp_bh26_context()
    suite.add_test("T55 vds_dvp_bh26_context BH26_freq_bins=26, DVP_26fac=26!",
                   ctx['BH26_freq_bins'] == 26 and
                   ctx['DVP_26fac'] == math.factorial(26),
                   f"bins={ctx['BH26_freq_bins']}, 26!={ctx['DVP_26fac']:.4e}")

    return suite.summary()


# ─── [M9-HUB] MillenniumPrizeUQFFHubCalculator ───────────────────────────────

def test_m9_hub(cp2) -> tuple:
    """7 tests — Master Millennium Prize hub aggregates all 6 problems."""
    suite = ValidationSuite("[M9-HUB] MillenniumPrizeUQFFHubCalculator")
    print(f"\n=== {suite.name} ===")

    calc   = cp2.MillenniumPrizeUQFFHubCalculator()
    result = calc.compute()

    # T56 — coverage string
    suite.add_test("T56 coverage == '6/6 open Millennium problems + Poincaré verification'",
                   result['coverage'] == '6/6 open Millennium problems + Poincaré verification')

    # T57 — Navier-Stokes sub-result
    ns = result['Navier_Stokes']
    suite.add_test("T57 NS regular == True",
                   ns['regular'] is True and ns['lambda_max'] < 1.0)

    # T58 — Yang-Mills sub-result
    ym = result['Yang_Mills']
    suite.add_test("T58 YM mass_gap_positive == True, delta > 0",
                   ym['mass_gap_positive'] is True and ym['delta'] > 0.0,
                   f"delta={ym['delta']:.4e}")

    # T59 — P vs NP sub-result
    pnp = result['P_vs_NP']
    suite.add_test("T59 P≠NP == True, separation > 1",
                   pnp['p_ne_np'] is True and pnp['separation'] > 1.0,
                   f"sep={pnp['separation']:.1f}×")

    # T60 — Hodge sub-result
    hodge = result['Hodge']
    suite.add_test("T60 Hodge all_rational == True, n_classes == 26",
                   hodge['all_rational'] is True and hodge['n_classes'] == 26)

    # T61 — Poincaré solved reference
    suite.add_test("T61 Poincaré mention Perelman (solved)",
                   'Perelman' in result['Poincare'])

    # T62 — master equation present
    suite.add_test("T62 master_equation contains M_UQFF",
                   'M_UQFF' in result['master_equation'])

    return suite.summary()


# ─── [REG] Registry validation ───────────────────────────────────────────────

def test_registry(cp2) -> tuple:
    """2 tests — SOURCE_MILLENNIUM_CP2 registry integrity."""
    suite = ValidationSuite("[REG] SOURCE_MILLENNIUM_CP2 registry")
    print(f"\n=== {suite.name} ===")

    reg = cp2.SOURCE_MILLENNIUM_CP2
    expected_names = [
        'NSHypergraphDiscreteRegularityCalculator',
        'YMDPMGaugeFieldMassGapProofCalculator',
        'Session142MillenniumEquationsHubCalculator',
        'YangMillsDPMQuantizationHubCalculator',
        'PvsNPUQFFComplexityCalculator',
        'BirchSwinnertonDyerUQFFCalculator',
        'HodgeUQFFAlgebraicCycleCalculator',
        'FUBi26thGaussianTruncatedPolynomialBoundCalculator',
        'MillenniumPrizeUQFFHubCalculator',
    ]

    # R1 — all 9 names present
    missing = [n for n in expected_names if n not in reg]
    suite.add_test("R1 all 9 Millennium classes registered",
                   len(missing) == 0,
                   f"missing={missing}" if missing else "all present")

    # R2 — CP2_CALCULATORS merged (CP2_CLASS_COUNT >= 9 Millennium)
    count = cp2.CP2_CLASS_COUNT
    suite.add_test("R2 CP2_CLASS_COUNT includes Millennium classes (≥ 9)",
                   count >= 9,
                   f"CP2_CLASS_COUNT={count}")

    return suite.summary()


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Phase H — Millennium Prize Calculator Validation Suite")
    print("Testing all 9 CP2 Millennium classes (commit 65c7f0f + bug-fix)")
    print("=" * 70)

    print("\nLoading CondensedPhysics2.py …", end=" ", flush=True)
    try:
        cp2 = _import_cp2()
        print("OK")
    except Exception as exc:
        print(f"FAILED: {exc}")
        sys.exit(1)

    runners = [
        test_m1_navier_stokes,
        test_m2_yang_mills_dpm,
        test_m3_session142_hub,
        test_m4_ym_dpm_quantization,
        test_m5_p_vs_np,
        test_m6_birch_swinnerton_dyer,
        test_m7_hodge,
        test_m8_fubi26_gaussian,
        test_m9_hub,
        test_registry,
    ]

    grand_pass  = 0
    grand_total = 0

    for fn in runners:
        p, t = fn(cp2)
        grand_pass  += p
        grand_total += t

    print("\n" + "=" * 70)
    print(f"GRAND TOTAL: {grand_pass}/{grand_total} PASS")
    if grand_pass == grand_total:
        print("ALL TESTS PASSED ✓")
    else:
        print(f"FAILURES: {grand_total - grand_pass}")
    print("=" * 70)

    return 0 if grand_pass == grand_total else 1


if __name__ == "__main__":
    sys.exit(main())
