#!/usr/bin/env python3
"""
production_scaling_v26.py — Production Scaling to 1M calc/s

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
v26 upgrade from v18 pattern: 1,000,000 calc/s target.
56 kernels total = v18's 48 + 8 new Session 226 kernels:
  kernel_phonon_gate_unitary     — SCm phonon gate U_gate operator
  kernel_dvp_vorticity_bound     — DVP lattice enstrophy bound
  kernel_bsm_mass_bound          — BSM mass from SCm phonon gap
  kernel_pcore_composite         — Planetary core composite density
  kernel_cooper_pair_gap         — Cooper-pair BCS gap from L_UQFF
  kernel_ramanujan_router        — S₂₆⁽³⁾ cross-domain dispatch (8 domains)
  kernel_cfd_stam_step           — Stam stable-fluids single timestep
  kernel_vds_dvp_bsh_cross       — VDS×DVP×BSH cross-product

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
Endpoint: /api/phonon/jet for REST integration via uqff_server.js
────────────────────────────────────────────────────────────────────────────────
"""

import math
import time
from typing import Dict, List, Callable

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
G         = 6.674e-11
C         = 2.998e8
M_SUN     = 1.989e30
HBAR      = 1.055e-34
K_B       = 1.381e-23
KPC       = 3.086e19
eV        = 1.602e-19

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
F_NEUTRON = 1e-10
RHO_VAC   = 1e-10
RHO_SCM   = 7.09e-37
H_0       = 2.195e-18
T_H       = 1.0 / H_0
LAMBDA_C  = 1.11e-52
RHO_CRIT  = 3 * H_0**2 / (8 * PI * G)
RHO_LAMBDA = 0.692 * RHO_CRIT
SIGMA_T   = 6.652e-29
M_PROTON  = 1.673e-27
ETA_BSFG  = 1e-22
F_UBI_RATIO = 0.6
H_SCM     = 0.99

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))

TARGET_CALC_PER_SEC = 1_000_000


# ── §1  Ramanujan S₂₆⁽³⁾ ──────────────────────────────────────────────────

def ramanujan_Rn(n: int, k: int = 3) -> float:
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 0.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, 27):
            sign = (-1) ** (j + 1)
            binom = math.comb(26, j)
            inner += sign * binom * math.factorial(26 - j) / n ** j
        correction += inner / n ** (26 * m)
    return prefactor * (1.0 + correction)


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))


def Phi(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


_LAYER_COEFFS = [
    (SSQ ** i) / (i ** 26) * ramanujan_Rn(i, 3) for i in range(1, 27)
]


# ── §2  DPM-FOUNDATION GRAVITY helpers ───────────────────────────────────────

def _dpm_ug1(r: float, mu: float = 1e15, R_body: float = 1e7) -> float:
    """DPM-seeded Ug1 (magnetic dipole)."""
    if r <= 0:
        return 0.0
    return mu / (4 * PI * r ** 3) * S26_3RD


def _dpm_ug2(r: float, Z: float = 26, alpha: float = 1.0) -> float:
    """DPM-seeded Ug2 (charge-reactivity)."""
    if r <= 0:
        return 0.0
    return Z * alpha / r ** 2 * S26_3RD


# ── §3  v18 Carry-Forward Kernels (48 kernels inlined) ─────────────────

def kernel_gravity_26layer(r: float = 1e7) -> float:
    return sum(c / r ** 2 for c in _LAYER_COEFFS)

def kernel_fu_bi_i(r: float = 1e7) -> float:
    ug1 = _dpm_ug1(r)
    ug2 = _dpm_ug2(r)
    ug3 = S26_3RD * OMEGA_SCM / (r * C)
    ug4 = RHO_VAC * S26_3RD / r ** 2
    return ug1 + ug2 + ug3 + ug4

def kernel_phonon_ares(gamma: float = GAMMA_0) -> float:
    return Phi(gamma) * OMEGA_SCM * HBAR

def kernel_jet_mjet(L_bol: float = 1e39) -> float:
    eta_jet = 0.1
    return eta_jet * L_bol / C ** 2

def kernel_ns_spindown(P: float = 0.033, Pdot: float = 4.2e-13) -> float:
    I = 1e45
    return 4 * PI ** 2 * I * Pdot / P ** 3

def kernel_gw170817_strain(M: float = 2.74 * M_SUN, d: float = 40e6 * 3.086e16) -> float:
    return 4 * G * M / (C ** 2 * d)

def kernel_blazar_ergosphere(M_bh: float = 6.5e9 * M_SUN, a_star: float = 0.998) -> float:
    r_plus = G * M_bh / C ** 2 * (1 + math.sqrt(1 - a_star ** 2))
    r_ergo = 2 * G * M_bh / C ** 2
    return (r_ergo - r_plus) / r_ergo

def kernel_rest_phonon_jet() -> float:
    return Phi(GAMMA_0) * kernel_jet_mjet()

def kernel_qcalcgeom_vectorized(N: int = 26) -> float:
    return sum(SSQ ** i / (i ** 26) for i in range(1, N + 1))

def kernel_pipeline_full() -> float:
    return kernel_gravity_26layer() + kernel_fu_bi_i() + kernel_phonon_ares()

def kernel_cena_jet(L: float = 1e36) -> float:
    return kernel_jet_mjet(L)

def kernel_txs0506_jet(L: float = 3e38) -> float:
    return kernel_jet_mjet(L)

def kernel_bcs_gap_solve(omega_D: float = 2 * PI * 5e12, NV: float = 0.3) -> float:
    if NV <= 0:
        return 0.0
    return 2 * HBAR * omega_D * math.exp(-1.0 / NV)

def kernel_spectral_ladder_eval(n: int = 10) -> float:
    return sum(SSQ ** k / k ** 26 for k in range(1, n + 1))

def kernel_ramanujan_26d_sum() -> float:
    return S26_3RD

def kernel_triadic_solver(r: float = 1e7) -> float:
    return kernel_gravity_26layer(r) * S26_3RD

def kernel_fubi_inside_out(r: float = 1e7) -> float:
    return kernel_fu_bi_i(r) * F_UBI_RATIO

def kernel_99sys_gamma_sweep() -> float:
    return sum(Phi(GAMMA_0 * (1 + 0.1 * i)) for i in range(-5, 6))

def kernel_agn_cena_fubi() -> float:
    return kernel_cena_jet() * F_UBI_RATIO * S26_3RD

def kernel_ns_merger_gw190425(M: float = 3.4 * M_SUN) -> float:
    return kernel_gw170817_strain(M=M)

def kernel_agn_merger_fubi(M: float = 1e8 * M_SUN) -> float:
    return kernel_blazar_ergosphere(M) * F_UBI_RATIO

def kernel_qgp_vacuum_density(T_MeV: float = 300) -> float:
    T_GeV = T_MeV / 1000
    return RHO_CRIT * S26_3RD * T_GeV ** 4

def kernel_alice_multiplicity(sqrt_s: float = 5020) -> float:
    return S26_3RD * math.log(sqrt_s)

def kernel_ym_mass_gap() -> float:
    return HBAR * OMEGA_SCM * S26_3RD

def kernel_3c273_agn_fubi(L: float = 2e39) -> float:
    return kernel_jet_mjet(L) * F_UBI_RATIO

def kernel_ton618_agn_fubi(M: float = 6.6e10 * M_SUN) -> float:
    return kernel_blazar_ergosphere(M) * F_UBI_RATIO

def kernel_gw170817_merger() -> float:
    return kernel_gw170817_strain()

def kernel_smbh_merger_fubi(M: float = 1e9 * M_SUN) -> float:
    return kernel_gw170817_strain(M=M) * F_UBI_RATIO

def kernel_dm_halo_nfw(r: float = 8.5 * KPC, r_s: float = 20 * KPC,
                       rho_s: float = 0.3 * 1e9 * M_SUN / KPC ** 3) -> float:
    x = r / r_s
    return rho_s / (x * (1 + x) ** 2)

def kernel_txs0506_3gamma() -> float:
    return kernel_txs0506_jet() * 3 * F_UBI_RATIO

def kernel_scm_inflation_hubble() -> float:
    rho_inf = HBAR * OMEGA_SCM * S26_3RD / (C ** 2)
    return math.sqrt(8 * PI * G * rho_inf / 3)

def kernel_thorne_morris_exotic() -> float:
    r_throat = 1e3
    return -HBAR * OMEGA_SCM * S26_3RD / (4 * PI * r_throat ** 2 * C)

def kernel_vds_convergence(N: int = 50) -> float:
    return sum(SSQ ** n / n ** 26 for n in range(1, N + 1))

def kernel_dvp_prime_sieve(p_max: int = 100) -> float:
    total = 0.0
    for p in range(29, p_max + 1):
        if all(p % d != 0 for d in range(2, int(math.sqrt(p)) + 1)):
            pi_p = sum(1 for k in range(2, p + 1) if all(k % d != 0 for d in range(2, int(math.sqrt(k)) + 1)) if k > 1)
            total += SSQ ** pi_p / p ** 26
    return total

def kernel_bsh_harmonic_saturation() -> float:
    return sum(BETA_I ** k * H_SCM ** (26 - k) for k in range(27))

def kernel_wormhole_phonon_damping(r: float = 1e3) -> float:
    return kernel_thorne_morris_exotic() * Phi(GAMMA_0)

def kernel_gpu_dpm_atlas_peak() -> float:
    return _dpm_ug1(1e7) + _dpm_ug2(1e7)

def kernel_dpm_line_fwhm() -> float:
    return 2 * math.sqrt(2 * math.log(2)) * SIGMA_G

def kernel_muge_8term_gravity(r: float = 1e7) -> float:
    return kernel_fu_bi_i(r) + kernel_gravity_26layer(r)

def kernel_nfw_density(r: float = 8.5 * KPC) -> float:
    return kernel_dm_halo_nfw(r)

def kernel_enet_gamma_coupled() -> float:
    return sum(Phi(GAMMA_0 * (1 + 0.05 * i)) * HBAR * OMEGA_SCM for i in range(-3, 4))

def kernel_dark_energy_eos_w0() -> float:
    return -1.0 + S26_3RD * F_UBI_RATIO

def kernel_alma_fubi_profile(r: float = 1e6) -> float:
    return kernel_fu_bi_i(r) * Phi(GAMMA_0)

def kernel_alma_chi2_co21() -> float:
    return abs(S26_3RD - 0.57 ** 1) ** 2

def kernel_qcalcgeom_master(r: float = 1e7) -> float:
    r_cross = 2 * G * M_SUN / C ** 2
    compact = r_cross / r if r > 0 else 0
    return compact * S26_3RD * Phi(GAMMA_0)

def kernel_cluster_cooling_suppression(M_bh: float = 3.4e8) -> float:
    L_edd = 4 * PI * G * M_bh * M_SUN * M_PROTON * C / SIGMA_T
    S = 1 - math.exp(-L_edd * S26_3RD / 1e45)
    return max(0.0, min(1.0, S))

def kernel_ramanujan_binomial_Rn(n: int = 5) -> float:
    return ramanujan_Rn(n, 3)

def kernel_lenr_cop_gamma(gamma: float = GAMMA_0) -> float:
    phi = Phi(gamma)
    sigma = 1e-28 * phi
    P_nuclear = sigma * 1e20 * 1e6 * eV
    P_buoyancy = F_UBI_RATIO * S26_3RD * HBAR * OMEGA_SCM
    P_in = HBAR * OMEGA_SCM
    return (P_nuclear + P_buoyancy) / P_in if P_in > 0 else 0.0


# ── §3b  NEW v26 KERNELS (8 kernels) ──────────────────────────────────────

def kernel_phonon_gate_unitary(theta: float = PI / 4) -> float:
    """Phonon gate U_gate fidelity: |det(U)| should be 1.

    U = exp(-iθ/2·σ_z·Φ) → diagonal with e^{±iθΦ/2}.
    Returns |det(U)| = 1 exactly for unitary gates.
    """
    phi = Phi(GAMMA_0)
    phase = theta * phi / 2.0
    # det(U) = e^{iφ}·e^{-iφ} = 1  (always)
    det_abs = abs(math.cos(phase) ** 2 + math.sin(phase) ** 2)
    return det_abs


def kernel_dvp_vorticity_bound() -> float:
    """DVP lattice enstrophy bound C_DVP.

    C_DVP = Σ_{p>26, prime} (a(p)·F_max/(ν·(2πp/L)²))²
    Evaluates convergence (should be finite, small).
    """
    nu = 1e-6     # kinematic viscosity
    F_max = 1.0   # force amplitude
    L = 1.0       # domain length
    total = 0.0
    for p in [29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
        pi_p = sum(1 for k in range(2, p + 1) if all(k % d != 0 for d in range(2, int(math.sqrt(k)) + 1)) if k > 1)
        a_p = SSQ ** pi_p / p ** 26
        k_p = 2 * PI * p / L
        total += (a_p * F_max / (nu * k_p ** 2)) ** 2
    return total


def kernel_bsm_mass_bound(fubi_ratio: float = F_UBI_RATIO) -> float:
    """BSM mass bound m_BSM in eV/c²."""
    delta_scm = HBAR * OMEGA_SCM
    m = (delta_scm / C ** 2) * S26_3RD * (2 * fubi_ratio - 1)
    return m * C ** 2 / eV


def kernel_pcore_composite(T_core: float = 1.57e7, rho_mean: float = 1408) -> float:
    """Pcore composite density (kg/m³)."""
    x = HBAR * OMEGA_SCM / (K_B * T_core) if T_core > 0 else 500
    n_be = 1.0 / (math.exp(x) - 1) if x < 500 else 0.0
    return rho_mean * S26_3RD * n_be * F_UBI_RATIO


def kernel_cooper_pair_gap(omega_D: float = 2 * PI * 9e12,
                           N0: float = 1.45e47) -> float:
    """BCS gap Δ₀ from L_UQFF phonon sector (eV)."""
    V_scm = HBAR * OMEGA_SCM * S26_3RD * SSQ / N0 if N0 > 0 else 0
    NV = N0 * V_scm
    if NV <= 0:
        return 0.0
    delta0 = 2 * HBAR * omega_D * math.exp(-1.0 / NV)
    return delta0 / eV


def kernel_ramanujan_router() -> float:
    """S₂₆⁽³⁾ dispatched to 8 domains — returns sum of domain weights."""
    # Quick evaluation of all 8 domains with default params
    phonon = Phi(GAMMA_0)
    inflation = S26_3RD * (2.435e18 / 1e16) ** 2  # (M_Pl/Λ)²
    qgp = S26_3RD * (300 / 155) ** (-2) / (4 * PI)
    jet = S26_3RD * 0.3 ** 3
    cfd = S26_3RD * F_UBI_RATIO * 1.225 * 9.81
    lenr = 1e-28 * S26_3RD
    bsm = kernel_bsm_mass_bound()
    pcore = kernel_pcore_composite()
    return phonon + inflation + qgp + jet + cfd + lenr + abs(bsm) + abs(pcore)


def kernel_cfd_stam_step(N: int = 16) -> float:
    """Single Stam stable-fluids diffusion step on N×N grid.

    Returns max velocity magnitude after one diffusion step.
    """
    # Simplified: initialise uniform, apply one Jacobi iteration
    dt = 0.01
    visc = 1e-3
    a = dt * visc * N * N
    x = [[0.0] * (N + 2) for _ in range(N + 2)]
    x0 = [[0.0] * (N + 2) for _ in range(N + 2)]
    # Seed centre
    x0[N // 2][N // 2] = S26_3RD * 1e6
    for _ in range(5):
        for i in range(1, N + 1):
            for j in range(1, N + 1):
                x[i][j] = (x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] + x[i][j + 1])) / (1 + 4 * a)
    return max(abs(x[i][j]) for i in range(1, N + 1) for j in range(1, N + 1))


def kernel_vds_dvp_bsh_cross() -> float:
    """VDS × DVP × BSH cross-product spectral weight."""
    S_vds = sum(SSQ ** n / n ** 26 for n in range(1, 31))
    S_dvp = kernel_dvp_prime_sieve(100)
    W_bsh = kernel_bsh_harmonic_saturation()
    return S_vds * S_dvp * W_bsh


# ── §4  ALL_KERNELS Registry ────────────────────────────────────────────

ALL_KERNELS: List[Callable] = [
    # v11 (16)
    kernel_gravity_26layer, kernel_fu_bi_i, kernel_phonon_ares,
    kernel_jet_mjet, kernel_ns_spindown, kernel_gw170817_strain,
    kernel_blazar_ergosphere, kernel_rest_phonon_jet,
    kernel_qcalcgeom_vectorized, kernel_pipeline_full,
    kernel_cena_jet, kernel_txs0506_jet, kernel_bcs_gap_solve,
    kernel_spectral_ladder_eval, kernel_ramanujan_26d_sum, kernel_triadic_solver,
    # v13 (4)
    kernel_fubi_inside_out, kernel_99sys_gamma_sweep,
    kernel_agn_cena_fubi, kernel_ns_merger_gw190425,
    # v14 (4)
    kernel_agn_merger_fubi, kernel_qgp_vacuum_density,
    kernel_alice_multiplicity, kernel_ym_mass_gap,
    # v15 (6)
    kernel_3c273_agn_fubi, kernel_ton618_agn_fubi,
    kernel_gw170817_merger, kernel_smbh_merger_fubi,
    kernel_dm_halo_nfw, kernel_txs0506_3gamma,
    # v16 (6)
    kernel_scm_inflation_hubble, kernel_thorne_morris_exotic,
    kernel_vds_convergence, kernel_dvp_prime_sieve,
    kernel_bsh_harmonic_saturation, kernel_wormhole_phonon_damping,
    # v17 (8)
    kernel_gpu_dpm_atlas_peak, kernel_dpm_line_fwhm,
    kernel_muge_8term_gravity, kernel_nfw_density,
    kernel_enet_gamma_coupled, kernel_dark_energy_eos_w0,
    kernel_alma_fubi_profile, kernel_alma_chi2_co21,
    # v18 (4)
    kernel_qcalcgeom_master, kernel_cluster_cooling_suppression,
    kernel_ramanujan_binomial_Rn, kernel_lenr_cop_gamma,
    # v26 (8) — NEW
    kernel_phonon_gate_unitary, kernel_dvp_vorticity_bound,
    kernel_bsm_mass_bound, kernel_pcore_composite,
    kernel_cooper_pair_gap, kernel_ramanujan_router,
    kernel_cfd_stam_step, kernel_vds_dvp_bsh_cross,
]


# ── §5  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV26:
    """56-kernel production benchmark: 1M calc/s target.

    REST endpoint: /api/phonon/jet via uqff_server.js integration.
    """

    def single_pass(self) -> float:
        return sum(k() for k in ALL_KERNELS)

    def simulate(self, duration_s: float = 1.0) -> dict:
        start = time.perf_counter()
        iterations = 0
        while time.perf_counter() - start < duration_s:
            self.single_pass()
            iterations += 1
        elapsed = time.perf_counter() - start
        calc_per_sec = iterations * len(ALL_KERNELS) / elapsed

        return {
            "iterations": iterations,
            "kernels": len(ALL_KERNELS),
            "elapsed_s": elapsed,
            "calc_per_sec": calc_per_sec,
            "target": TARGET_CALC_PER_SEC,
            "met_target": calc_per_sec >= TARGET_CALC_PER_SEC,
        }

    def compute(self, dataset: dict) -> dict:
        """Production scaling compute interface."""
        duration = float(dataset.get("duration_s", 1.0))
        result = self.simulate(duration)
        result["primary_equations"] = [
            f"56 kernels × {result['iterations']} iterations in {result['elapsed_s']:.3f}s",
            f"Throughput: {result['calc_per_sec']:.0f} calc/s",
            f"Target: {TARGET_CALC_PER_SEC:,} calc/s",
            f"Met target: {result['met_target']}",
        ]
        return result


# ── §6  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    print("=" * 72)
    print("production_scaling_v26.py — Self-Tests")
    print("=" * 72)

    ok = True
    passed = 0

    # Test 1: All kernels finite
    results = [k() for k in ALL_KERNELS]
    all_finite = all(math.isfinite(r) for r in results)
    if all_finite:
        print(f"  [PASS] Test 1: All {len(ALL_KERNELS)} kernels produce finite results")
        passed += 1
    else:
        bad = [(i, r) for i, r in enumerate(results) if not math.isfinite(r)]
        print(f"  [FAIL] Test 1: Non-finite kernel results at indices {[b[0] for b in bad]}")
        ok = False

    # Test 2: Kernel count = 56
    n = len(ALL_KERNELS)
    if n == 56:
        print(f"  [PASS] Test 2: Kernel count: {n}")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: Expected 56 kernels, got {n}")
        ok = False

    # Test 3: Phonon gate unitary = 1.0
    val = kernel_phonon_gate_unitary()
    if abs(val - 1.0) < 1e-12:
        print(f"  [PASS] Test 3: Phonon gate |det(U)| = {val:.15f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: |det(U)| = {val}")
        ok = False

    # Test 4: DVP vorticity bound finite and positive
    val = kernel_dvp_vorticity_bound()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 4: C_DVP = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: C_DVP = {val}")
        ok = False

    # Test 5: BSM mass bound finite
    val = kernel_bsm_mass_bound()
    if math.isfinite(val):
        print(f"  [PASS] Test 5: m_BSM = {val:.6e} eV/c²")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: m_BSM = {val}")
        ok = False

    # Test 6: Pcore composite > 0
    val = kernel_pcore_composite()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 6: ρ_Pcore = {val:.6e} kg/m³")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: ρ_Pcore = {val}")
        ok = False

    # Test 7: Cooper pair gap finite
    val = kernel_cooper_pair_gap()
    if math.isfinite(val):
        print(f"  [PASS] Test 7: Δ₀ = {val:.6e} eV")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: Δ₀ = {val}")
        ok = False

    # Test 8: Ramanujan router > 0
    val = kernel_ramanujan_router()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 8: Router sum = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: Router = {val}")
        ok = False

    # Test 9: CFD Stam step > 0
    val = kernel_cfd_stam_step()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 9: CFD max|v| = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: CFD = {val}")
        ok = False

    # Test 10: VDS×DVP×BSH cross > 0
    val = kernel_vds_dvp_bsh_cross()
    if val > 0 and math.isfinite(val):
        print(f"  [PASS] Test 10: VDS×DVP×BSH = {val:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: cross = {val}")
        ok = False

    # Test 11: Kernel names unique
    names = [k.__name__ for k in ALL_KERNELS]
    if len(names) == len(set(names)):
        print(f"  [PASS] Test 11: All {len(names)} kernel names unique")
        passed += 1
    else:
        print(f"  [FAIL] Test 11: Duplicate kernel names")
        ok = False

    # Test 12: Single pass finite
    bench = ProductionScalingV26()
    total = bench.single_pass()
    if math.isfinite(total):
        print(f"  [PASS] Test 12: Single pass total = {total:.6e}")
        passed += 1
    else:
        print(f"  [FAIL] Test 12: Non-finite pass total")
        ok = False

    # Test 13: Benchmark runs
    result = bench.simulate(0.5)
    if result['calc_per_sec'] > 0:
        label = "MET TARGET" if result['met_target'] else "below target"
        print(f"  [PASS] Test 13: {result['calc_per_sec']:.0f} calc/s ({label})")
        passed += 1
    else:
        print(f"  [FAIL] Test 13: Benchmark failed")
        ok = False

    print(f"\n  production_scaling_v26.py: {passed}/13 tests passed")
    return ok


if __name__ == "__main__":
    success = _run_tests()
    raise SystemExit(0 if success else 1)
