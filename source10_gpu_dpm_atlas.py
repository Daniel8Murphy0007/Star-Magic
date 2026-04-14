#!/usr/bin/env python3
"""
source10_gpu_dpm_atlas.py — GPU-Vectorized DPM S₂₆⁽³⁾ Spectral Atlas

Session 224 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
GPU-batched Dipole Moment (DPM) spectral atlas engine.
Computes S₂₆⁽³⁾ spectral profiles across 26 quantum states using torch
tensors with automatic CUDA fallback to CPU.

Core operation:
    DPM(ω, i) = [SSq]^i / i^26 · R_n(i,3) · exp(-(ω - ω_SCm)² / 2σ²)
    Atlas(ω)  = Σ_{i=1}^{26} DPM(ω, i)     [batched matmul across states]

ALMA Cycle 12 targets: Orion M42, Lagoon M8, Eagle M16, Carina, Trifid M20,
                        Omega M17, Rosette NGC 2237, Flame NGC 2024

Architecture: Pure calculator. No hardcoded systems. Tier 2 compute.
              Falls back to NumPy/pure-Python if torch unavailable.
────────────────────────────────────────────────────────────────────────────────
"""

import math
import time
from typing import Any, Dict, List, Optional, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
M_SUN     = 1.989e30

OMEGA_SCM = 2 * PI * 1.25e12       # rad/s  SCm phonon resonance
SSQ       = 0.57                    # string sector coupling
BETA_I    = 0.603                   # buoyancy coupling
GAMMA_0   = 2 * PI * 0.1e12        # rad/s  linewidth center
SIGMA_G   = 0.08 * 2 * PI * 1e12   # rad/s  linewidth sigma
N_LAYERS  = 26                      # 26D quantum states

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def _ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor R_n^{(26,k)}."""
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


S26_3RD = sum((SSQ ** n) / (n ** 26) * _ramanujan_Rn(n, 3) for n in range(1, 27))


# Pre-compute per-layer coefficients c_i = [SSq]^i / i^26 · R_n(i,3)
_LAYER_COEFFS = [
    (SSQ ** i) / (i ** 26) * _ramanujan_Rn(i, 3) for i in range(1, N_LAYERS + 1)
]


# ── §1  BACKEND ABSTRACTION ───────────────────────────────────────────────

_BACKEND = None  # 'torch', 'numpy', or 'pure'
_torch = None
_np = None


def _detect_backend():
    """Detect best available compute backend."""
    global _BACKEND, _torch, _np
    if _BACKEND is not None:
        return _BACKEND
    try:
        import torch as _t
        _torch = _t
        _BACKEND = 'torch'
    except ImportError:
        pass
    if _BACKEND is None:
        try:
            import numpy as _n
            _np = _n
            _BACKEND = 'numpy'
        except ImportError:
            _BACKEND = 'pure'
    return _BACKEND


def _get_device() -> str:
    """Return best torch device string."""
    if _BACKEND == 'torch' and _torch is not None:
        return 'cuda' if _torch.cuda.is_available() else 'cpu'
    return 'cpu'


# ── §2  GPU/CPU DPM SPECTRAL CALCULATOR ───────────────────────────────────

class DPMSpectralAtlas:
    """
    GPU-vectorized Dipole Moment spectral atlas across 26 quantum states.

    For each frequency bin ω_j in the input grid, computes:

        DPM(ω_j) = Σ_{i=1}^{26} c_i · Φ_i(ω_j)

    where:
        c_i      = [SSq]^i / i^26 · R_n(i,3)         [layer coefficient]
        Φ_i(ω_j) = exp(-(ω_j - ω_SCm)² / 2σ_i²)     [Gaussian profile]
        σ_i      = σ_G · (1 + 0.02·(i-1))             [slight broadening per layer]

    The matmul formulation:
        C = [c_1, ..., c_26]                   shape (26,)
        G = [[Φ_1(ω_1), ..., Φ_1(ω_N)],       shape (26, N_freq)
             [Φ_26(ω_1), ..., Φ_26(ω_N)]]
        Atlas = C @ G                          shape (N_freq,)

    Variables (from dataset):
        freq_min_THz:   minimum frequency (THz, default 0.5)
        freq_max_THz:   maximum frequency (THz, default 2.0)
        n_freq_bins:    number of frequency bins (default 1024)
        omega_scm:      center frequency (rad/s, default OMEGA_SCM)
        sigma_g:        base linewidth (rad/s, default SIGMA_G)
        broadening:     per-layer broadening factor (default 0.02)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        backend = _detect_backend()
        device = _get_device()

        freq_min = dataset.get('freq_min_THz', 0.5)
        freq_max = dataset.get('freq_max_THz', 2.0)
        n_bins   = dataset.get('n_freq_bins', 1024)
        omega_c  = dataset.get('omega_scm', OMEGA_SCM)
        sig_base = dataset.get('sigma_g', SIGMA_G)
        broadening = dataset.get('broadening', 0.02)

        t_start = time.perf_counter()

        if backend == 'torch':
            result = self._compute_torch(
                freq_min, freq_max, n_bins, omega_c, sig_base, broadening, device
            )
        elif backend == 'numpy':
            result = self._compute_numpy(
                freq_min, freq_max, n_bins, omega_c, sig_base, broadening
            )
        else:
            result = self._compute_pure(
                freq_min, freq_max, n_bins, omega_c, sig_base, broadening
            )

        elapsed_ms = (time.perf_counter() - t_start) * 1000

        omega_grid = result['omega_grid']
        atlas = result['atlas']
        peak_idx = max(range(len(atlas)), key=lambda i: atlas[i])

        return {
            'primary_equations': [
                f"DPM(ω) = Σ_{{i=1}}^{{26}} ([SSq]^i / i^26) · R_n(i,3) · Φ_i(ω)",
                f"Atlas(ω) = C @ G   [matmul: (26,) × (26, {n_bins}) → ({n_bins},)]",
                f"σ_i = σ_G · (1 + {broadening}·(i-1))  [per-layer broadening]",
            ],
            'atlas_values': atlas,
            'omega_grid_rad_s': omega_grid,
            'peak_omega_rad_s': omega_grid[peak_idx],
            'peak_freq_THz': omega_grid[peak_idx] / (2 * PI * 1e12),
            'peak_amplitude': atlas[peak_idx],
            'S26_3RD': S26_3RD,
            'n_layers': N_LAYERS,
            'n_freq_bins': n_bins,
            'backend': backend,
            'device': device,
            'elapsed_ms': elapsed_ms,
            'layer_coefficients': list(_LAYER_COEFFS),
        }

    # ── Torch backend ─────────────────────────────────────────────────────

    def _compute_torch(self, fmin, fmax, n_bins, omega_c, sig_base, broad, device):
        torch = _torch
        omega = torch.linspace(
            fmin * 2 * PI * 1e12, fmax * 2 * PI * 1e12, n_bins,
            device=device, dtype=torch.float64
        )
        coeffs = torch.tensor(_LAYER_COEFFS, device=device, dtype=torch.float64)

        # Build Gaussian matrix (26, n_bins)
        gauss = torch.zeros(N_LAYERS, n_bins, device=device, dtype=torch.float64)
        for i in range(N_LAYERS):
            sigma_i = sig_base * (1.0 + broad * i)
            gauss[i] = torch.exp(-((omega - omega_c) ** 2) / (2 * sigma_i ** 2))

        # Batched matmul: (26,) @ (26, n_bins) → (n_bins,)
        atlas = torch.matmul(coeffs, gauss)

        return {
            'omega_grid': omega.cpu().tolist(),
            'atlas': atlas.cpu().tolist(),
        }

    # ── NumPy backend ─────────────────────────────────────────────────────

    def _compute_numpy(self, fmin, fmax, n_bins, omega_c, sig_base, broad):
        np = _np
        omega = np.linspace(fmin * 2 * PI * 1e12, fmax * 2 * PI * 1e12, n_bins)
        coeffs = np.array(_LAYER_COEFFS)

        gauss = np.zeros((N_LAYERS, n_bins))
        for i in range(N_LAYERS):
            sigma_i = sig_base * (1.0 + broad * i)
            gauss[i] = np.exp(-((omega - omega_c) ** 2) / (2 * sigma_i ** 2))

        atlas = coeffs @ gauss
        return {'omega_grid': omega.tolist(), 'atlas': atlas.tolist()}

    # ── Pure Python backend ───────────────────────────────────────────────

    def _compute_pure(self, fmin, fmax, n_bins, omega_c, sig_base, broad):
        omega = [
            fmin * 2 * PI * 1e12 + j * (fmax - fmin) * 2 * PI * 1e12 / (n_bins - 1)
            for j in range(n_bins)
        ]
        atlas = [0.0] * n_bins
        for j in range(n_bins):
            total = 0.0
            for i in range(N_LAYERS):
                sigma_i = sig_base * (1.0 + broad * i)
                phi = math.exp(-((omega[j] - omega_c) ** 2) / (2 * sigma_i ** 2))
                total += _LAYER_COEFFS[i] * phi
            atlas[j] = total
        return {'omega_grid': omega, 'atlas': atlas}


# ── §3  MULTI-SYSTEM BATCH ATLAS ──────────────────────────────────────────

class DPMBatchAtlas:
    """
    Batch DPM atlas across multiple astrophysical systems.

    Accepts a list of system datasets, each with their own physical
    parameters (mass, distance, SFR, etc.) plus optional frequency
    overrides. Computes DPM atlas for each and aggregates.

    Variables (from dataset):
        systems:  list of dicts, each containing:
            name:          system name (e.g. "Orion M42")
            freq_min_THz:  minimum frequency (THz)
            freq_max_THz:  maximum frequency (THz)
            n_freq_bins:   number of bins
            M:             mass (kg, for gravity coupling)
            r:             distance (m, for flux scaling)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        systems = dataset.get('systems', [])
        if not systems:
            return {'primary_equations': [], 'error': 'No systems provided'}

        engine = DPMSpectralAtlas()
        results = []
        t_start = time.perf_counter()

        for sys in systems:
            res = engine.compute(sys)
            res['system_name'] = sys.get('name', 'unnamed')
            if 'M' in sys and 'r' in sys:
                # Gravity coupling weight
                g_surf = G * sys['M'] / sys['r'] ** 2
                res['gravity_coupling_m_s2'] = g_surf
                # Scale atlas by gravity factor normalized to solar
                g_solar = G * M_SUN / (6.96e8) ** 2
                scale = g_surf / g_solar if g_solar > 0 else 1.0
                res['scaled_peak'] = res['peak_amplitude'] * scale
            results.append(res)

        elapsed_ms = (time.perf_counter() - t_start) * 1000

        return {
            'primary_equations': [
                "Batch DPM: Atlas_k(ω) = (g_k/g_☉) · Σ_{i=1}^{26} c_i · Φ_i(ω)",
                f"Systems processed: {len(results)}",
            ],
            'system_results': results,
            'n_systems': len(results),
            'total_elapsed_ms': elapsed_ms,
        }


# ── §4  SPECTRAL LINE IDENTIFIER ──────────────────────────────────────────

class DPMLineIdentifier:
    """
    Identify DPM spectral lines (peaks) in the atlas above a threshold.

    Scans the atlas output for local maxima exceeding a given fraction
    of the global peak. Returns line positions, widths, and amplitudes.

    Variables (from dataset):
        atlas_values:     list of atlas amplitudes
        omega_grid_rad_s: list of frequency values (rad/s)
        threshold_frac:   minimum fraction of peak to count as line (default 0.01)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        atlas = dataset.get('atlas_values', [])
        omega = dataset.get('omega_grid_rad_s', [])
        threshold = dataset.get('threshold_frac', 0.01)

        if len(atlas) < 3 or len(omega) < 3:
            return {'primary_equations': [], 'lines': [], 'error': 'Insufficient data'}

        peak_val = max(atlas)
        cutoff = threshold * peak_val

        lines = []
        for j in range(1, len(atlas) - 1):
            if atlas[j] >= atlas[j - 1] and atlas[j] > atlas[j + 1] and atlas[j] > cutoff:
                # Estimate FWHM by half-max crossing
                half_max = atlas[j] / 2.0
                left = j
                while left > 0 and atlas[left] > half_max:
                    left -= 1
                right = j
                while right < len(atlas) - 1 and atlas[right] > half_max:
                    right += 1
                fwhm_rad = omega[right] - omega[left] if right > left else 0.0
                lines.append({
                    'omega_rad_s': omega[j],
                    'freq_THz': omega[j] / (2 * PI * 1e12),
                    'amplitude': atlas[j],
                    'fwhm_rad_s': fwhm_rad,
                    'fwhm_GHz': fwhm_rad / (2 * PI * 1e9),
                    'relative_strength': atlas[j] / peak_val if peak_val > 0 else 0.0,
                })

        return {
            'primary_equations': [
                "Line detection: local maxima above threshold fraction of peak",
                f"Threshold: {threshold:.3f} × peak ({peak_val:.4e})",
                f"Lines found: {len(lines)}",
            ],
            'lines': lines,
            'n_lines': len(lines),
            'peak_amplitude': peak_val,
        }


# ── §5  ALMA TARGET PROFILE GENERATOR ─────────────────────────────────────

class ALMATargetProfileGenerator:
    """
    Generate predicted UQFF DPM spectral profiles for ALMA Cycle 12 targets.

    Each target receives a DPM atlas computation at standard ALMA Band 6
    (211-275 GHz / 0.211-0.275 THz) and Band 7 (275-373 GHz) resolution.

    Variables (from dataset):
        targets:     list of target dicts with {name, M, r, z}
        band:        'B6' (211-275 GHz), 'B7' (275-373 GHz), or 'full' (default)
        n_channels:  spectral channels (default 2048)
    """

    ALMA_BANDS = {
        'B6': (0.211, 0.275),   # THz
        'B7': (0.275, 0.373),   # THz
        'B3': (0.084, 0.116),   # THz
        'B4': (0.125, 0.163),   # THz
        'full': (0.084, 0.950), # THz  (full ALMA range)
    }

    def compute(self, dataset: dict) -> Dict[str, Any]:
        targets = dataset.get('targets', [])
        band = dataset.get('band', 'full')
        n_ch = dataset.get('n_channels', 2048)

        if band not in self.ALMA_BANDS:
            return {'primary_equations': [], 'error': f'Unknown band: {band}'}

        fmin, fmax = self.ALMA_BANDS[band]

        batch_engine = DPMBatchAtlas()
        systems = []
        for t in targets:
            systems.append({
                'name': t.get('name', 'unnamed'),
                'freq_min_THz': fmin,
                'freq_max_THz': fmax,
                'n_freq_bins': n_ch,
                'M': t.get('M', M_SUN),
                'r': t.get('r', 6.96e8),
            })

        result = batch_engine.compute({'systems': systems})

        return {
            'primary_equations': [
                f"ALMA Band {band}: {fmin:.3f}–{fmax:.3f} THz",
                f"Channels: {n_ch}",
                f"DPM atlas per target with gravity-scaled amplitude",
            ],
            'band': band,
            'freq_range_THz': (fmin, fmax),
            'n_channels': n_ch,
            'target_profiles': result.get('system_results', []),
            'n_targets': len(targets),
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  SELF-TESTS
# ══════════════════════════════════════════════════════════════════════════════

def _run_tests():
    print("=" * 72)
    print("source10_gpu_dpm_atlas.py — Self-Tests")
    print("=" * 72)

    passed = 0
    failed = 0

    # Test 1: Backend detection
    backend = _detect_backend()
    assert backend in ('torch', 'numpy', 'pure'), f"Bad backend: {backend}"
    print(f"  [PASS] Test 1: Backend detected = {backend}")
    passed += 1

    # Test 2: Single atlas computation
    engine = DPMSpectralAtlas()
    res = engine.compute({'n_freq_bins': 512})
    assert 'atlas_values' in res
    assert len(res['atlas_values']) == 512
    assert res['peak_amplitude'] > 0
    print(f"  [PASS] Test 2: Single atlas, peak = {res['peak_amplitude']:.6e}")
    passed += 1

    # Test 3: Peak near 1.25 THz
    peak_THz = res['peak_freq_THz']
    assert 1.0 < peak_THz < 1.5, f"Peak at {peak_THz} THz, expected ~1.25"
    print(f"  [PASS] Test 3: Peak frequency = {peak_THz:.4f} THz")
    passed += 1

    # Test 4: S26_3RD consistency
    assert abs(res['S26_3RD'] - S26_3RD) < 1e-20
    print(f"  [PASS] Test 4: S26_3RD = {S26_3RD:.6e}")
    passed += 1

    # Test 5: Batch atlas
    batch = DPMBatchAtlas()
    bres = batch.compute({'systems': [
        {'name': 'Orion M42', 'M': 2000 * M_SUN, 'r': 1.3e16},
        {'name': 'Lagoon M8', 'M': 5000 * M_SUN, 'r': 4.0e19},
    ]})
    assert bres['n_systems'] == 2
    assert len(bres['system_results']) == 2
    print(f"  [PASS] Test 5: Batch atlas, {bres['n_systems']} systems")
    passed += 1

    # Test 6: Gravity coupling computed
    for sr in bres['system_results']:
        assert 'gravity_coupling_m_s2' in sr
        assert sr['gravity_coupling_m_s2'] > 0
    print(f"  [PASS] Test 6: Gravity coupling present for all systems")
    passed += 1

    # Test 7: Line identifier
    lid = DPMLineIdentifier()
    lres = lid.compute({
        'atlas_values': res['atlas_values'],
        'omega_grid_rad_s': res['omega_grid_rad_s'],
        'threshold_frac': 0.01,
    })
    assert lres['n_lines'] >= 1
    print(f"  [PASS] Test 7: Line identifier found {lres['n_lines']} lines")
    passed += 1

    # Test 8: ALMA target profiles
    alma = ALMATargetProfileGenerator()
    ares = alma.compute({
        'targets': [
            {'name': 'Orion M42', 'M': 2000 * M_SUN, 'r': 1.3e16},
            {'name': 'Eagle M16', 'M': 8000 * M_SUN, 'r': 5.5e19},
        ],
        'band': 'full',
        'n_channels': 256,
    })
    assert ares['n_targets'] == 2
    assert len(ares['target_profiles']) == 2
    print(f"  [PASS] Test 8: ALMA profiles for {ares['n_targets']} targets")
    passed += 1

    # Test 9: Band 6 specific range
    ares6 = alma.compute({
        'targets': [{'name': 'test', 'M': M_SUN, 'r': 6.96e8}],
        'band': 'B6',
        'n_channels': 128,
    })
    assert ares6['freq_range_THz'] == (0.211, 0.275)
    print(f"  [PASS] Test 9: ALMA Band 6 freq range correct")
    passed += 1

    # Test 10: Primary equations present
    assert len(res['primary_equations']) >= 1
    assert len(bres['primary_equations']) >= 1
    assert len(lres['primary_equations']) >= 1
    print(f"  [PASS] Test 10: Primary equations in all outputs")
    passed += 1

    # Test 11: Performance (should be < 5s even on CPU pure Python)
    assert res['elapsed_ms'] < 5000, f"Too slow: {res['elapsed_ms']:.1f} ms"
    print(f"  [PASS] Test 11: Atlas computed in {res['elapsed_ms']:.2f} ms")
    passed += 1

    print("-" * 72)
    print(f"Results: {passed}/{passed + failed} passed, {failed} failed")
    return passed, failed


if __name__ == '__main__':
    _run_tests()
