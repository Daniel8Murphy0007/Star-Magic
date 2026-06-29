# UQFF Performance Profile (Tier-3 G10)

**Iterations per surface:** 50
**Python:** 3.10.12
**Platform:** linux
**Calculator:** uqff_pure_calculator.py (2.54 MB, 794 dispatch keys)

## Headline numbers

| Metric | Value |
|---|---|
| Cold import (median of 3 subprocesses) | **0.54 s** |
| Warm import (already-loaded module) | **1453.61 ms** |
| Module memory footprint (allocated) | **6914 KB** |
| Median per-dispatch closure call | **3.41 µs** |
| p95 per-dispatch closure call | **7.02 µs** |
| Slowest single dispatch call (sampled) | **10.01 µs** (`coronal_heating`) |

## Per-surface latency (slowest first)

| Surface | Median | p95 | Min | Max |
|---|---:|---:|---:|---:|
| `calculate_status_report` | 9.70 ms | 16.60 ms | 4.84 ms | 19.59 ms |
| `calculate_paradox` | 6.58 ms | 8.86 ms | 3.70 ms | 10.66 ms |
| `calculate_vds_dvp_bh26` | 434.93 µs | 581.49 µs | 250.21 µs | 667.25 µs |
| `calculate_vacuum_ledger` | 160.61 µs | 301.55 µs | 93.00 µs | 464.04 µs |
| `calculate_scm` | 158.38 µs | 534.92 µs | 82.68 µs | 801.05 µs |
| `calculate_resonant_adpm` | 151.30 µs | 301.79 µs | 86.73 µs | 881.52 µs |
| `calculate_f_u_bi_i` | 134.74 µs | 219.36 µs | 128.37 µs | 233.85 µs |
| `calculate_cosmology` | 104.55 µs | 147.02 µs | 66.74 µs | 160.38 µs |
| `calculate_triadic_g` | 102.48 µs | 196.91 µs | 91.39 µs | 583.73 µs |
| `calculate_whitepaper` | 93.54 µs | 179.13 µs | 90.63 µs | 197.60 µs |
| `calculate_f_u_bi` | 93.12 µs | 299.58 µs | 80.15 µs | 570.59 µs |
| `calculate_analytic_closures` | 85.35 µs | 138.12 µs | 78.88 µs | 889.94 µs |
| `calculate_particle_physics` | 48.65 µs | 91.40 µs | 39.87 µs | 223.60 µs |
| `calculate_astrophysics` | 34.46 µs | 65.44 µs | 33.72 µs | 128.32 µs |
| `calculate_gw_events` | 24.69 µs | 53.17 µs | 14.07 µs | 98.29 µs |
| `calculate_caduceus` | 14.87 µs | 29.01 µs | 13.45 µs | 54.50 µs |
| `calculate_agn_jet` | 14.82 µs | 71.11 µs | 12.90 µs | 124.75 µs |
| `calculate_lenr_full` | 13.55 µs | 15.63 µs | 11.52 µs | 64.96 µs |
| `calculate_higgs_precision` | 10.16 µs | 11.76 µs | 8.11 µs | 37.64 µs |
| `calculate_bsm_constraints` | 8.63 µs | 11.27 µs | 8.46 µs | 33.60 µs |
| `calculate_high_energy_astro` | 8.29 µs | 9.64 µs | 5.16 µs | 63.31 µs |
| `calculate_qgp` | 7.42 µs | 7.96 µs | 6.85 µs | 32.18 µs |
| `calculate_f_u_zero` | 6.20 µs | 7.62 µs | 5.11 µs | 54.98 µs |
| `calculate_lenr` | 4.95 µs | 7.17 µs | 3.30 µs | 40.97 µs |
| `calculate_negative_time_dual_existence` | 4.90 µs | 6.66 µs | 2.49 µs | 31.90 µs |
| `calculate_nuclear_magic` | 4.65 µs | 5.62 µs | 2.99 µs | 23.90 µs |
| `calculate_shell_orbital` | 2.65 µs | 3.15 µs | 2.22 µs | 18.08 µs |
| `calculate_ua_layers` | 1.96 µs | 2.01 µs | 1.90 µs | 12.71 µs |
| `calculate_si_derivations` | 1.73 µs | 2.60 µs | 1.68 µs | 19.39 µs |
| `calculate_quantum_gravity` | 1.65 µs | 2.03 µs | 1.45 µs | 48.20 µs |
| `calculate_bsd_rank_cohomology` | 1.23 µs | 1.43 µs | 1.17 µs | 37.07 µs |
| `calculate_black_hole` | 1.14 µs | 1.36 µs | 1.11 µs | 15.02 µs |
| `calculate_universal_inertial_operator` | 920.00 ns | 980.04 ns | 879.98 ns | 3.56 µs |
| `calculate_dpm_grinding` | 419.97 ns | 480.04 ns | 379.98 ns | 5.03 µs |

## 10 slowest dispatched closures (50 sampled at random)

| Paradox key | Median per call |
|---|---:|
| `coronal_heating` | 10.01 µs |
| `final_parsec` | 8.32 µs |
| `info_paradox_full_qg` | 7.02 µs |
| `dm_particle` | 6.36 µs |
| `faber_jackson_relation` | 5.97 µs |
| `infinity_paradox` | 5.42 µs |
| `faint_young_sun_paradox` | 5.39 µs |
| `smooth_poincare_4d` | 4.98 µs |
| `collatz` | 4.85 µs |
| `cmb_cold_spot` | 4.75 µs |

## Interpretation

- **Cold import** is the time from `import uqff_pure_calculator` (no .pyc cache) to module ready. The single 48k-line / 2.66 MB file is parsed once. Subsequent imports in the same process are sub-millisecond.
- **Per-dispatch closure call** is the round-trip through `calculate_paradox({'paradox': key})` -> `_paradox_proof` -> the named closure function -> JSON-serializable dict return.
- **Bucket surfaces** (cosmology, particle_physics, etc.) are the slowest single-call surfaces because they iterate their full observable list (60+ observables for cosmology) on each call.
- **Memory footprint** is the module's resident allocation. The 794-key dispatch table dominates; closure function bodies are amortized at module load.

## Tier-3 implications

- Cold import < 5s: acceptable for one-off CLI invocations; would be the bottleneck for a server with many short-lived workers. The REST API (uqff_api.py) only imports once at startup, so this is a non-issue for /predict endpoints.
- Per-dispatch median around 50-500 µs: excellent for interactive use. A single client can issue >2,000 closure lookups per second.
- The slowest bucket surfaces (a few ms each) are still well within interactive latency. Tier-3 H1 (modular refactor) could amortize the cold-import cost by lazy-loading surfaces on first access.
