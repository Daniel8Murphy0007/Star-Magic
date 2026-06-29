# UQFF Star-Magic — Architecture Overview

A bird's-eye view of how the UQFF codebase is structured, from physics
primitives to public API to test infrastructure.

For the rules and constraints governing the codebase, see `CLAUDE.md`.
For the user-facing API, see `INPUT_DOMAINS.md`.

---

## High-level layout

```
Star-Magic/
├── uqff_pure_calculator.py        # The 2.66 MB single-file calculator (48,405 lines)
├── uqff_cli.py                    # CLI entry point: `uqff predict/list/status/...`
├── uqff_fidelity_tests.py         # 857-test fidelity gate (regression + structural)
├── uqff_exact_closures.cpp        # 368 C++ reference functions (no Python dep)
│
├── pyproject.toml                 # PEP 621 package metadata + build config
├── setup.py                       # Legacy setup for optional C++ pybind11 extension
├── MANIFEST.in                    # Files to ship in sdist
├── LICENSE                        # Dual AGPL-3.0 + commercial license notice
├── CITATION.cff                   # Canonical citation form (CFF 1.2.0)
│
├── CLAUDE.md                      # Project rules (read FIRST every session)
├── SESSION_LOG.md                 # Append-only history
├── PRODUCTION_ROADMAP.md          # 4-tier roadmap to commercial-grade
├── README.md                      # User-facing project intro
│
├── scripts/
│   ├── ci_smoke.py                # CI smoke contract (called by .github/workflows/ci.yml)
│   └── ci_strip_nulls.py          # Strip trailing null bytes from test file
│
├── docs/                          # Sphinx documentation source
│   ├── conf.py
│   ├── index.rst
│   └── api/, etc.
│
├── whitepapers/                   # 1,795 PAPER_XXXX.md files
│
└── .github/workflows/
    ├── ci.yml                     # Test matrix (3 OS × 4 Python) + coverage + build
    └── release.yml                # PyPI publish on tag push (Trusted Publishing OIDC)
```

---

## The calculator — internal structure

`uqff_pure_calculator.py` is intentionally a single file (CLAUDE.md Rule 6).
Inside, it follows this layered structure top-to-bottom:

### 1. Imports + canonical primitive constants (top of file, ~lines 1-200)

The 9 truly-independent primitives + 2 derivative + scaffolding constants
are defined as module-level uppercase names: `RHO_SCM`, `D_PHYS`, `D_CRIT`,
`SO_FIVE`, `A_FIVE`, `N_CH`, `TRZ`, `K_MEX`, `PHI_RESONANCE`, `BETA_I`,
etc. These are imported by anyone using the calculator.

### 2. Internal helper functions (~lines 200-9000)

Named with leading underscore: `_derive_constant`, `_paradox_proof`,
`_widom_larsen_*`, `_kozima_*`, `_meson_cascade`, etc. Not part of the
public API; subject to change.

### 3. The ~616 unique closure functions (~lines 9000-39000)

Each closure is a top-level `def _l96_uqff_axiom_<name>_closure() -> dict`
that takes no arguments and returns a dict with the closure's value(s),
residual, status tier, and primary source. Example:

```python
def _l96_uqff_axiom_hubble_tension_closure() -> Dict[str, Any]:
    uqff = TRZ**2 * D_CRIT * 67.4
    target = 67.4
    residual = abs(uqff - target) / abs(target) * 100.0
    return {
        'hubble_tension_target': target,
        'UQFF_formula_value': uqff,
        'residual_pct': residual,
        'status_tier': 'TIER0_DIRECT_LOCK',
        ...
    }
```

### 4. PARADOX_TO_CLOSURE dispatch table (~lines 39000-48000)

A dict mapping 794 lowercase paradox names to the closure functions above.
Multiple names can map to the same function (aliases): e.g.,
`olbers_paradox`, `alders_olbers`, `aldors_paradox` all map to
`_l96_uqff_axiom_olbers_paradox_closure`.

```python
PARADOX_TO_CLOSURE = {
    "olbers_paradox": _l96_uqff_axiom_olbers_paradox_closure,
    "alders_olbers": _l96_uqff_axiom_olbers_paradox_closure,  # alias
    ...
    "hubble_tension": _l96_uqff_axiom_hubble_tension_closure,
    ...
}
```

### 5. The 34 public `calculate_*` surfaces (bottom of file, ~lines 48000-48405)

Each takes a `dataset: dict` and returns `{'value': X}`. Most surfaces are
thin wrappers around either:

- a direct computation (e.g., `calculate_vacuum_ledger`)
- the dispatcher (`calculate_paradox` routes via `_paradox_proof`)
- a bucket aggregator (e.g., `calculate_cosmology` returns a list of 56 observables)

See `INPUT_DOMAINS.md` for the full per-surface API contract.

---

## The 9-sector Lagrangian

Physics-wise, UQFF organizes itself into 9 Lagrangian sectors:

```
L_UQFF = L_EH      (General Relativity in 26D)
       + L_YM      (Yang-Mills gauge; mass gap = 1.736 GeV)
       + L_Dirac   (Fermion / LENR; Kozima TNCF)
       + L_SCm     (SCm manifold; V(φ₀) = −ρ_SCm)
       + L_mag     (U_m magnetism; Heaviside amplifier)
       + L_buoy    (F_U_Bi_i buoyancy)
       + L_aether  (Two-component vacuum)
       + L_LENR    (Nuclear transmutation; COP parametric)
       + L_KK      (Kaluza-Klein 26D compactification)
```

See `GLOSSARY.md` for the canonical paper backing each sector.

---

## Buckets A-K — public surface organization

Closures are grouped into 11 "buckets" exposed as 9 public surfaces:

| Bucket | Surface | Closure count |
|---|---|---|
| A (Millennium) | `calculate_paradox` (routes to 8 millennium derivations) | 8 |
| B (paradoxes) | `calculate_paradox` (general dispatch) | ~80 |
| C (cosmology) | `calculate_cosmology` | 56 |
| D (particle physics) | `calculate_particle_physics` | 48 |
| E (GW events) | `calculate_gw_events` | 22 |
| F (AGN/jet) | `calculate_agn_jet` | 22 |
| G (astrophysics) | `calculate_astrophysics` | 36 |
| H (high-energy astro) | `calculate_high_energy_astro` | 10 |
| I (QGP) | `calculate_qgp` | 9 |
| J (Higgs precision) | `calculate_higgs_precision` | 13 |
| K (BSM constraints) | `calculate_bsm_constraints` | 17 |

Plus standalone surfaces: `calculate_universal_inertial_operator`,
`calculate_nuclear_magic`, `calculate_lenr`, `calculate_lenr_full`,
`calculate_f_u_zero`, `calculate_ua_layers`, `calculate_dpm_grinding`,
`calculate_caduceus`, `calculate_shell_orbital`,
`calculate_si_derivations`, `calculate_quantum_gravity`,
`calculate_black_hole`, `calculate_vds_dvp_bh26`,
`calculate_bsd_rank_cohomology`,
`calculate_negative_time_dual_existence`,
`calculate_status_report`, `calculate_whitepaper`,
`calculate_resonant_adpm`, `calculate_scm`, `calculate_f_u_bi`,
`calculate_f_u_bi_i`, `calculate_triadic_g`, `calculate_vacuum_ledger`,
`calculate_analytic_closures`.

---

## The fidelity gate

`uqff_fidelity_tests.py` is a single-file, no-pytest, no-dependencies
test harness. Runs in seconds. 857 tests organized into 60+ "blocks":

```
Block #1   — structural integrity (imports + 34 surfaces return correct shape)
Block #2-N — per-closure EXACT regression pins
Block #28  — PAPER_1521/1522 LANDMARK pins (D_BSFG, K_MEX EXACT)
Block #58  — legacy_freeform sweep (calls all 794 paradoxes, asserts non-None)
Block #59  — public surface sweep (calls all 34 surfaces with empty dataset)
Block #60  — bucket observable exercise (iterates observables lists)
```

Each test is a single-line `_exact()` or `check()` call. Pass/fail
recorded; non-zero exit if any failure.

Run locally: `python uqff_fidelity_tests.py`
Run via CLI: `uqff gate`

---

## The CLI (`uqff` command)

After `pip install uqff`, the `uqff` command exposes:

```
uqff predict <name>        → fetch a closure by lowercase name
uqff list [--filter X]     → list dispatch keys
uqff status [--json]       → production status summary
uqff surfaces [--json]     → list public calculate_* surfaces
uqff version               → version + headline metrics
uqff gate                  → run the fidelity gate
```

See `uqff_cli.py` for the implementation. Registered as a console script
in `pyproject.toml`:

```toml
[project.scripts]
uqff = "uqff_cli:main"
```

---

## CI/CD pipeline

### `ci.yml` runs on every push + PR:

1. **smoke** (ubuntu × py3.12) — diagnostics + null-strip + gate + smoke contract
2. **fidelity-gate** (3 OS × 4 Python = 12 jobs) — full gate, all combos
3. **coverage** — coverage.py + XML artifact + step summary
4. **build-package** — sdist + wheel + twine check + artifact upload
5. **smoke-test-install** (3 OS) — fresh-venv `pip install` + smoke contract

### `release.yml` runs on tag push (`v*`):

1. **build** — sdist + wheel + twine check
2. **publish-pypi** (if tag) — Trusted Publishing OIDC to PyPI
3. **github-release** — `gh release create` with curated notes + binaries

---

## Edit-tool safety

Per CLAUDE.md, the `uqff_pure_calculator.py` file is large enough (2.66 MB)
that simple text editors can truncate it mid-write. Documented mitigations:

- Use Python heredoc + `replace()` pattern instead of Edit/Write
- Always run with `PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/tmp`
- Verify with `wc -l` after every edit
- 6 truncation incidents recovered this project via Python splice from
  the protected backup files

---

## Performance characteristics

| Operation | Time | Notes |
|---|---|---|
| `import uqff_pure_calculator` (cold) | 2-5 s | One-time on first import |
| `import uqff_pure_calculator` (warm) | < 1 s | After .pyc cache |
| `calculate_paradox({"paradox": k})` | < 1 ms | Pure Python dispatch + arithmetic |
| `calculate_status_report({})` | ~50 ms | Iterates 263 schema-tagged closures |
| Full fidelity gate (857 tests) | ~3-5 s | On modern hardware |
| `pip install uqff` | ~10-15 s | 539 KB wheel |
| Resident memory after import | ~25 MB | Dominated by dispatch table |

---

## How to extend

To add a new closure (full workflow):

1. Read `CLAUDE.md` rules (especially 2, 3, 4, 5, 8, 9)
2. Author whitepaper: `whitepapers/PAPER_<NNNN>_<description>.md`
3. Add closure function to `uqff_pure_calculator.py`:
   ```python
   def _l96_uqff_axiom_<your_name>_closure() -> Dict[str, Any]:
       uqff = <your_formula_using_primitives>
       target = <measured_value>
       residual = abs(uqff - target) / abs(target) * 100.0
       return {
           '<your_name>_target': target,
           'UQFF_formula_value': uqff,
           'residual_pct': residual,
           'status_tier': '<EXACT|sub_0.01|...>',
           'description': '<one line>',
           'primary_source': 'PAPER_<NNNN>_<slug>'
       }
   ```
4. Register dispatch entry (lowercase only):
   ```python
   PARADOX_TO_CLOSURE = {
       ...,
       "<your_name>": _l96_uqff_axiom_<your_name>_closure,
   }
   ```
5. Add regression pin to `uqff_fidelity_tests.py`:
   ```python
   _exact("<Your Label> PAPER_<NNNN>", <python_expression>, <expected_value>, tol=1e-6)
   ```
6. Run gate: `python uqff_fidelity_tests.py` — must exit 0
7. Append SESSION_LOG.md entry
8. Submit PR

See `docs/contributing.rst` for the full guide and PR template.

---

## What's deliberately NOT in the architecture

Per CLAUDE.md rules, the calculator deliberately omits:

- **Docstrings** (Rule 3 — purity requirement)
- **Comments** explaining what code does (Rule 3)
- **datetime / file-writes / classes / `__main__`** (Rule 6)
- **Standard-Model named constants** (Rule 4 — no SM contamination)
- **Provenance metadata strings** inside closures (Rule 5 — `{'value': X}` only)

These restrictions ARE deliberate. External tooling (this ARCHITECTURE.md,
the whitepapers, the Sphinx docs) provides the explanatory layer that the
calculator source intentionally omits.
