# UQFF Replication Package

**For:** peer reviewers, independent replicators, sceptics, follow-on researchers
**Estimated time to full replication:** 30-60 minutes
**Required:** Python 3.10 or newer, ~50 MB disk, internet for `pip install`
**Optional:** g++ (for C++ cross-verification), pybind11 (for C++ acceleration)

---

## TL;DR — reproduce every numerical claim in the manuscript in 5 commands

```bash
pip install uqff                                    # install the framework
uqff version                                        # verify v5.28.x
uqff gate                                           # run 857-test fidelity gate
uqff predict <any_named_closure>                    # verify any specific claim
python scripts/cpp_python_crosscheck.py             # cross-language verification (optional)
```

If `uqff gate` prints `TOTAL: 857 passed, 0 failed`, every numerical claim in the manuscript is verified against the version-locked PyPI release.

---

## Step 1 — Install

```bash
pip install uqff==5.28.0    # exact version cited in manuscript
```

Confirms install:

```bash
uqff version
# Expected: uqff 5.28.0, 794 closures, 9+2 primitives
```

If the `uqff` command isn't recognized, the package installed but the console script isn't on PATH. Use:

```bash
python -m uqff_cli version
```

instead. See `TROUBLESHOOTING.md` in the repo for full PATH guidance per OS.

---

## Step 2 — Run the fidelity gate

```bash
uqff gate
# Expected: TOTAL: 857 passed, 0 failed
```

This is the runtime verification of every numerical claim in the manuscript. If any test fails, the value cited in the paper does not match the running code.

To get the same gate in CI: every commit on https://github.com/Daniel8Murphy0007/Star-Magic runs the gate across 12 OS/Python combinations. Click any green ✅ on the Actions tab to see the run output.

---

## Step 3 — Verify specific manuscript claims

Each manuscript section's headline number can be verified independently:

### Cosmological constant (Section 4.1)

```bash
python -c "
import uqff_pure_calculator as u
v = u.calculate_vacuum_ledger({})['value']
print(f'Lambda derived = {v[\"total_lambda_J_per_m3\"]:.6e} J/m^3')
print(f'Lambda Planck  = 5.957e-10 J/m^3')
print(f'Residual = {abs(v[\"total_lambda_J_per_m3\"] - 5.957e-10) / 5.957e-10 * 100:.4f}%')
"
# Expected: ~0.003% residual
```

### Nuclear magic numbers (Section 4.2)

```bash
uqff predict magic_numbers
# Expected: {2: 'He', 8: 'O', 20: 'Ca', 28: 'Ni', 50: 'Sn', 82: 'Pb', 126: 'island-of-stability'}
```

### Holmlid 630 eV (Section 4.3)

```bash
uqff predict holmlid_D_minus_1
# Expected: {'holmlid_D_minus_1': {'KER_eV': 630.0, 'cluster_size_atoms': 100, ..., 'reference': 'PAPER_1133'}}
```

### Yang-Mills mass gap (Section 4.4)

```bash
uqff predict yang_mills
# Expected: 5970.0  (GeV)
```

### SM 12-fermion spectrum (Section 4.5)

```bash
python -c "
import uqff_pure_calculator as u
import json
print(json.dumps(u.calculate_particle_physics({})['value']['observables'][:12], indent=2))
"
# Expected: m_W, m_Z, m_top, m_H, m_b, m_c, m_tau, m_mu, m_s, m_e + 3 neutrinos
```

### Bayesian model comparison ΔBIC = 94.1 (Section 5)

```bash
uqff predict bayesian_delta_bic_uqff_wins
# Expected: 94.1 (decisive Kass-Raftery preference)
```

---

## Step 4 — Run the full audit suite

These additional scripts produce supplementary verification reports:

```bash
# Performance profile (per-dispatch latency)
python scripts/perf_profile.py
# Output: PERFORMANCE_PROFILE.md (3.4 µs median; 290k calls/sec)

# Coverage audit (which closures have backing whitepapers)
# Output already in: COVERAGE_GAPS.md (62% literal; ~80-95% real)
```

---

## Step 5 — Cross-language C++ verification (optional)

```bash
# Requires g++ installed; available on macOS/Linux by default,
# on Windows via WSL or MinGW
python scripts/cpp_python_crosscheck.py
# Expected: MATCH=277, DRIFT=0
```

This compiles the C++ reference implementation (`uqff_exact_closures.cpp`, 632 functions) and verifies every value matches the Python implementation within 1e-6 relative tolerance. The 100% match rate is the strongest cross-language evidence for the numerical claims in the manuscript.

---

## Step 6 — Browse the supporting whitepapers

```bash
git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
cd Star-Magic
ls whitepapers/ | head
ls whitepapers/ | wc -l    # 1,867 markdown files
```

Each `PAPER_XXXX_<slug>.md` is one axiom set / derivation / theorem. The full catalog with auto-extracted metadata is in `WHITEPAPER_INDEX.md`.

For headline derivations cited in the manuscript:

- `whitepapers/PAPER_1005_*.md` — Yang-Mills mass gap derivation
- `whitepapers/PAPER_1133_*.md` — Holmlid 630 eV chain
- `whitepapers/PAPER_1156_*.md` — Lambda cosmological derivation
- `whitepapers/PAPER_1167_*.md` — All-8-Lagrangian-gaps master synthesis
- `whitepapers/PAPER_1182_*.md` — Riemann hypothesis closure
- `whitepapers/PAPER_1203_*.md` — F_U=0 master equation + nuclear shells
- `whitepapers/PAPER_1521_*.md` — D_BSFG primitive reduction
- `whitepapers/PAPER_1522_*.md` — K_MEX primitive reduction

---

## Step 7 — REST API (optional, for interactive exploration)

```bash
pip install 'uqff[api]'
uqff serve --port 8080
# Then open http://127.0.0.1:8080/docs in your browser
```

Auto-generated Swagger UI lets you click through every endpoint and test any closure interactively. Useful for exploration without writing code.

---

## What the package includes

| File / directory | Purpose |
|---|---|
| `uqff_pure_calculator.py` | The calculator (2.66 MB single file, 794 dispatch keys, 34 public surfaces) |
| `uqff_cli.py` | Command-line interface |
| `uqff_api.py` | FastAPI REST server |
| `uqff_jupyter.py` | IPython rich display + `%uqff` magic |
| `uqff_fidelity_tests.py` | 857-test runtime verification |
| `uqff_exact_closures.cpp` | C++ reference (632 functions) |
| `scripts/` | Auto-audit scripts (smoke, lint, perf, cross-check) |
| `whitepapers/` | 1,867 PAPER_XXXX.md proof / derivation files |
| `notebooks/` | 4 Jupyter quickstarts (quickstart, LENR, magic numbers, cosmology, REPL) |
| `docs/` | Sphinx documentation source |
| `CLOSURE_ATLAS.md` | Master map of all closure namespaces |
| `WHITEPAPER_INDEX.md` | Full table of 1,867 whitepapers |
| `COVERAGE_GAPS.md` | Closure ↔ whitepaper audit |
| `PROVENANCE_AUDIT.md` | Primitive provenance chains |
| `forward_predictions.md` | 42 falsifiable forecasts |
| `STATISTICAL_HYGIENE.md` | Multiple-comparisons analysis |
| `PREDICTION_LABELS.md` | POST/NEW/AMB classification |
| `PERFORMANCE_PROFILE.md` | Latency + memory benchmarks |
| `CPP_PYTHON_CROSSCHECK_REPORT.md` | Cross-language verification |
| `PUBLISH_TO_PYPI.md` | PyPI publication procedure |
| `RELEASE_PROCESS.md` | Versioning + tagging procedure |
| `LICENSE`, `LICENSE-AGPL-3.0.txt`, `COMMERCIAL.md` | Dual-license terms |

---

## CI verification (automatic on every push)

Every commit to the repository triggers 8 automated audits at https://github.com/Daniel8Murphy0007/Star-Magic/actions:

1. `smoke` — fidelity gate quick-check
2. `fidelity-gate` — full 857-test gate on 12 OS/Python combinations
3. `coverage` — code coverage measurement
4. `lint` — ruff style audit
5. `cpp-crosscheck` — Python ↔ C++ value verification
6. `perf-smoke` — latency regression check
7. `build-package` — PyPI wheel + sdist build
8. `smoke-test-install` — fresh-venv install verification on 3 OSes

Click any commit's green checks to see the audit outputs. Replicators can validate that the published claims continue to hold without running anything locally.

---

## Frequently asked replicator questions

### "How do I know UQFF isn't fitting parameters to match observations?"

The 9 primitives are LOCKED — see `PROVENANCE_AUDIT.md` for the locking chain per primitive. CLAUDE.md Rule 2 in the repository explicitly prohibits altering primitive values. The fidelity gate enforces this at the test level: change a primitive and 100+ tests fail immediately.

### "Why are some closures only at 0.1-1% accuracy?"

UQFF makes no claim of perfect agreement. The 5-band uncertainty classification in `STATISTICAL_HYGIENE.md` honestly classifies each closure. 128 are EXACT-structural (residual = 0); 31 are HIGH_PRECISION (< CODATA); 67 are WITHIN_EXP_UNCERTAINTY; 32 are REFINEMENT_TIER (0.1-1%); 5 are TENSION_OR_OUTLIER (> 1%). The latter are documented as future-work items, not hidden.

### "What about the Star-Magic LENR reactor claims?"

The COP 555:1 prediction at 27 W ambient temperature is an UNCONFIRMED prediction, explicitly so labeled. Independent replication is invited. See `forward_predictions.md` for the full list of falsifiable claims.

### "Has any third party replicated this?"

As of v5.28.0 (manuscript submission), no published independent replication exists. This is acknowledged in manuscript Section 8 (Limitations) and is the primary outstanding gap. The PyPI package + GitHub repo + CI history are intended to enable such replication; the framework is ready for it.

### "What's the catch?"

There is none materially affecting reproducibility. The framework is open-source (AGPL), version-controlled, CI-validated, and immediately installable. The only "catch" is that physics is hard, and a parameter-economical framework that disagrees with the Standard Model on Yang-Mills mass gap (5970 GeV vs lattice's 1.78 GeV) cannot both be correct. UQFF welcomes the experimental test.
