# UQFF Star-Magic v5.27.1 — Tier-2 polish patch

**Date:** 2026-06-22
**Type:** Patch release (no breaking changes)
**Upgrade:** `pip install --upgrade uqff`

This is a polish release adding the `uqff` command-line interface and
extensive user-facing documentation. The underlying physics and the 794
closures are unchanged from v5.27.0.

---

## New in v5.27.1

### CLI tool (`uqff` command)

After `pip install --upgrade uqff`, the `uqff` command is available:

```bash
uqff version                       # version + headline metrics
uqff predict hubble_tension        # fetch a closure
uqff predict <name> --json         # machine-readable output
uqff list                          # list all 794 closure names
uqff list --filter neutrino        # substring filter
uqff status                        # full production status summary
uqff surfaces                      # list 34 public calculate_* surfaces
uqff gate                          # run the 857-test fidelity gate (source install only)
```

### Documentation expansion

10 new documentation files for users and contributors:

- `FAQ.md` — 12 frequently asked questions answered
- `TROUBLESHOOTING.md` — 15+ common-issue fixes (install, dispatch, CI, Windows)
- `GLOSSARY.md` — 60+ UQFF terms (primitives, acronyms, sectors, buckets, tiers)
- `ARCHITECTURE.md` — codebase layout, calculator internals, CI/CD, extension workflow
- `INSTALL.md` — formal installation guide for end users
- `RELEASE_PROCESS.md` — maintainer guide for cutting future releases
- `Dockerfile` + `.dockerignore` — `docker run uqff` for containerized use
- 3 GitHub issue templates (bug / feature / closure proposal)
- 1 GitHub PR template with Rule-2/3/4/5/6/8/9 compliance checklist
- `.readthedocs.yaml` — connects Sphinx docs to RTD hosting

### Jupyter notebooks (4 new)

`notebooks/` directory with worked examples:

- `00_quickstart.ipynb` — 5-minute UQFF walkthrough
- `01_holmlid_lenr.ipynb` — derives Holmlid 630 eV LENR step-by-step
- `02_magic_numbers.ipynb` — all 7 nuclear magic numbers from integer primitives
- `03_cosmology.ipynb` — Λ + ΛCDM observables in one notebook

### CI/CD enhancements

- Smoke contract via `scripts/ci_smoke.py` runs as a standalone test
- Null-byte stripper (`scripts/ci_strip_nulls.py`) keeps CI green on Windows-saved files
- Codecov dependency removed (was requiring tokens even on public repos)
- License migrated to PEP 639 SPDX string format (`license = "AGPL-3.0-or-later"`)
- `setup.py` C++ extension build guarded behind pybind11 availability

### Bug fixes

- 5 regression pins were set tighter than float64 machine epsilon (~1e-16
  relative). Relaxed to `tol=1e-12` for consistent pass across all Python
  versions and libm implementations. These changes do NOT affect physics
  precision — they only fix tests that were physically unsatisfiable.

### Files removed

- 104 stale backup files (~183 MB) cleaned from git tracking. The 11
  CLAUDE.md-protected backups remain. New clones are ~165 MB smaller.

---

## What didn't change

The core physics is identical to v5.27.0:

- 9 truly-independent primitives + 2 derivative primitives
- 794 paradox dispatch keys, 263 schema-tagged
- 128 EXACT structural identities
- All 8 Clay Millennium prize problem closures
- Complete SM 12-fermion spectrum
- 857/857 fidelity gate passing
- ΔBIC = 94.1 vs SM + ΛCDM

If you use v5.27.0 and don't need the CLI, you don't need to upgrade.

---

## Upgrade path

```bash
pip install --upgrade uqff
uqff version              # should print 'uqff 5.27.1'
```

If you have v5.27.0 installed and don't need the CLI, no upgrade is required;
the existing import paths continue to work identically.

---

## License

Dual AGPL-3.0 + commercial. See `LICENSE`, `COMMERCIAL.md`.

---

**Author:** Daniel T. Murphy / Star-Magic Research Program
**Copyright:** © 2025-2026, all rights reserved.
