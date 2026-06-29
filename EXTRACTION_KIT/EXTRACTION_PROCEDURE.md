# Aetheric-Propulsion Extraction Procedure (Phase Step 7)

**Author:** Daniel T. Murphy
**Source:** Star-Magic v5.30.0+ (this repo)
**Target:** https://github.com/Daniel8Murphy0007/Aetheric-Propulsion (to be created)

This document is the procedure for extracting the assimilation geometry from Star-Magic
into a standalone `aetheric-propulsion` pip package. It can be executed at any time after
the EXTRACTION_KIT/ is in place — the kit stays in sync with Star-Magic via a refresh script.

---

## What gets extracted

Bundled into the new package (`aetheric_propulsion/`):

| File in Star-Magic | Becomes in Aetheric-Propulsion |
|---|---|
| `uqff_pure_calculator.py` | `aetheric_propulsion/calculator.py` |
| `assimilation_dispatch.py` | `aetheric_propulsion/assimilation_dispatch.py` |
| `qcalcgeom_solver.py` | `aetheric_propulsion/qcalcgeom_solver.py` |
| `provenance_recorder.py` | `aetheric_propulsion/provenance_recorder.py` |
| `_build_overdetermination_views.py` | `aetheric_propulsion/_build_overdetermination_views.py` |
| `geometry_backends/` | `aetheric_propulsion/geometry_backends/` |
| `numeric_backends/` | `aetheric_propulsion/numeric_backends/` |

Data files bundled (read-only docs that ship with the wheel):

| File in Star-Magic | Becomes in Aetheric-Propulsion |
|---|---|
| `OVERDETERMINATION_MAP.csv` | `data/OVERDETERMINATION_MAP.csv` |
| `OVERDETERMINATION_WIDE.csv` | `data/OVERDETERMINATION_WIDE.csv` |
| `OVERDETERMINATION_MAP.md` | `data/OVERDETERMINATION_MAP.md` |
| `ASSIMILATION_GEOMETRY_ATLAS.md` | `data/ASSIMILATION_GEOMETRY_ATLAS.md` |

What stays in Star-Magic: **everything**. The extraction is additive at this stage —
both repos ship the same core science during the peer-review phase (per the
"Maximum access for academic peer-review and NASA-Roses grant panels" mandate).

---

## Step-by-step procedure

### Step 7.1 — Create the GitHub repo

1. Go to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion (URL reserved).
2. Initialize as a public repository. Match the Star-Magic license setup:
   - Add a LICENSE file (will be overwritten in step 7.4 with AGPL-3.0)
   - No README (will be written in step 7.4)
   - .gitignore: Python template

### Step 7.2 — Clone locally and refresh the kit

```bash
# Clone the new Aetheric-Propulsion repo
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007
git clone https://github.com/Daniel8Murphy0007/Aetheric-Propulsion.git

# Refresh the EXTRACTION_KIT/new_repo_layout/ from current Star-Magic state
cd Star-Magic
python3 EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py
# (no target dir argument -> refresh only)
```

### Step 7.3 — Migrate files to the new repo

```bash
# Run migration script with target dir
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
python3 EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py \
        /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Aetheric-Propulsion
# Expected output: Migrated N files to /c/Users/tmsjd/...
```

### Step 7.4 — Verify, commit, push

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Aetheric-Propulsion

# Verify the package builds + smoke test passes
python3 -m pip install -e .
python3 tests/test_smoke.py
# Expected: SMOKE: PASS

# Commit + push
git add .
git commit -m "Aetheric-Propulsion v0.1.0 — extracted from Star-Magic v5.30.0"
git push origin main
```

### Step 7.5 — Configure PyPI Trusted Publishing

1. Go to https://pypi.org/manage/account/publishing/
2. Add a new pending publisher:
   - **PyPI project name:** `aetheric-propulsion`
   - **Owner:** `Daniel8Murphy0007`
   - **Repository:** `Aetheric-Propulsion`
   - **Workflow:** `release.yml`
   - **Environment:** (leave blank)
3. Save.

### Step 7.6 — First PyPI release

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Aetheric-Propulsion
git tag -a v0.1.0 -m "v0.1.0 — first extract from Star-Magic v5.30.0"
git push origin v0.1.0
# release.yml triggers, publishes to pypi.org/project/aetheric-propulsion/
```

### Step 7.7 (optional, deferred) — Star-Magic cross-link

After Aetheric-Propulsion is live on PyPI, update Star-Magic README to add a "Related
packages" section linking to it. Do NOT remove any files from Star-Magic during
peer-review phase.

---

## Future re-syncs

When Star-Magic's assimilation dispatch grows (new rounds inject more observables),
re-run the migration script to keep Aetheric-Propulsion in sync:

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
python3 EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py \
        /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Aetheric-Propulsion
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Aetheric-Propulsion
# Run tests, commit, tag, push for the next release
```

This staged-kit architecture keeps the two repos loosely coupled while maintaining
single-source-of-truth in Star-Magic.

---

## Verification harness

`test_extraction_kit.py` in Star-Magic verifies the kit's contents are well-formed and
the migration script runs idempotently. Run after any change to the kit:

```bash
python3 test_extraction_kit.py
# Expected: PHASE STEP 7 SUCCESS CRITERION MET.
```
