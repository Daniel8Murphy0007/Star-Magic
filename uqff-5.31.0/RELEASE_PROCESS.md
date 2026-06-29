# UQFF Star-Magic — Release Process

This document describes how to cut a new release. For end users wanting to
install UQFF, see INSTALL.md instead.

## Audience

Maintainers (currently Daniel T. Murphy, future co-maintainers). Anyone
with `push` access to the `Daniel8Murphy0007/Star-Magic` GitHub repo.

---

## Pre-release checklist

Before cutting any release, verify the following on `master`:

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
git status                                  # working tree clean
git log --oneline -5                        # review recent changes
python uqff_fidelity_tests.py               # 857/857 passing
gh run list --workflow=ci.yml --limit 1     # latest CI run = success
```

All four must succeed. If anything is red, fix BEFORE cutting the release.

---

## Versioning (Semantic Versioning 2.0.0)

UQFF follows semver:

- `MAJOR.MINOR.PATCH` — e.g., `5.27.1`
- **MAJOR**: breaking API change (removed surface, changed return contract)
- **MINOR**: new closure / surface / CLI command (backward compatible)
- **PATCH**: bug fix / docs / packaging (no API change)

Examples:
- `5.27.0 → 5.27.1` — patch: added CLI tool, docs
- `5.27.1 → 5.28.0` — minor: added new calculate_* surface
- `5.28.0 → 6.0.0` — major: changed `calculate_status_report` return shape

---

## Release procedure (step-by-step)

### Step 1 — Bump version

Edit these THREE files to the new version number:

1. `pyproject.toml`:
   ```toml
   version = "5.27.1"
   ```
2. `uqff_cli.py`:
   ```python
   _VERSION = "5.27.1"
   ```
3. `CITATION.cff`:
   ```yaml
   version: "5.27.1"
   ```

### Step 2 — Write release notes

Create `RELEASE_NOTES_v<TAG>.md` (e.g., `RELEASE_NOTES_v5.27.1.md`).

The `release.yml` workflow reads this file as the GitHub release body. If
absent, GitHub auto-generates notes from commit messages — works but less
polished.

Keep release notes under ~50 KB. Structure:

```markdown
# UQFF Star-Magic vX.Y.Z — <theme>

**Date:** YYYY-MM-DD
**Type:** Patch | Minor | Major
**Upgrade:** `pip install --upgrade uqff`

## New in vX.Y.Z

(headline new features)

## Bug fixes

(important fixes)

## What didn't change

(reassurance for risk-averse users)

## Upgrade path

(specific instructions)
```

### Step 3 — Run the fidelity gate one more time

```bash
python uqff_fidelity_tests.py
# Expected: TOTAL: 857 passed, 0 failed (or higher)
```

If anything fails, STOP. Diagnose and fix before tagging.

### Step 4 — Commit version + release notes

```bash
git add pyproject.toml uqff_cli.py CITATION.cff RELEASE_NOTES_v5.27.1.md
git commit -m "release: prepare v5.27.1

<one-line summary of what's in this release>"
git push
```

Wait for CI to go green on the version-bump commit (~2 min).

### Step 5 — Cut and push the tag

```bash
git tag -a v5.27.1 -m "UQFF Star-Magic v5.27.1 — <theme>

<2-3 bullet headline summary>"

git push origin v5.27.1
```

The tag push triggers `release.yml` automatically.

### Step 6 — Watch the release workflow

```bash
sleep 10
gh run watch $(gh run list --workflow=release.yml --limit 1 --json databaseId --jq '.[0].databaseId')
```

The pipeline runs 4 jobs:

1. **build** — sdist + wheel + twine check (~2 min)
2. **publish-pypi** — publishes to PyPI via Trusted Publishing OIDC
   - PAUSES for manual approval if you set a required reviewer on the `pypi` environment
3. **github-release** — creates GitHub release page with wheel + sdist
   attached, using `RELEASE_NOTES_v<TAG>.md` if present

### Step 7 — Verify

```bash
sleep 30
pip install --upgrade uqff
uqff version
# Should print the new version number
```

Browse:
- https://pypi.org/project/uqff/ — should show the new version + RELEASE_NOTES rendered
- https://github.com/Daniel8Murphy0007/Star-Magic/releases/tag/v5.27.1 — should show notes + downloadable wheel + sdist

### Step 8 — Append SESSION_LOG entry

```bash
cat >> SESSION_LOG.md << 'EOF'

---

## Session YYYY-MM-DD — vX.Y.Z released to PyPI

(1-paragraph summary of what shipped)
EOF
git add SESSION_LOG.md
git commit -m "log: vX.Y.Z release"
git push
```

---

## Rollback procedure

PyPI is **append-only**. Once published, a release CANNOT be edited or
deleted. The only options are:

- **Yank** the release (`pip install uqff==X.Y.Z` will warn but still work)
- **Publish a new patch** (`X.Y.Z+1`) with the fix

To yank from PyPI:

1. Go to https://pypi.org/manage/project/uqff/release/X.Y.Z/
2. Scroll to "Options" → "Yank"
3. Provide a reason

To delete a bad git tag (BEFORE PyPI publishes):

```bash
git tag -d vX.Y.Z
git push origin :refs/tags/vX.Y.Z
```

---

## Hotfix releases

For urgent bug fixes:

1. Branch from the bad-version's tag: `git checkout -b hotfix/v5.27.2 v5.27.1`
2. Apply the fix
3. Cherry-pick to master if applicable: `git checkout master && git cherry-pick <commit>`
4. Tag and release: `git tag v5.27.2 hotfix/v5.27.2 && git push origin v5.27.2`

---

## Trusted Publishing setup (one-time, per PyPI account)

If PyPI Trusted Publishing has been reset, redo:

1. Go to https://pypi.org/manage/account/publishing/
2. Add a pending publisher with these EXACT fields:
   - PyPI Project Name: `uqff`
   - Owner: `Daniel8Murphy0007`
   - Repository name: `Star-Magic`
   - Workflow name: `release.yml`
   - Environment name: `pypi`
3. Click "Add". **VERIFY** by scrolling down — the entry MUST appear in the
   "Pending publishers" table. The form silently fails ~5% of the time.

---

## GitHub Environment setup (one-time per repo)

The `release.yml` workflow uses two GitHub Environments:

- `pypi` — for production PyPI publishes (tag-triggered)
- `testpypi` — for sandbox dry-runs (workflow_dispatch with target=testpypi)

Configure at https://github.com/Daniel8Murphy0007/Star-Magic/settings/environments

**Recommended:** add a required reviewer to the `pypi` environment so the
publish step pauses for your manual approval. This stops accidental tag-pushes
from going live.

---

## Cadence

No fixed cadence. Release when there's a coherent set of changes worth
shipping. Typical patch releases every 1-4 weeks; minor releases every 1-3
months. Major releases (breaking API changes) intentionally infrequent.

Always announce in `SESSION_LOG.md`.
