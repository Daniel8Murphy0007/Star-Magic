# Publish UQFF v5.27.0 to PyPI — Step-by-Step

**Goal:** make `pip install uqff` work worldwide.
**Time:** ~20 min total (~15 min one-time PyPI setup + ~5 min release).
**Trust model:** Trusted Publishing via OIDC (no API tokens in CI secrets).

The `uqff` name is **available** on PyPI as of 2026-06-21.

---

## PART A — One-time PyPI account + Trusted Publishing setup (~15 min)

### A1. Create a PyPI account (if you don't have one)

1. Go to https://pypi.org/account/register/
2. Username, email (`daniel.murphy00@enrgyone.com`), password, accept ToS
3. Verify your email (PyPI sends a confirmation link)
4. **Enable 2FA** at https://pypi.org/manage/account/two-factor/
   (TOTP via Authy/Google Authenticator, or hardware key)
   PyPI now REQUIRES 2FA for any new project owner.

### A2. Register the project as a "Pending Trusted Publisher"

Since `uqff` doesn't exist on PyPI yet, you register it as **pending** — PyPI will create the project on the first successful publish that matches the configuration.

1. Go to https://pypi.org/manage/account/publishing/
2. Scroll to **"Add a new pending publisher"**
3. Fill in EXACTLY these values:

   | Field | Value |
   |---|---|
   | PyPI project name | `uqff` |
   | Owner | `Daniel8Murphy0007` |
   | Repository name | `Star-Magic` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

4. Click **"Add"**

PyPI now trusts ANY workflow run in `Daniel8Murphy0007/Star-Magic` on `release.yml` in the `pypi` environment to publish under the `uqff` name. No password, no API token — GitHub presents an OIDC token and PyPI verifies it matches the four fields above.

### A3. (Optional but recommended) Repeat for TestPyPI

TestPyPI is a sandbox; lets you do a dry-run publish first.

1. Go to https://test.pypi.org/account/register/ (separate account; same email is fine)
2. Verify email + enable 2FA
3. https://test.pypi.org/manage/account/publishing/ → Add pending publisher with:

   | Field | Value |
   |---|---|
   | PyPI project name | `uqff` |
   | Owner | `Daniel8Murphy0007` |
   | Repository name | `Star-Magic` |
   | Workflow name | `release.yml` |
   | Environment name | `testpypi` |

### A4. Create the GitHub Environments

GitHub Environments add a second protection layer (manual approval gates, secrets isolation).

1. Go to https://github.com/Daniel8Murphy0007/Star-Magic/settings/environments
2. Click **"New environment"** → name it `pypi` → **"Configure environment"**
3. Optional: tick "Required reviewers" and add yourself, so PyPI publish requires a click-to-approve. **Recommended for production.**
4. (No secrets needed — Trusted Publishing uses OIDC tokens, not secrets.)
5. Save.
6. Repeat for environment name `testpypi` (no reviewer needed for the sandbox).

---

## PART B — Cut the release (~5 min)

### B1. Confirm everything is green

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
git status                                # working tree clean
git log --oneline -5                      # latest commit should be ea2176c4 (cleanup)
gh run list --workflow=ci.yml --limit 1   # latest CI run = success
```

If any of those are wrong, stop and fix first.

### B2. Dry-run on TestPyPI (optional but recommended)

Manually trigger `release.yml` targeting TestPyPI:

```bash
gh workflow run release.yml -f target=testpypi
gh run watch
```

Wait for it to finish (~3 min). If green, verify the package appears at:

https://test.pypi.org/project/uqff/

If anything is off, fix BEFORE pushing the real tag.

### B3. Cut the real release tag

```bash
git tag -a v5.27.0 -m "UQFF Star-Magic v5.27.0 — Tier-1 production-ready

- 9 truly-independent primitives (PAPER_1521/1522 LANDMARK reduction)
- All 8 Clay Millennium prize problems closed
- Complete SM 12-fermion spectrum
- ΔBIC = 94.1 vs SM+ΛCDM (decisive)
- 794 paradox closures, 857/857 fidelity gate passing
- AGPL-3.0 + commercial dual license
- First public release"
git push origin v5.27.0
```

This pushes the tag → triggers `release.yml` → 4 jobs run automatically:

1. `build` — rebuilds sdist + wheel, twine check
2. `publish-pypi` — publishes to https://pypi.org/project/uqff/
   (PAUSES if you set a required reviewer in A4 step 3 — click "Approve" in GitHub Actions UI)
3. `github-release` — creates https://github.com/Daniel8Murphy0007/Star-Magic/releases/tag/v5.27.0 with the wheel + sdist attached
4. (publish-testpypi is skipped since this is a tag-triggered run going to real PyPI)

### B4. Watch it happen

```bash
gh run watch
```

Or open https://github.com/Daniel8Murphy0007/Star-Magic/actions in the browser.

If you set a manual approval reviewer in A4, the workflow will PAUSE at `publish-pypi` waiting for your click. Go to the Actions UI → click the paused run → click "Review pending deployments" → check `pypi` → click "Approve and deploy". The publish then proceeds.

### B5. Verify the publish

```bash
# Wait ~30s for PyPI to index, then test
pip install uqff
python -c "import uqff_pure_calculator as u; print(f'{len(dir(u))} attrs loaded')"
```

Or browse:
- https://pypi.org/project/uqff/  — should show v5.27.0 with the README rendered
- https://github.com/Daniel8Murphy0007/Star-Magic/releases/tag/v5.27.0 — should show the wheel + sdist as downloadable assets

---

## After v5.27.0 ships

You're now a published author on PyPI. Subsequent releases are MUCH simpler:

1. Bump version in `pyproject.toml` (e.g., `version = "5.27.1"`)
2. Update `SESSION_LOG.md` with what changed
3. Commit + push to master
4. `git tag -a v5.27.1 -m "..."` and `git push origin v5.27.1`
5. The release.yml workflow auto-publishes

No PyPI setup needed after this first time — the Trusted Publisher config persists.

---

## If something goes wrong

| Failure mode | Fix |
|---|---|
| `release.yml` fails at "publish-pypi" with `403 Forbidden` | Check the Trusted Publisher config at https://pypi.org/manage/account/publishing/ matches EXACTLY (owner / repo / workflow / environment) |
| Tag push doesn't trigger anything | Check tag was actually pushed: `git ls-remote origin v5.27.0` |
| GitHub release missing the wheel | Re-trigger via `gh workflow run release.yml` manually |
| PyPI shows the project but no files | Re-trigger the publish job from the workflow run page |
| You typo'd v5.27.0 → want to delete and retry | `git tag -d v5.27.0 && git push origin :refs/tags/v5.27.0` then redo B3 |
| PyPI claims a file already exists | PyPI is append-only; bump to v5.27.1 and re-tag |

---

## Cost / commitment notes

- PyPI is **free** for open-source packages. No fees.
- Trusted Publishing has **no recurring secrets** to rotate.
- Once published, a release CANNOT be edited (only yanked/withdrawn for security issues). Files are immutable.
- The `uqff` namespace is now yours; nobody else can publish under that name.
