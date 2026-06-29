# Zenodo DOI Setup (Tier-4 O3)

**Purpose:** Get a permanent citable DOI for the UQFF GitHub repository.
**Time:** ~15 minutes one-time setup; subsequent releases auto-get a DOI.
**Cost:** Free.
**Why:** A DOI (Digital Object Identifier) gives your work a permanent,
citable address. Journals, citation managers (Zotero, Mendeley), and academic
databases (Web of Science, Google Scholar) recognize DOIs. The peer-review
submission cover letter benefits greatly from including a DOI.

---

## What you get

After this 15-min setup, every GitHub release of the Star-Magic repository
automatically gets a permanent DOI on Zenodo (operated by CERN). Example:

```
DOI: 10.5281/zenodo.XXXXXXX
Citable as: Murphy, D. T. (2026). UQFF Star-Magic v5.28.0. Zenodo.
            https://doi.org/10.5281/zenodo.XXXXXXX
```

Add this DOI to:
- `CITATION.cff` (the `identifier` field — already partly set up)
- `README.md` (top of file, as a badge)
- `pyproject.toml` (`project.urls` section)
- Manuscript title page
- arXiv preprint metadata
- All future release notes

---

## Step-by-step setup

### Step 1 — Create a Zenodo account

1. Go to https://zenodo.org/signup
2. Choose **"Sign up with GitHub"** (easiest — uses your existing GitHub identity)
3. Authorize Zenodo to read your GitHub repositories
4. Verify your email if prompted

### Step 2 — Enable Zenodo integration for Star-Magic

1. Go to https://zenodo.org/account/settings/github/
2. You'll see a list of all your GitHub repositories
3. Find **`Daniel8Murphy0007/Star-Magic`** in the list
4. Toggle the switch from **OFF** to **ON**
5. That's it. Zenodo is now watching for new GitHub releases on this repo.

### Step 3 — Trigger the first DOI

Zenodo only creates DOIs for *new* GitHub releases AFTER the toggle was switched on. Your existing v5.27.0, v5.27.1, v5.27.2, v5.28.0 releases will NOT auto-get DOIs.

To trigger one, cut a new tag:

```bash
cd /c/Users/tmsjd/source/repos/Daniel8Murphy0007/Star-Magic
git tag -a v5.28.1 -m "v5.28.1 — Zenodo DOI activation"
git push origin v5.28.1
```

This triggers two things in parallel:
1. The existing release.yml workflow publishes v5.28.1 to PyPI
2. Zenodo detects the new GitHub release and mints a DOI

Within ~5 minutes, the DOI appears on the Zenodo page for the repo.

### Step 4 — Find your DOI

1. Go to https://zenodo.org/account/settings/github/
2. Click on **`Daniel8Murphy0007/Star-Magic`**
3. You'll see a "Releases" list with each release's DOI
4. The most recent release (v5.28.1) will show its DOI like `10.5281/zenodo.XXXXXXX`
5. There's also a **concept DOI** (one for ALL versions of the repo combined) — this is the one you usually want in citations because it always resolves to the latest version

Copy both DOIs.

### Step 5 — Add the DOI to your project

Update `CITATION.cff`:

```yaml
identifiers:
  - description: This is the archived snapshot of version 5.28.1 of UQFF.
    type: doi
    value: 10.5281/zenodo.XXXXXXX        # version-specific DOI
  - description: This DOI represents all versions, and will always resolve to the latest one.
    type: doi
    value: 10.5281/zenodo.YYYYYYY        # concept DOI (parent)
```

Add a Zenodo badge to the top of `README.md`:

```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.YYYYYYY.svg)](https://doi.org/10.5281/zenodo.YYYYYYY)
```

Add to `pyproject.toml`:

```toml
[project.urls]
DOI = "https://doi.org/10.5281/zenodo.YYYYYYY"
```

Commit and push these changes.

---

## After setup — auto-DOI on every release

Once configured, every future `git push origin vX.Y.Z` tag automatically:

1. Triggers `release.yml` → publishes to PyPI
2. Triggers GitHub release creation
3. Zenodo detects the GitHub release → mints a new DOI within minutes

You don't need to do anything per-release. Just tag and push.

---

## How citations look once DOI is in place

### BibTeX

```bibtex
@software{murphy_uqff_2026,
  author       = {Murphy, Daniel T.},
  title        = {UQFF Star-Magic: A vacuum-first unified physics framework},
  year         = 2026,
  version      = {5.28.1},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.YYYYYYY},
  url          = {https://doi.org/10.5281/zenodo.YYYYYYY},
}
```

### Plain text (APA-like)

```
Murphy, D. T. (2026). UQFF Star-Magic: A vacuum-first unified physics framework
(Version 5.28.1) [Software]. Zenodo. https://doi.org/10.5281/zenodo.YYYYYYY
```

### Manuscript footnote

```
The framework is open-source under AGPL-3.0 + commercial dual license at
github.com/Daniel8Murphy0007/Star-Magic. Permanent version-specific archive
at doi.org/10.5281/zenodo.YYYYYYY. Live PyPI package: pip install uqff.
```

---

## Why this matters for peer review

When you submit to Foundations of Physics or any other journal:

- The cover letter cites: arXiv preprint + Zenodo DOI + PyPI version
- Reviewers can permanently link to the EXACT version of the code being reviewed
- If the code changes during peer review, the v5.28.1-specific DOI doesn't change — preserving review continuity
- The journal's published reference list will include the DOI, making the work permanently findable

Without a DOI:
- Citations point to git commit hashes (uglier; not recognized by standard citation tools)
- No guaranteed long-term URL (Zenodo is CERN-backed; GitHub could theoretically go down)
- No standardized metadata for academic databases

---

## Costs

- Zenodo: FREE forever (CERN-funded scientific archive)
- GitHub releases: FREE (already in your workflow)
- Setup: 15 minutes one-time
- Per-release: 0 minutes (automatic)

---

## Cross-references

- `CITATION.cff` — update with DOI fields after Step 4
- `README.md` — add Zenodo badge
- `pyproject.toml` — add DOI to project.urls
- `PEER_REVIEW_SUBMISSION_PLAN.md` — reference Zenodo DOI in cover letter
- `COVER_LETTER_TEMPLATE.md` — add `[DOI: 10.5281/zenodo.YYYYYYY]` to signature block
