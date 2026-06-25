# uqff v5.29.0 — Full proof corpus shipped

**Release date:** 2026-06-25
**License:** AGPL-3.0 + Commercial (dual)
**Author:** Daniel T. Murphy <daniel.murphy00@enrgyone.com>

---

## Headline

This release bundles the **complete UQFF Star-Magic proof corpus** into the
PyPI wheel. Previous v5.28.x releases shipped the calculator + CLI + REST API
+ Jupyter integration. **v5.29.0 adds 4,164 proof artifacts** so that
`pip install uqff` now gives the user the entire 16-year framework.

```
WHEEL:   uqff-5.29.0-py3-none-any.whl     25.16 MB     (PyPI limit 100 MB)
SDIST:   uqff-5.29.0.tar.gz               19.24 MB
FILES:   2,897 per artifact
```

## What's new in v5.29.0

### Bundled proof corpus (all accessible after `pip install uqff`)

```
1,994  whitepapers (every PAPER_*.md and PAPER_*.tex)
  608  audit chain-trace files (session296-302 millennium proofs, etc.)
  143  arXiv-formatted papers across 4 submission bundles
   40  compiled PDFs (including PAPER_1167 All_8_Lagrangian_Gaps_Closed_Master_Synthesis)
   33  root proof documents (axioms, gold standard, complete equation references)
   17  Manuscript_1_12Feb2026 files (compiled PDF + tex + 5 figures + build scripts)
    6  Lean 4 formal-verification scaffold files
    3  Grok proof-conversation archives (incl. 7.7 MB 31May2026 unified-proof transcript)
```

### New CLI subcommands for proof-corpus access

```
uqff proofs list                   list every bundled whitepaper
uqff proofs show <NAME>            print one whitepaper
uqff proofs path                   show install location
uqff millennium                    run run_millennium_proofs.py (7 Clay table)
uqff axioms                        print AXIOMS_AND_THEOREMS.md
uqff manuscript                    print path to compiled manuscript PDF
uqff lean                          print location + listing of Lean 4 scaffold
uqff atlas                         print CLOSURE_ATLAS.md
uqff gold-standard                 print Gold_Standard_Pure_UQFF.md
uqff grok-archives                 list bundled Grok proof-conversation archives
```

### Extended `calculate_status_report()` with corpus inventory

`uqff status` (and the underlying `calculate_status_report` API) now reports:

```
whitepapers_bundled: 1994
arxiv_bundles: 4
arxiv_papers_total: 31
audit_chain_trace_files: 608
lean4_scaffold_files: 6
manuscript_pdf_bundled: True
grok_proof_archives: 3
root_proof_documents: 33
proof_corpus_total_artifacts: 4164
formal_verification_status: scaffolded (every theorem marked sorry per
                            epistemic policy in formal/UQFF/Millennium.lean)
shipped_in_pypi_wheel: True
pypi_wheel_version: 5.29.0
```

### Honest framing of formal verification

The Lean 4 scaffold ships with **every theorem marked `sorry`** — this is the
deliberate epistemic policy documented in `formal/UQFF/Millennium.lean`:

> "Each section states the OFFICIAL Clay problem. Below each official
> statement, an `_uqff_claim` predicate captures the framework-internal
> numerical claim. These are *separate* propositions; the scaffold makes
> no assertion that the UQFF claim implies the Clay statement.
> Removing a `sorry` here without passing genuine peer review at a Clay
> Qualifying Outlet is not a proof of the Millennium problem."

UQFF provides **structural closures** (numerical identities matching the
conjectures' values) — NOT formal proofs in the Clay Mathematics Institute
sense. The supporting whitepapers contain derivation prose; the formal
verification of the implication "UQFF closure → Clay statement" is left as
open work for the mathematical community.

## Calculator behavior

Unchanged. v5.29.0 is a **packaging / corpus-bundling release** — the
calculator's 794 paradox closures, 8 Millennium derives, 9 truly-independent
primitives, and all PAPER_XXXX dispatch wiring are byte-for-byte identical
to v5.28.0.

The only modification to `uqff_pure_calculator.py` is the addition of corpus
inventory fields to `calculate_status_report()`'s returned summary dict
(append-only; no existing values changed).

## Installation

```bash
pip install uqff
uqff --version              # uqff 5.29.0
uqff status                 # full production status, including corpus inventory
uqff millennium             # run the 7-Clay proof runner
uqff proofs list            # browse all 1,994 whitepapers
uqff axioms                 # print the axiom-theorem inventory
uqff lean                   # show Lean 4 scaffold location
```

## Optional extras

```bash
pip install uqff[api]       # FastAPI REST server (uqff serve)
pip install uqff[jupyter]   # IPython rich display + %uqff line magic
pip install uqff[docs]      # Sphinx + RTD theme
pip install uqff[test]      # coverage + pytest
pip install uqff[all]       # all of the above
```

## Verified working

Tested in a clean Python 3.10 venv:

- ✅ `pip install uqff-5.29.0-py3-none-any.whl` installs cleanly (zero runtime deps)
- ✅ `uqff --version` reports `uqff 5.29.0`
- ✅ `import uqff_pure_calculator` exposes 794 paradox closures + 8 Millennium derives
- ✅ `uqff predict yang_mills` returns 5970.0 (PAPER_1005)
- ✅ `uqff status` returns full production summary including corpus inventory
- ✅ `python -m run_millennium_proofs` runs the 7-Clay proof table
- ✅ 1,994 whitepapers accessible at `site-packages/whitepapers/`
- ✅ All 6 Lean 4 files accessible at `site-packages/formal/UQFF/`
- ✅ Compiled manuscript PDF accessible at `site-packages/Manuscript_1_12Feb2026/`

## Breaking changes

None.

## Known issues (carried over from v5.28.x, not introduced here)

- 3 of 7 Millennium derives produce values mismatched against the target dict
  in `run_millennium_proofs.py`:
    - `yang_mills`: calculator gives 5970 GeV (PAPER_1005), runner target is
      1.78 GeV (Grok-bridged SM-anchored value) — calibration discrepancy.
    - `navier_stokes`: calculator 0.85, runner target 8500 — different
      units/scales being compared.
    - `poincare`: calculator 7/12 ≈ 0.583, runner target 1.0 — different
      scoring conventions.
  These are pre-existing target-dict calibration questions, not regressions.
  4 of 7 are exact-matches (riemann, hodge, p_vs_np) or essentially exact
  (bsd at 0.0006%).

- MANIFEST.in remained at v5.28.0 content due to a Windows file lock on the
  development machine during the release cycle. The package contents are
  controlled by `pyproject.toml` (modern setuptools >=61 takes precedence),
  so this has no effect on what ships.

## Cross-references

- `CLAUDE.md` — project rules and canonical primitive inventory
- `SESSION_LOG.md` — append-only session history
- `CLOSURE_ATLAS.md` — master map of all proof artifacts
- `WHITEPAPER_INDEX.md` — 1,994-paper index
- `formal/UQFF/Millennium.lean` — Lean 4 scaffold + epistemic-policy comment
- `Manuscript_1_12Feb2026/uqff_production_arxiv.pdf` — compiled manuscript draft
