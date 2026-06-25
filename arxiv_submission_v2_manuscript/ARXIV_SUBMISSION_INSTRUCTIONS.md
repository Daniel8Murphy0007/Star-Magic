# arXiv submission instructions — UQFF Star-Magic v2 manuscript

**Manuscript title:** UQFF Star-Magic: A vacuum-first unified physics framework deriving over one thousand measured observables from nine truly-independent primitives
**Author:** Daniel T. Murphy
**Source files:** This directory (`arxiv_submission_v2_manuscript/`)
**Compiled PDF:** `uqff_manuscript_v2.pdf` (61 pages, 591 KB)

---

## What this directory contains

```
abstract.tex                              ← \input'd by main
section_01_introduction.tex
section_02_primitives.tex
section_03_closure_architecture.tex
section_04_1_cosmological_constant.tex
section_04_2_magic_numbers.tex
section_04_3_holmlid_lenr.tex
section_04_4_yang_mills.tex
section_04_5_sm_spectrum.tex
section_04_6_forward_predictions.tex
section_05_statistical_hygiene.tex
section_06_provenance.tex
section_07_reproducibility.tex
section_08_limitations.tex
section_09_conclusions.tex
references.bib
uqff_manuscript_v2.bbl                    ← arXiv-compatible pre-compiled bibliography
uqff_manuscript_v2.tex                    ← main knit file
uqff_manuscript_v2.pdf                    ← reference PDF (do NOT upload to arXiv)
ARXIV_SUBMISSION_INSTRUCTIONS.md          ← this file (do NOT upload)
```

---

## Step-by-step arXiv submission

### Step 1 — Create a .tar.gz of the source files

```bash
cd C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic
tar -czf uqff_arxiv_submission_v2.tar.gz \
    arxiv_submission_v2_manuscript/abstract.tex \
    arxiv_submission_v2_manuscript/section_*.tex \
    arxiv_submission_v2_manuscript/references.bib \
    arxiv_submission_v2_manuscript/uqff_manuscript_v2.bbl \
    arxiv_submission_v2_manuscript/uqff_manuscript_v2.tex
```

The resulting `uqff_arxiv_submission_v2.tar.gz` is what you upload to arXiv.

### Step 2 — Go to arXiv and log in

URL: https://arxiv.org/submit

If you do not yet have an arXiv account, create one. You may need an endorsement from an existing arXiv author to submit to physics.gen-ph; see https://arxiv.org/help/endorsement for the procedure.

### Step 3 — Start a new submission

- **License:** arXiv non-exclusive (compatible with your AGPL-3.0 dual license; no commercial conflict)
- **Primary archive:** physics
- **Primary subject class:** physics.gen-ph (General Physics) — the appropriate category for parameter-economy frameworks
- **Cross-list:** math-ph (Mathematical Physics) — given the Lean 4 scaffold and BIC analysis
- **(Optional)** hep-ph (Phenomenology) — given the Standard Model spectrum closures

### Step 4 — Upload the .tar.gz

Drag-drop or browse-and-upload `uqff_arxiv_submission_v2.tar.gz`.

arXiv's processing pipeline will:
1. Untar the source files
2. Run pdflatex + bibtex + pdflatex + pdflatex (automatic)
3. Produce a PDF and a HTML abstract page
4. Show you the preview

Verify the preview matches your local PDF (61 pages, with the boxed BIC equation in §5).

### Step 5 — Metadata

- **Title:** UQFF Star-Magic: A vacuum-first unified physics framework deriving over one thousand measured observables from nine truly-independent primitives
- **Authors:** Daniel T. Murphy
- **Abstract:** Copy verbatim from `abstract.tex` (strip the `\begin{abstract}` / `\end{abstract}` and `\textbf` markup, keep the keywords line at the bottom)
- **Comments:** "61 pages, 13 sections, AGPL-3.0 + commercial dual-licensed source at github.com/Daniel8Murphy0007/Star-Magic ; PyPI: pip install uqff (v5.29.1+); 866-test fidelity gate."
- **MSC class:** (optional) 81-04 (Computational software, Mathematical physics)
- **PACS:** (optional) 04.50.Kd (Modified theories of gravity), 11.10.-z (Field theory), 12.10.-g (Unified field theories)
- **DOI:** (leave blank — Zenodo DOI gets added later via journal cross-reference)

### Step 6 — Submit + record the arXiv ID

After submission, arXiv typically assigns the new submission an ID within a few hours (e.g., `arXiv:2606.XXXXX`). You will receive an email confirmation.

**Save the arXiv ID** --- you will need to:
- Update CITATION.cff with the arXiv identifier
- Add the arXiv link to README.md
- Include the arXiv ID in the cover letter to *Foundations of Physics*
- Update the manuscript's title page footnote (if revising before Foundations of Physics submission)

### Step 7 — Follow-up

- The arXiv preprint becomes citable immediately after announcement (usually 1-2 business days).
- *Foundations of Physics* accepts arXiv-preprinted manuscripts (most physics journals do); the cover letter should cite the arXiv ID.
- After Foundations of Physics acceptance (if accepted), the journal version is the canonical citation; the arXiv version remains as the working preprint.

---

## Endorsement note (if needed)

arXiv requires endorsement for first-time submitters in physics.gen-ph. If you do not have an arXiv-author co-author who can endorse, the standard route is:

1. Email an arXiv author who has previously published in physics.gen-ph (e.g., one of the suggested reviewers from the cover letter, if appropriate). Request endorsement, briefly describing the work.
2. They submit an endorsement code via the arXiv endorsement system.
3. You enter the code at submission time.

For a self-contained framework with a published PyPI package, this endorsement step is usually straightforward. Allow 1-2 weeks lead time.

---

## What NOT to include in the arXiv submission

- The compiled `.pdf` (arXiv compiles its own from your `.tex` sources)
- The `.aux`, `.log`, `.out`, `.toc`, `.synctex.gz` files (auto-generated)
- The `.gitignore`, `README.md`, this `ARXIV_SUBMISSION_INSTRUCTIONS.md` file
- Any backup files (`*.backup`, `*.orig`, `*~`)
- The full Star-Magic repository (only the manuscript source goes to arXiv; the code stays at the GitHub URL cited inside the paper)

---

## Verification before submission

```bash
# Confirm tar contents
tar -tzf uqff_arxiv_submission_v2.tar.gz
# Should show 16 files: abstract + 13 sections + references + bbl + main

# Confirm size
ls -la uqff_arxiv_submission_v2.tar.gz
# Should be ~100-200 KB (arXiv max is 50 MB; we're well under)

# (Optional) Test-compile locally to mimic arXiv's pipeline
mkdir test_compile && cd test_compile
tar -xzf ../uqff_arxiv_submission_v2.tar.gz --strip-components=1
pdflatex uqff_manuscript_v2.tex
bibtex uqff_manuscript_v2
pdflatex uqff_manuscript_v2.tex
pdflatex uqff_manuscript_v2.tex
# Should produce uqff_manuscript_v2.pdf identical to your reference PDF
```

