# pdf2/ — arXiv-Compliant UQFF Whitepaper PDF Corpus

This folder contains PDF renderings of every UQFF whitepaper following the arXiv publishing rules.

## Build

Run from repo root:

```bash
python _build_pdf2_arxiv_compliant.py
```

Or a subset:

```bash
python _build_pdf2_arxiv_compliant.py --limit 20                 # first 20
python _build_pdf2_arxiv_compliant.py --pattern "PAPER_10*"      # PAPER_100-109 range
python _build_pdf2_arxiv_compliant.py --pattern "PAPER_11??"     # PAPER_1100-1199
python _build_pdf2_arxiv_compliant.py --jobs 4                   # 4 parallel pandoc procs
python _build_pdf2_arxiv_compliant.py --force                    # rebuild every PDF
```

## arXiv publishing rules honored

| Rule                              | Implementation                                             |
|-----------------------------------|------------------------------------------------------------|
| Text-searchable output            | pandoc → LaTeX → pdflatex/xelatex/lualatex (no raster)     |
| Embedded fonts                    | fontspec + DejaVu via lualatex (pdf_header.tex)            |
| Standard paper size               | letterpaper (0.9-in margins)                               |
| Metadata (title/author/subject)   | hyperref \pdftitle etc., pulled from YAML frontmatter      |
| Reasonable file size              | ~100–500 KB for short papers                               |
| Source-available                  | .md and .tex sources under whitepapers/ preserved          |
| No unusual PDF features           | pandoc default output, no JS/video/forms                   |
| Consistent template               | shared pdf_header.tex header (Star-Magic project standard) |

## Idempotency

If `pdf2/PAPER_NNN.pdf` exists and is newer than `whitepapers/PAPER_NNN.md`, the file is skipped. Delete a PDF to force a rebuild, or run with `--force` for all.

## Log

Each run writes `pdf2/_build_log.txt` recording pass/fail for every paper. Failures capture the last ~800 chars of pandoc/latex stderr for diagnosis.

## Failures to expect

Some source files have pre-existing content issues that block any pipeline:

- Mismatched math delimiters (e.g., `$foo$ = bar $$` in one file — double-dollar closes an already-open dollar-mode)
- Non-ASCII characters not covered by DejaVu font subset
- Malformed YAML frontmatter
- Extremely long math expressions that overflow default `\displaystyle`

The script logs these and continues. Fix the source markdown in `whitepapers/` and rerun — only the failed papers rebuild.

## Sibling folders

- `pdf/` — the previous production PDF corpus (1,296 files), preserved as historical reference
- `pdf_backup_pandoc_2026-05-11/`, `pdf_backup_pandoc_2026-05-12/`, `pdf_backup_pandoc_2026-06-05/` — earlier snapshots

`pdf2/` is the fresh, arXiv-rule-driven regeneration.

## Author

Daniel T. Murphy · Star-Magic Research Program · daniel.murphy00@enrgyone.com

## License

Follows parent repository dual license (AGPL-3.0-or-later / Commercial).
