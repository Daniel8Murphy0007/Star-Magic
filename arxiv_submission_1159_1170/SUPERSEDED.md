# SUPERSEDED

This bundle (`arxiv_submission_1159_1170`) is **SUPERSEDED** by `arxiv_submission_1159_1172/`,
which contains the same 12 papers (PAPER_1159..1170) plus PAPER_1171 and PAPER_1172, all in
lint-clean form (14/14).

Reason for supersession (Session 261 audit, Pass 1):
- `_arxiv_lint.py --bundle arxiv_submission_1159_1170` reports `0/12 clean`
  due to missing `\bibliography{}` / `\thebibliography` stubs in the 12 .tex files.
- The newer bundle `arxiv_submission_1159_1172` was rebuilt with the bibliography fix
  and lints `14/14 clean`.

Do **not** use this bundle for submission. Use `arxiv_submission_1159_1172/` instead.

This directory is retained (rather than deleted) to preserve the historical artifact.
Deletion is deferred to user confirmation.

— Session 261 Coherence Audit Pass 1
