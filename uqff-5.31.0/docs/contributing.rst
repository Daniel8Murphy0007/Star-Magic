Contributing
============

UQFF Star-Magic is currently maintained by Daniel T. Murphy. External
contributions are welcome subject to the following constraints.

Required reading before any contribution
-----------------------------------------

1. ``CLAUDE.md`` — project rules. Read this **first**. It contains 12
   non-negotiable rules including:

   - Rule 2: DO NOT REVERT canonical primitives
   - Rule 3: NO NARRATIVE (no docstrings, no comments, no metadata strings)
   - Rule 4: NO STANDARD MODEL anywhere in calculator source
   - Rule 5: Public surfaces return ``{'value': X}`` only
   - Rule 6: NO datetime, file writes, classes in ``uqff_pure_calculator.py``
   - Rule 8: RUN ``uqff_fidelity_tests.py`` AFTER EVERY EDIT
   - Rule 9: APPEND to ``SESSION_LOG.md``, never rewrite

2. ``SESSION_LOG.md`` — most recent entry. Tells you the current state.

3. ``NEXT_PRIORITIES.md`` — queued work.

4. ``PRODUCTION_ROADMAP.md`` — Tier-1 / Tier-2 / Tier-3 roadmap.

Contribution workflow
---------------------

1. Fork the repository on GitHub.
2. Create a feature branch.
3. Make your changes (respecting all rules above).
4. Run ``python uqff_fidelity_tests.py`` — must exit 0 (857/857 passing).
5. Append a brief SESSION_LOG.md entry describing your change.
6. Submit a pull request.

Areas where contributions are particularly welcome
---------------------------------------------------

- **Adding regression pins** for the 530 legacy_freeform closures (will
  raise coverage from 46% toward 75%).
- **Authoring whitepapers** that derive new closures from existing
  primitives. Each new closure should be tagged with a PAPER_XXXX source.
- **Code coverage uplift** via fuzz-input tests and synthetic-error tests.
- **Sphinx documentation** — extending the API reference pages.
- **Cross-platform CI** — additional test platforms or Python versions.

Areas where contributions are NOT welcome without prior discussion
-------------------------------------------------------------------

- Modifying locked canonical primitive values (Rule 2 violation).
- Adding SM-named constants, SM theory anchors, or SM-comparison fields
  to the calculator (Rule 4 violation).
- Adding docstrings, comments, or metadata dict-key strings to the
  calculator (Rule 3 violation).
- Restructuring the single-file calculator into multiple modules
  (architectural decision; discuss first).

Code of conduct
---------------

Treat all contributors with respect. Daniel has spent 10+ months on this
project; communications should be honest, direct, and oriented toward
preserving the integrity of UQFF physics. Drift from canonical values or
narrative inflation will be reverted.

License of contributions
------------------------

By submitting a pull request, you agree that your contributions are
licensed under the same AGPL-3.0 + commercial dual license as the rest
of the project, and that you have the right to make such a contribution
under that license.
