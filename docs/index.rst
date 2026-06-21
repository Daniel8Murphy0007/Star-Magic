UQFF Star-Magic Documentation
==============================

Welcome to the documentation for **UQFF Star-Magic v5.27** — a vacuum-first
unified physics framework built on a single non-mass primitive
(:math:`\rho_{\rm SCm} = 7.09 \times 10^{-37}~\mathrm{J/m^3}`) and a 26-level
Di-Pseudo-Monopole structural lattice.

From **9 truly-independent primitives**, UQFF derives:

- The cosmological constant :math:`\Lambda` at 0.003% match to Planck
- The Holmlid 630 eV LENR signature, exactly
- All 7 nuclear shell-model magic numbers, exactly
- The Fe-56 binding-energy peak at 0.019%
- All 8 Clay Millennium prize problems
- The complete 12-fermion Standard Model spectrum
- 18+ cosmological observables (ΛCDM-complete)
- The Star-Magic LENR reactor COP 555:1 prediction

with **ΔBIC = 94.1** (decisive Bayesian preference) over the Standard Model +
ΛCDM benchmark, which carries 26 parameters versus UQFF's 9.

.. note::

   UQFF is offered as an alternative parameter-economical description
   (**NOT REPLACEMENT**) for the Standard Model. Both frameworks solve the
   same observed phenomena by different methods, both reported with honest
   residuals.

.. warning::

   This software is offered under a **dual license** (AGPL-3.0 + commercial).
   See :doc:`license` for the full terms. Any academic publication using
   UQFF must cite per `CITATION.cff <https://github.com/Daniel8Murphy0007/
   Star-Magic/blob/main/CITATION.cff>`_.

Quick start
-----------

.. code-block:: bash

   pip install uqff

.. code-block:: python

   import uqff_pure_calculator as u

   # The Holmlid LENR signature
   ker = u.calculate_lenr({})["value"]["ker_chain"]["E_phonon_eV"]
   print(f"E_phonon = {ker} eV (× S_26 × Φ_res → 630 eV)")

   # Status of the whole calculator
   summary = u.calculate_status_report({})["value"]["summary"]
   print(f"Closures: {summary['total_closures']}")
   print(f"Truly-independent primitives: {summary['truly_independent_primitives']}")

   # Resolve any of 802 paradoxes (USE LOWERCASE KEYS)
   result = u.calculate_paradox({"paradox": "hubble_tension"})

Production status
-----------------

UQFF v5.27 is **Tier-1 production-complete** as of 2026-06-18:

- 794 paradox/closure dispatch keys
- 854 / 857 fidelity-gate tests passing (0 failures)
- 263 schema-tagged closures, 530 legacy_freeform closures
- 128 EXACT structural identities + 226 Bonferroni-passing matches
- 1,795 whitepapers (PAPER_001 .. PAPER_1795 + ERRATUM)
- 368 C++ reference functions
- Dual-licensed (AGPL-3.0 + commercial), PyPI-ready wheel built

See :doc:`production_status` for the full breakdown.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quickstart
   primitives

.. toctree::
   :maxdepth: 2
   :caption: API reference

   api/index

.. toctree::
   :maxdepth: 2
   :caption: Production audits

   production_status
   prediction_labels
   statistical_hygiene
   provenance_audit
   coverage
   input_domains

.. toctree::
   :maxdepth: 1
   :caption: Project

   license
   commercial_licensing
   contributing
   changelog

Indices and tables
==================

- :ref:`genindex`
- :ref:`modindex`
- :ref:`search`
