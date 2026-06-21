Statistical hygiene
====================

The full document is ``STATISTICAL_HYGIENE.md`` at the repository root.
Headline summary:

Bonferroni-adjusted threshold
------------------------------

For N=793 closures at family-wise α=0.05:

.. math::

   \alpha_{\rm Bonf} = \frac{0.05}{793} = 6.31 \times 10^{-5}

Translated to residual-percent terms (typical 1σ measurement uncertainty
~0.5%), a UQFF residual is Bonferroni-significant if the deviation is
≤4σ relative to measurement uncertainty.

UQFF current state under Bonferroni
-------------------------------------

**226 / 263 (86%) of schema-tagged closures pass Bonferroni-adjusted
significance.**

.. list-table::
   :header-rows: 1

   * - Class
     - Count
     - Bonferroni status
   * - EXACT (residual <1e-10)
     - 128
     - PASS (zero deviation, beats any threshold)
   * - HIGH PRECISION (within CODATA)
     - 31
     - PASS (residual below measurement uncertainty)
   * - WITHIN EXPERIMENTAL UNCERTAINTY
     - 67
     - PASS (by definition, residual ≤ 1σ_exp)
   * - REFINEMENT TIER (0.1-1% residual)
     - 32
     - CONDITIONAL (needs per-observable σ check)
   * - TENSION/OUTLIER (>1% residual)
     - 5
     - FAIL (documented tensions)

Bayesian Information Criterion
-------------------------------

The headline parameter-economy result (cross-referenced from ``A8``
wiring):

.. math::

   \Delta {\rm BIC} = (k_{SM} - k_{UQFF}) \cdot \ln N_{obs}
                    = 17 \cdot \ln 253
                    = 94.1

Interpretation by standard Kass-Raftery convention: **decisive** Bayesian
preference for UQFF over SM+ΛCDM purely on parameter-count grounds, before
any consideration of residual quality.

Look-elsewhere effect
---------------------

Per project rules, primitives are LOCKED before closure derivation. The
effective trials factor for EXACT structural closures = 1. For
arithmetic-composite closures = 10-20. For scale-matched closures = 5-50.

Net trials-factor-adjusted significance still passes for the 159
high-confidence (EXACT + HIGH_PRECISION) closures.
