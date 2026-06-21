Primitive provenance audit
==========================

The full document is ``PROVENANCE_AUDIT.md`` at the repository root. See
also :doc:`primitives` for the 9 truly-independent primitives table.

Provenance grades
-----------------

.. list-table::
   :header-rows: 1

   * - Grade
     - Meaning
   * - **A++**
     - Mathematical necessity (algebra dimension, group order, structural derivation)
   * - **A**
     - First-principles paper with explicit derivation chain
   * - **B+ / B**
     - Paper derives from other UQFF quantities + one measurement anchor
   * - **C+ / C**
     - Paper specifies the value with limited derivation; awaits further analysis

PAPER_1521/1522 LANDMARK reduction
-----------------------------------

In session 2026-06-18, two primitives previously listed as "locked
canonical" were proven derivative:

- ``D_BSFG = D_crit − 2·SO_5 = 26 − 20 = 6 EXACT`` (PAPER_1521)
- ``K_MEX = Φ_5/6 · SO_5 / D_phys = (5/6)·10/4 = 25/12 EXACT`` (PAPER_1522)

Net effect: the "11 locked canonical primitives" count reduced to **9
truly independent + 2 derivative**, strengthening the ΔBIC parameter-
economy comparison against SM+ΛCDM by 2·ln(N_obs).

Resolved open questions
------------------------

**Q3 (resolved 2026-06-18):** Can SSq be derived from F_TRZ + Φ_res?

**Answer:** No. Exhaustive numerical search over all 2-element and
3-element combinations of the other primitives with operators
{+, −, ×, ÷, 1−x, composite} yields NO match within 0.3% of SSq=0.57. The
closest accidental near-matches (4/7 = 0.5714, 7/12 = 0.5833) have no
UQFF structural interpretation. SSq is **truly independent** at the
rational-arithmetic level. PAPER_1154 may derive it from a transcendental
fixed-point relation, but this does not reduce the primitive count.

Open questions (queued for Tier-1B / Tier-2)
---------------------------------------------

- **Q4:** Can β_i = 0.6029 be derived?
- **Q5:** Is N_CH = 9 derivable from D_BSFG + D_phys?
- Per-paper consistency table showing which papers cite which primitive at which precision.
- Numerical sensitivity ∂(closure RMS residual)/∂(primitive value) per primitive.
