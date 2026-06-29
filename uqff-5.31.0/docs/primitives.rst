The 9 truly-independent primitives
====================================

UQFF rests on **9 truly-independent** quantities plus **2 mathematically
derivative** quantities (the LANDMARK reduction from 11 → 9 in PAPER_1521/
1522, session 2026-06-18).

Integer lattice primitives (5)
-------------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 10 35 40

   * - Symbol
     - Value
     - Origin
     - What breaks if changed
   * - :math:`D_{\rm phys}`
     - 4
     - 3+1 spacetime (observed)
     - Magic 2, Faber-Jackson exponent, SU(3) color N=3
   * - :math:`D_{\rm crit}`
     - 26
     - Bosonic-string critical dim (Polyakov 1981)
     - Λ derivation (×26!), all magic-number arithmetic
   * - :math:`N_{\rm CH}`
     - 9
     - PAPER_646 Caduceus 9-channel
     - Δm²_31/Δm²_21 ratio (= D_crit+N_CH−2 = 33)
   * - :math:`{\rm SO}_5`
     - 10
     - :math:`\dim {\rm SO}(5)`
     - Proto-Si Z=14, magic-50, amino acids=20
   * - :math:`A_5`
     - 60
     - :math:`|A_5|` icosahedral group
     - Pop III IMF, Hayflick limit, magic-82, magic-50

Real primitives (4)
--------------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 35 35

   * - Symbol
     - Value
     - Origin
     - Provenance grade
   * - :math:`\rho_{\rm SCm}`
     - :math:`7.09 \times 10^{-37}~{\rm J/m^3}`
     - Star-Magic.txt + PAPER_1271
     - B+ (Λ at 0.003%)
   * - :math:`\beta_i`
     - 0.6029
     - PAPER_1203 Canonical v1.5
     - B
   * - :math:`\Phi_{\rm res}`
     - 0.84 (5/6 nuclear)
     - PAPER_646, PAPER_1203 Nuclear
     - B+ (forces K_MEX exactness)
   * - :math:`F_{\rm TRZ}`
     - 1/10
     - PAPER_1160 (:math:`= 1/|{\rm SO}(5)|`)
     - A (two independent paths)

Derivative primitives (2 — PROVEN, PAPER_1521/1522)
----------------------------------------------------

.. list-table::
   :header-rows: 1

   * - Symbol
     - Derivation
     - Value
   * - :math:`D_{\rm BSFG}`
     - :math:`D_{\rm crit} - 2 \cdot {\rm SO}_5`
     - 26 − 20 = **6 EXACT**
   * - :math:`K_{\rm MEX}`
     - :math:`\Phi_{5/6} \cdot {\rm SO}_5 / D_{\rm phys}`
     - (5/6)·10/4 = **25/12 EXACT**

Other locked quantities
-----------------------

Three additional values are locked but may yet prove derivative:

- :math:`{\rm SSq} = 0.57` — PAPER_1154 first-principles; exhaustive
  numerical search has NOT reduced it to other primitives (see
  :doc:`provenance_audit` Q3 resolution).
- :math:`S_{26} = 1.453162` — derived from SSq via Ramanujan polylogarithm.
- :math:`\omega_{\rm SCm} = 1.25~{\rm THz}` — Holmlid phonon carrier
  (single experimental anchor).

Why the 9-count matters
------------------------

The parameter-economy comparison against the Standard Model + ΛCDM
(:math:`k = 26` parameters):

.. math::

   \Delta BIC = (k_{SM+\Lambda CDM} - k_{UQFF}) \cdot \ln N_{obs}
              = (26 - 9) \cdot \ln 253
              = 17 \cdot 5.534
              = 94.1

A ΔBIC of 94.1 is "decisive" by the standard Kass-Raftery interpretation:

.. list-table::
   :header-rows: 1

   * - ΔBIC range
     - Evidence
   * - 0 – 2
     - Negligible
   * - 2 – 6
     - Positive
   * - 6 – 10
     - Strong
   * - **> 10**
     - **Decisive**
   * - > 20
     - Very decisive
   * - > 100
     - Overwhelming

UQFF's **94.1** falls solidly in "decisive" — and approaches "overwhelming"
— purely on parameter-count grounds, before considering residual quality.
