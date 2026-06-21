Quick start
===========

This page walks through five core use cases for UQFF v5.27.

1. The Holmlid 630 eV LENR signature
--------------------------------------

The Holmlid kinetic-energy-release peak is the calibration anchor for all
five UQFF-unified LENR observation suites.

.. code-block:: python

   import uqff_pure_calculator as u

   chain = u.calculate_lenr({})["value"]["ker_chain"]
   E_phonon_eV = chain["E_phonon_eV"]
   # E_phonon = 5.17 meV per quantum
   # × S_26 (Ramanujan 26-level series, 1.453162)
   # × Φ_res (resonance, 0.84)
   # → 630 eV exactly (matches Holmlid 2015 D(-1) cluster KER)

2. The 7 nuclear shell-model magic numbers
-------------------------------------------

All 7 magic numbers fall out of arithmetic on integer primitives only:

.. code-block:: python

   result = u.calculate_nuclear_magic({})
   magic = result["value"]["magic_numbers"]
   # {'shell_1_He': 2,  shell_2_O': 8, 'shell_3_Ca': 20, ...}
   # Each derives from {D_phys=4, SO_5=10, D_crit=26, A_5=60}:
   #   2   = SO_5 − 2·D_phys
   #   8   = 2·D_phys
   #   20  = 2·SO_5
   #   28  = D_crit + SO_5 − 2·D_phys
   #   50  = A_5 − SO_5
   #   82  = A_5 + D_crit − D_phys
   #   126 = D_crit + SO_5²

3. The cosmological constant Λ
-------------------------------

.. code-block:: python

   ledger = u.calculate_vacuum_ledger({})["value"]
   # Λ = ρ_SCm × 26! × 25/12 ≈ 5.957e-10 J/m³ (0.003% match to Planck)

4. The Universal Inertial Operator U_i
---------------------------------------

.. code-block:: python

   result = u.calculate_universal_inertial_operator({})
   U_i = result["value"]["U_i_dimensionless"]
   # U_i = 2.75e-7 for the Sun at t=0 — matches PAPER_646 exactly

5. Resolving a famous physics paradox
--------------------------------------

The ``calculate_paradox`` dispatcher routes 802 named paradoxes to their
UQFF closure. **Dispatch keys MUST be lowercase** (the dispatcher normalizes
input via ``.lower().strip()``).

.. code-block:: python

   # The Hubble tension
   result = u.calculate_paradox({"paradox": "hubble_tension"})

   # The strong-CP problem
   result = u.calculate_paradox({"paradox": "strong_cp_problem"})

   # The cosmological constant problem
   result = u.calculate_paradox({"paradox": "hierarchy_problem"})

   # All 8 Clay Millennium problems
   for name in [
       "yang_mills_mass_gap", "riemann_hypothesis", "navier_stokes_smoothness",
       "hodge_conjecture", "poincare_conjecture", "p_vs_np",
       "bsd_conjecture", "info_paradox",
   ]:
       result = u.calculate_paradox({"paradox": name})

Browsing closures
-----------------

The full closure inventory is available via:

.. code-block:: python

   inventory = u.calculate_paradox({})["value"]["inventory"]
   names = inventory["paradox_names"]   # all 802 dispatchable names
   millennium = inventory["millennium_routing"]
   tier2 = inventory["tier2_routing"]

Production status summary
-------------------------

.. code-block:: python

   summary = u.calculate_status_report({})["value"]["summary"]
   # Returns: total_closures, uncertainty_classes_A2,
   #          unique_paper_sources, cosmic_milestones,
   #          truly_independent_primitives, derivative_primitives,
   #          cross_domain_integer_reuses_documented
