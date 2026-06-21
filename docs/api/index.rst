API Reference
=============

The UQFF calculator exposes **34 public** ``calculate_*`` surfaces, all
following the stateless contract:

.. code-block:: python

   result: dict = calculate_<name>(dataset: dict) -> {"value": X}

See :doc:`/input_domains` for the full dataset key catalog for each surface.

.. note::

   Per CLAUDE.md Rule 3, the calculator source has **no docstrings**. Sphinx
   autodoc will show signatures and source links only. The functional
   documentation for each surface is in :doc:`/input_domains`.

Module overview
---------------

.. automodule:: uqff_pure_calculator
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: PARADOX_TO_CLOSURE, PARADOX_TO_MILLENNIUM

Public surfaces by category
---------------------------

Primitive computation
~~~~~~~~~~~~~~~~~~~~~

- ``calculate_resonant_adpm`` — ω · cos(π·t_n) · Φ_res
- ``calculate_scm`` — SCm 26-level density with t_n modulation
- ``calculate_f_u_bi`` — Universal Buoyancy F_UBi
- ``calculate_f_u_bi_i`` — F_U_Bi_i 4-layer master integral
- ``calculate_triadic_g`` — w_C·g_comp + w_R·g_res + w_B·g_buoy
- ``calculate_vacuum_ledger`` — 4-term Λ derivation
- ``calculate_analytic_closures`` — symbolic dispatch hub

PAPER-646/1141/1203 canonical
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``calculate_universal_inertial_operator`` — U_i = 2.75e-7 (Sun, t=0)
- ``calculate_nuclear_magic`` — 7 magic numbers + BE/A + deuteron + α
- ``calculate_lenr`` — single-variant LENR (holmlid/parkhomov/pons/mizuno/rossi)
- ``calculate_f_u_zero`` — F_U=0 master equation + r_hz root-find

Phase-2 canonical
~~~~~~~~~~~~~~~~~

- ``calculate_ua_layers`` — UA', UA'', UA''', UA''''
- ``calculate_dpm_grinding`` — 5-step grinding sequence
- ``calculate_caduceus`` — 26 pinch points + π decimals
- ``calculate_shell_orbital`` — Mayer-Jensen shell occupancy
- ``calculate_lenr_full`` — full per-reactor report

BUCKET 0 (loop-closure)
~~~~~~~~~~~~~~~~~~~~~~~

- ``calculate_negative_time_dual_existence``
- ``calculate_si_derivations`` — c at 0.13%, G at 0.08%
- ``calculate_quantum_gravity`` — 26D unification
- ``calculate_black_hole`` — r_min via 26! finite bound
- ``calculate_vds_dvp_bh26`` — VDS + DVP + BH26 spine
- ``calculate_bsd_rank_cohomology`` — BSD 0.30598 (0.005%)

BUCKET A-K (observable suites)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``calculate_paradox`` — 802 paradox/dispatch closure router
- ``calculate_cosmology`` — 56 cosmological observables
- ``calculate_particle_physics`` — 48 particle physics observables
- ``calculate_gw_events`` — 22 GW events
- ``calculate_agn_jet`` — 22 AGN/jet systems
- ``calculate_astrophysics`` — 36 astrophysical systems
- ``calculate_high_energy_astro`` — 10 TeV/PeV sources
- ``calculate_qgp`` — 9 QGP observables
- ``calculate_higgs_precision`` — 13 Higgs observables
- ``calculate_bsm_constraints`` — 17 BSM constraints

Production utilities
~~~~~~~~~~~~~~~~~~~~

- ``calculate_status_report`` — full uncertainty + milestone summary
- ``calculate_whitepaper`` — whitepaper catalog summary
