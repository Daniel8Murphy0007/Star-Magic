Input domain documentation
===========================

The complete per-surface dataset-key catalog is in ``INPUT_DOMAINS.md`` at
the repository root.

Universal contract
------------------

Every public surface follows the stateless contract:

.. code-block:: python

   result: dict = calculate_<name>(dataset: dict) -> {"value": X}

- ``dataset`` is always a Python ``dict``. May be empty ``{}`` to get
  default-parameter behavior.
- Return is always ``{"value": X}`` where X may be a number, dict, or
  list-of-observables.
- No exceptions raised for malformed input; instead ``{"value": None}``
  is returned.
- All surfaces stateless and side-effect-free (CLAUDE.md Rule 6).
- All numerical values use SI units unless explicitly noted.

Common pitfalls
---------------

1. **Mixed-case dispatch keys silently fail** for ``calculate_paradox``.
   The dispatcher lowercases input but does NOT raise on miss. Always
   use lowercase keys: ``"hubble_tension"``, NOT ``"Hubble_Tension"``.

2. **Empty dataset returns the canonical-primitive default**, not None.
   To get a "no data" signal, set explicit ``target=None`` or
   ``paradox=None``.

3. **Time ``t_n`` is dimensionless in [0,1]**, not in seconds. Multiply
   by the canonical period if you need a physical-time interpretation.

4. **All masses in kg, distances in m, frequencies in Hz**. There is no
   internal unit-conversion layer.

5. **Returned ``'anchor'`` value is the literal measured observable**,
   not a SM-theory prediction. Per Rule 4, the calculator never compares
   to SM-derived anchors; only to direct experimental measurements.

For per-surface detail
-----------------------

See ``INPUT_DOMAINS.md`` for the full dataset key tables for all 34
surfaces.
