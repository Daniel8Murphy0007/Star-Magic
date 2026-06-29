Code coverage
==============

The full document is ``COVERAGE.md`` at the repository root. Headline:

Current state
-------------

- **Tool:** ``coverage.py`` 7.14.2
- **Measurement:** 45.68% (8,344 / 18,265 statements covered)
- **Fidelity gate:** 857 / 857 tests passing during measurement

Reproducing
-----------

.. code-block:: bash

   pip install coverage
   python -m coverage run --source=uqff_pure_calculator uqff_fidelity_tests.py
   python -m coverage report --include='*uqff_pure_calculator*'

.. note::

   The ``uqff_fidelity_tests.py`` file may have trailing null bytes from
   previous repair operations. Strip them before measurement:

   .. code-block:: bash

      python -c "import pathlib; p=pathlib.Path('uqff_fidelity_tests.py'); d=p.read_bytes(); p.write_bytes(d.replace(b'\x00',b''))"

Why 46%?
--------

The calculator has 18,265 statements, of which ~5,000 are the
function-body lines of the 616 unique dispatch closures. Each closure is
~8.9 executable lines and is exercised when ``calculate_paradox`` is
called with its name (covered by the new block #58 sweep added in
session 2026-06-18).

The uncovered 54% breaks down as:

- ~5,500 statements in the 530 legacy_freeform closure bodies (each
  needs a regression pin to fully measure)
- ~2,000 statements in unused dispatch-key aliases (multiple keys
  pointing to the same function)
- ~1,000 statements in defensive try/except branches and isinstance guards
- ~1,000 statements in helper functions only called by unwired reactor
  variants
- ~400 statements in unreachable defensive code paths

Path to higher coverage
-----------------------

.. list-table::
   :header-rows: 1

   * - Phase
     - Target
     - Method
   * - Tier-1B
     - 75%
     - Add regression pins for the 530 legacy_freeform closures
   * - Tier-2 v1.0
     - 85%
     - Synthetic-error-input fuzz tests + per-reactor LENR exercise
   * - Tier-2 v1.1
     - 90%
     - Mark unreachable branches with ``# pragma: no cover``
   * - Tier-3
     - 95%+
     - Property-based testing (Hypothesis) over numerical inputs

Comparable projects
-------------------

For context:

- SciPy: ~85% (mature, decades of work)
- NumPy: ~90% (dedicated coverage CI)
- Astropy: ~85%
- pandas: ~92%
- **UQFF (today): 46%** (first measurement, ~6 weeks of intensive wiring)
