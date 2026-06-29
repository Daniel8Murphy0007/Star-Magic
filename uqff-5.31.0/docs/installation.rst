Installation
============

Requirements
------------

- Python 3.10 or newer (tested on 3.10, 3.11, 3.12, 3.13)
- No runtime dependencies (UQFF uses only Python standard library)
- Optional: ``coverage`` for code coverage measurement
- Optional: ``sphinx`` + ``sphinx-rtd-theme`` for building these docs

From PyPI
---------

.. code-block:: bash

   pip install uqff

This installs the ``uqff_pure_calculator`` module as a single-file
distribution (~2.6 MB).

From source
-----------

.. code-block:: bash

   git clone https://github.com/Daniel8Murphy0007/Star-Magic.git
   cd Star-Magic
   pip install .

For development installation:

.. code-block:: bash

   pip install -e ".[test,docs]"

Verifying the install
---------------------

After install, run the fidelity gate:

.. code-block:: bash

   python -c "import uqff_pure_calculator as u; print('UQFF', u.__name__ ,'loaded')"
   python uqff_fidelity_tests.py

Expected output: ``TOTAL: 857 passed, 0 failed``.

License
-------

UQFF is distributed under a **DUAL LICENSE** model. You must choose one:

- **AGPL-3.0** (free for academic/research/non-commercial use)
- **Commercial license** (for proprietary deployments; contact
  ``daniel.murphy00@enrgyone.com``)

See the :doc:`license` and :doc:`commercial_licensing` pages for full terms.
