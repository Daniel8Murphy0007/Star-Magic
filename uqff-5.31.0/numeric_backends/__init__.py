"""
numeric_backends/__init__.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  uqff_pure_calculator (locked primitives + closures)
Dependencies (external):  sympy>=1.12

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

This package exposes the three independent numeric backends:
    - symbolic   (sympy algebraic identities)
    - numerical  (Python float / sympy.N)
    - discrete   (integer-primitive + 26-state hypergraph per PAPER_1080)

All three implement the same closure set (currently: the 8 Clay Millennium
derivations) and must produce agreeing results within their numeric precision.
The convergence check is the cross-validation harness in
    test_3numeric_millennium_crosscheck.py
"""

from . import symbolic
from . import numerical
from . import discrete

BACKENDS = {
    "symbolic":  symbolic,
    "numerical": numerical,
    "discrete":  discrete,
}

__all__ = ["symbolic", "numerical", "discrete", "BACKENDS"]
