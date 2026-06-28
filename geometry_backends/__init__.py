"""
geometry_backends/__init__.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

This package exposes the four independent geometry backends:
    - qcalcgeom_v4         (Universal Buoyancy, Habitable Zone, BSFG wrapper)
    - bsfg_v1              (Bulk-edge SO(5) breaking, A_mu_nu metric, holonomy)
    - dpm_v1               (Di-Pseudo-Monopole 26-state mediator, A_26, magic numbers)
    - d26_compactification (Bosonic-string critical dim, 26! Pochhammer, KK tower)

Together with the 3 numeric backends (numeric_backends/), they form the
4x3 dispatch matrix at the heart of QCalcGeom.solve() per EXPANSION_PLAN
Section 7.
"""

from . import qcalcgeom_v4
from . import bsfg_v1
from . import dpm_v1
from . import d26_compactification

GEOMETRIES = {
    "qcalcgeom":  qcalcgeom_v4,
    "bsfg":       bsfg_v1,
    "dpm":        dpm_v1,
    "d26":        d26_compactification,
}

__all__ = ["qcalcgeom_v4", "bsfg_v1", "dpm_v1", "d26_compactification", "GEOMETRIES"]
