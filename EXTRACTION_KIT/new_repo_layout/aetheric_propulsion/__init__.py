"""Aetheric-Propulsion — UQFF assimilation geometry standalone package.

Sibling to `uqff`; bundles the calculator + dispatch + solver bus + 4 geometry +
3 numeric backends as a fully self-contained pip package.

Public API (re-exported from internal modules):
    calculate_analytic_closures(dataset)     -> {'value': X}
    calculate_qcalcgeom_compute_FUBi(dataset)
    calculate_qcalcgeom_compute_FUBii(dataset)
    calculate_qcalcgeom_compute_F_U(dataset)
    calculate_qcalcgeom_solve_habitable_zone(dataset)
    calculate_qcalcgeom_compute_emergent_mass(dataset)
    calculate_3numeric_decomposition(dataset)
    calculate_geometry_decomposition(dataset)
    calculate_overdetermination(dataset)
    solve(observable, geometry='auto', numeric='numerical', decompose=False)
    DISPATCH               # the 116-observable catalog
    domains()              # the 10 dispatch domains
"""
__version__ = "0.1.0"

# Re-export the calculator surface
from .calculator import (
    calculate_analytic_closures,
    calculate_qcalcgeom_compute_FUBi,
    calculate_qcalcgeom_compute_FUBii,
    calculate_qcalcgeom_compute_F_U,
    calculate_qcalcgeom_solve_habitable_zone,
    calculate_qcalcgeom_compute_emergent_mass,
    calculate_3numeric_decomposition,
    calculate_geometry_decomposition,
    calculate_overdetermination,
)

# Re-export the solver bus + dispatch
from .qcalcgeom_solver import solve
from .assimilation_dispatch import DISPATCH, domains

__all__ = [
    "__version__",
    "calculate_analytic_closures",
    "calculate_qcalcgeom_compute_FUBi",
    "calculate_qcalcgeom_compute_FUBii",
    "calculate_qcalcgeom_compute_F_U",
    "calculate_qcalcgeom_solve_habitable_zone",
    "calculate_qcalcgeom_compute_emergent_mass",
    "calculate_3numeric_decomposition",
    "calculate_geometry_decomposition",
    "calculate_overdetermination",
    "solve",
    "DISPATCH",
    "domains",
]
