# Aetheric-Propulsion — UQFF Assimilation Geometry Standalone

[![PyPI version](https://img.shields.io/pypi/v/aetheric-propulsion.svg)](https://pypi.org/project/aetheric-propulsion/)
[![License: AGPL-3.0 + Commercial](https://img.shields.io/badge/License-AGPL--3.0%20%2B%20Commercial-blue.svg)](LICENSE-AGPL-3.0.txt)

**Standalone bundle of the UQFF Phase E/F/G assimilation geometry.**

This package is a sibling to [`uqff`](https://github.com/Daniel8Murphy0007/Star-Magic).
It bundles the calculator, 4 geometry backends, 3 numeric backends, the 116-observable
assimilation dispatch catalog, and the BAO + Cabibbo dual-closure infrastructure as
a fully self-contained pip package. **No runtime dependency on `uqff`.**

## What this is for

`uqff` is the canonical academic-access framework (free, AGPL-3.0). `aetheric-propulsion`
is the staging ground for commercial-tier extensions: hardware reactor control,
proprietary API endpoints, custom calibration suites. During the peer-review phase
(three universities + 8 NASA-Roses grant panels), both repositories ship the same
core science; only the long-term extension paths differ.

## Quick start

```bash
pip install aetheric-propulsion
python -c "
from aetheric_propulsion import calculate_analytic_closures
r = calculate_analytic_closures({'qcalcgeom_solve': {'observable': 'alpha_inverse'}})
print(r)
# -> {'value': 137.0}
"
```

## What's included (from Star-Magic v5.30.0 extraction)

| Module | Contents |
|---|---|
| `aetheric_propulsion.calculator` | Full UQFF pure calculator (42 public `calculate_*` surfaces) |
| `aetheric_propulsion.assimilation_dispatch` | 116-observable dispatch catalog across 10 domains |
| `aetheric_propulsion.qcalcgeom_solver` | 4 x 3 dispatch matrix solver bus |
| `aetheric_propulsion.geometry_backends` | qcalcgeom_v4, bsfg_v1, dpm_v1, d26_compactification |
| `aetheric_propulsion.numeric_backends` | symbolic, numerical, discrete |
| `aetheric_propulsion.provenance_recorder` | per-closure provenance chain assembly |

Plus bundled documentation: `OVERDETERMINATION_MAP.csv` (long format),
`OVERDETERMINATION_WIDE.csv`, `OVERDETERMINATION_MAP.md`, `ASSIMILATION_GEOMETRY_ATLAS.md`.

## Multi-path corroboration — same as Star-Magic

Both BAO and Cabibbo demonstrate the framework's evidence framework: two structurally-
independent UQFF closures sharing only one or two primitives converging on the same
observable at varying precision. Joint probability of random agreement at <0.03% is
below 10^-6.

Reference: `PAPER_1156` Appendix A (in the Star-Magic whitepapers/ directory).

## License

AGPL-3.0-or-later for academic, research, personal, and non-commercial use.
Commercial license available separately — see `COMMERCIAL.md`.

## Provenance

This package is mechanically extracted from
[Star-Magic](https://github.com/Daniel8Murphy0007/Star-Magic) via the
`EXTRACTION_KIT/_step7_migrate_to_aetheric_propulsion.py` script. The Star-Magic
SESSION_LOG.md provides the round-by-round audit trail for every closure.
