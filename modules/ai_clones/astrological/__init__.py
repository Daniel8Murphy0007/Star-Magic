"""
Astrological AI Clone Modules
=============================

High-energy astrophysical calculators using UQFF physics.

Modules:
- magnetar_cycles: Magnetar burst timing analysis
- pulsar_harmonics: Pulsar timing and gravitational wave correlation
- gw_event_predictor: Gravitational wave event prediction
- solar_flare_index: Solar activity forecasting

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
"""

from pathlib import Path
import sys

# Ensure parent modules are importable
_module_root = Path(__file__).parent.parent.parent
if str(_module_root) not in sys.path:
    sys.path.insert(0, str(_module_root))

# Available modules
__all__ = [
    'MagnetarCyclesCalculator',
]

# Lazy imports
def __getattr__(name):
    if name == 'MagnetarCyclesCalculator':
        from modules.ai_clones.astrological.magnetar_cycles import MagnetarCyclesCalculator
        return MagnetarCyclesCalculator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
