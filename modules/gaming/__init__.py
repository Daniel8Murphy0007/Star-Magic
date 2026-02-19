"""
Gaming Modules
==============

Interactive physics-based gaming modules using UQFF framework.

Categories:
- wormhole_traveler: Morris-Thorne wormhole navigation
- cosmic_mapping: Galaxy visualization and exploration
- physics_paradigm: UQFF vs GR physics comparison

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
"""

from pathlib import Path
import sys

# Ensure parent modules are importable
_module_root = Path(__file__).parent.parent
if str(_module_root) not in sys.path:
    sys.path.insert(0, str(_module_root))

# Available modules
__all__ = [
    'WormholeTrajectoryCalculator',
    'CosmicMapper',
]

# Lazy imports
def __getattr__(name):
    if name == 'WormholeTrajectoryCalculator':
        from modules.gaming.wormhole_traveler.trajectories import WormholeTrajectoryCalculator
        return WormholeTrajectoryCalculator
    if name == 'CosmicMapper':
        from modules.gaming.cosmic_mapping.galaxy_render import CosmicMapper
        return CosmicMapper
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
