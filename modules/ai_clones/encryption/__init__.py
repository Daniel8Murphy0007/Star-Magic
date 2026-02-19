"""
Encryption AI Clone Modules
===========================

UQFF-based encryption and cryptographic calculators.

Modules:
- aes_256_uqff: AES-256 with UQFF field-derived keys
- rsa_8192_cosmic: RSA-8192 with cosmic parameter seeding
- lattice_crypto: Lattice-based post-quantum cryptography

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
    'AES256UQFFCalculator',
]

# Lazy imports
def __getattr__(name):
    if name == 'AES256UQFFCalculator':
        from modules.ai_clones.encryption.aes_256_uqff import AES256UQFFCalculator
        return AES256UQFFCalculator
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
