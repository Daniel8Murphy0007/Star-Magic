#!/usr/bin/env python3
"""
ramanujan_cryptography_pi_math.py — PImath Encryption via S₂₆⁽³⁾ Key Generation

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
Implements PImath cryptographic primitives using Ramanujan S₂₆⁽³⁾ spectral
weight as the key expansion seed.

Gap closed:
  - K_PImath = S₂₆⁽³⁾ · π^(n/26) mod 113   (key generation)
  - 26-layer encryption with π-cycle digit-root mapping
  - encrypt()/decrypt() round-trip functions
  - Key expansion from S₂₆⁽³⁾ spectral weight

Physics:
  The 26D compactification layers of UQFF provide 26 independent key
  layers. Each layer n uses π^(n/26) mod 113 as the rotation basis,
  with S₂₆⁽³⁾ as the master seed. The modulus 113 is chosen because
  26! mod 113 = 12 (a UQFF-significant remainder).

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
import hashlib
from typing import Dict, List, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
OMEGA_SCM = 2 * PI * 1.25e12
GAMMA_0   = 2 * PI * 0.1e12
SSQ       = 0.57
BETA_I    = 0.603
H_SCM     = 0.99
F_UBI_RATIO = 0.6
N_LAYERS  = 26
MODULUS   = 113              # 26! mod 113 = 12

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor R_n^{(26,k)}."""
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 0.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, 27):
            sign = (-1) ** (j + 1)
            binom = math.comb(26, j)
            inner += sign * binom * math.factorial(26 - j) / n ** j
        correction += inner / n ** (26 * m)
    return prefactor * (1.0 + correction)


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))


# ── §1  KEY GENERATION ─────────────────────────────────────────────────────

class PIMathKeyGenerator:
    """Generate 26-layer encryption keys from S₂₆⁽³⁾ and π.

    K_PImath(n) = floor(S₂₆⁽³⁾ · π^(n/26) · 10^12) mod 113

    Each of 26 layers provides an independent key byte via the
    transcendental digits of π raised to a fractional 26D power.
    """

    def __init__(self, modulus: int = MODULUS):
        self.modulus = modulus
        self.master_seed = S26_3RD

    def layer_key(self, n: int) -> int:
        """Compute key for layer n (1..26)."""
        raw = self.master_seed * (PI ** (n / 26.0))
        # Extract integer from high-precision float
        scaled = abs(raw) * 1e12
        return int(scaled) % self.modulus

    def full_key(self) -> List[int]:
        """Generate all 26 layer keys."""
        return [self.layer_key(n) for n in range(1, N_LAYERS + 1)]

    def key_hash(self) -> str:
        """SHA-256 hash of the full 26-layer key for verification."""
        key_bytes = bytes(self.full_key())
        return hashlib.sha256(key_bytes).hexdigest()

    def digit_root(self, n: int) -> int:
        """π-cycle digit root: iterated digit sum of π^(n/26) integer part."""
        val = int(abs(PI ** (n / 26.0)) * 1e6) % 10000
        while val >= 10:
            val = sum(int(d) for d in str(val))
        return val

    def digit_root_map(self) -> List[int]:
        """Full 26-layer digit root mapping."""
        return [self.digit_root(n) for n in range(1, N_LAYERS + 1)]

    def compute(self, dataset: dict) -> dict:
        """Full key generation computation."""
        keys = self.full_key()
        dr_map = self.digit_root_map()
        key_sum = sum(keys)

        return {
            "keys": keys,
            "key_sum": key_sum,
            "key_hash_sha256": self.key_hash(),
            "modulus": self.modulus,
            "master_seed": self.master_seed,
            "digit_root_map": dr_map,
            "primary_equations": [
                f"K_PImath(n) = floor(S₂₆⁽³⁾·π^(n/26)·10¹²) mod {self.modulus}",
                f"S₂₆⁽³⁾ = {self.master_seed:.6e}",
                f"Key[1..26] = {keys}",
                f"Digit roots = {dr_map}",
                f"Key sum = {key_sum}",
            ],
        }


# ── §2  26-LAYER ENCRYPT/DECRYPT ──────────────────────────────────────────

class PIMathCipher:
    """26-layer PImath encryption using S₂₆⁽³⁾ key expansion.

    Encryption: For each byte b at position i,
        c_i = (b + K[i mod 26] + DR[i mod 26]) mod 256

    Decryption reverses the operation:
        b = (c_i - K[i mod 26] - DR[i mod 26]) mod 256

    The 26-layer structure mirrors the 26D UQFF compactification,
    with each layer providing independent rotational mixing.
    """

    def __init__(self, modulus: int = MODULUS):
        self.keygen = PIMathKeyGenerator(modulus)
        self.keys = self.keygen.full_key()
        self.dr_map = self.keygen.digit_root_map()

    def encrypt_bytes(self, plaintext: bytes) -> bytes:
        """Encrypt raw bytes through 26-layer rotation."""
        result = bytearray(len(plaintext))
        for i, b in enumerate(plaintext):
            layer = i % N_LAYERS
            result[i] = (b + self.keys[layer] + self.dr_map[layer]) % 256
        return bytes(result)

    def decrypt_bytes(self, ciphertext: bytes) -> bytes:
        """Decrypt raw bytes (inverse of encrypt)."""
        result = bytearray(len(ciphertext))
        for i, b in enumerate(ciphertext):
            layer = i % N_LAYERS
            result[i] = (b - self.keys[layer] - self.dr_map[layer]) % 256
        return bytes(result)

    def encrypt(self, plaintext: str) -> List[int]:
        """Encrypt a string, return list of cipher bytes."""
        raw = plaintext.encode('utf-8')
        enc = self.encrypt_bytes(raw)
        return list(enc)

    def decrypt(self, ciphertext: List[int]) -> str:
        """Decrypt a list of cipher bytes back to string."""
        raw = bytes(ciphertext)
        dec = self.decrypt_bytes(raw)
        return dec.decode('utf-8')

    def compute(self, dataset: dict) -> dict:
        """Encrypt/decrypt round-trip demonstration."""
        message = dataset.get("message", "UQFF Star Magic 26D")
        encrypted = self.encrypt(message)
        decrypted = self.decrypt(encrypted)
        roundtrip_ok = (decrypted == message)

        return {
            "plaintext": message,
            "ciphertext": encrypted,
            "decrypted": decrypted,
            "roundtrip_ok": roundtrip_ok,
            "cipher_len": len(encrypted),
            "keys_used": self.keys,
            "primary_equations": [
                f"Encrypt: c_i = (b + K[i%26] + DR[i%26]) mod 256",
                f"Decrypt: b = (c_i - K[i%26] - DR[i%26]) mod 256",
                f"Round-trip verified: {roundtrip_ok}",
                f"Message length: {len(message)} → cipher length: {len(encrypted)}",
            ],
        }


# ── §3  KEY EXPANSION ANALYSIS ────────────────────────────────────────────

class PIMathKeyAnalysis:
    """Analyze S₂₆⁽³⁾-based key expansion properties.

    Evaluates:
    - Key entropy estimate (bits)
    - Layer independence (correlation between adjacent keys)
    - Modular distribution uniformity
    """

    def compute(self, dataset: dict) -> dict:
        modulus = int(dataset.get("modulus", MODULUS))
        keygen = PIMathKeyGenerator(modulus)
        keys = keygen.full_key()

        # Entropy estimate: H = -Σ p_i log₂(p_i)
        counts = [0] * modulus
        for k in keys:
            counts[k] += 1
        probs = [c / N_LAYERS for c in counts if c > 0]
        entropy = -sum(p * math.log2(p) for p in probs)

        # Adjacent key correlation
        diffs = [abs(keys[i+1] - keys[i]) for i in range(len(keys) - 1)]
        mean_diff = sum(diffs) / len(diffs) if diffs else 0

        # Unique keys
        unique = len(set(keys))

        return {
            "keys": keys,
            "modulus": modulus,
            "entropy_bits": entropy,
            "unique_keys": unique,
            "mean_adjacent_diff": mean_diff,
            "max_key": max(keys),
            "min_key": min(keys),
            "S26_3RD": S26_3RD,
            "primary_equations": [
                f"K_PImath(n) = floor(S₂₆⁽³⁾·π^(n/26)·10¹²) mod {modulus}",
                f"Entropy = {entropy:.4f} bits",
                f"Unique keys: {unique}/{N_LAYERS}",
                f"Mean adjacent diff: {mean_diff:.2f}",
            ],
        }


# ── §4  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("PImath Cryptography — S₂₆⁽³⁾ Key Generation")
    print("=" * 72)

    # Key generation
    keygen = PIMathKeyGenerator()
    kg_result = keygen.compute({})
    print(f"\nS₂₆⁽³⁾ = {kg_result['master_seed']:.6e}")
    print(f"Keys [1..26]: {kg_result['keys']}")
    print(f"Digit roots:  {kg_result['digit_root_map']}")
    print(f"SHA-256: {kg_result['key_hash_sha256']}")

    # Encrypt/decrypt
    cipher = PIMathCipher()
    msg = "UQFF Star Magic 26D"
    cr = cipher.compute({"message": msg})
    print(f"\nPlaintext:  '{cr['plaintext']}'")
    print(f"Ciphertext: {cr['ciphertext']}")
    print(f"Decrypted:  '{cr['decrypted']}'")
    print(f"Round-trip: {cr['roundtrip_ok']}")

    # Key analysis
    analysis = PIMathKeyAnalysis()
    ar = analysis.compute({})
    print(f"\nEntropy: {ar['entropy_bits']:.4f} bits")
    print(f"Unique keys: {ar['unique_keys']}/{N_LAYERS}")

    print("\n✓ All PImath cryptography calculations complete")
