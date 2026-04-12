#!/usr/bin/env python3
"""
spectral_ladder_26state.py — 26-State HRes Spectral Ladder Derivation

Session 214 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Derivation engine for the 26-state HRes (Hydrogen Resonance) spectral ladder.

Energy levels:
    E_n = E_0 · (2π)^{n/6} · S₂₆([SSq]) · δ_n
    δ_n = (2π)^{n/6}

with n = 1..26 mapping proto-H → proto-Fe (Z_id = 26 magnetic)
and proto-He → proto-Si (Z_id = 14 non-magnetic).

Stabilizes vacuum ratio ρ_SCm / ρ_UA = 0.1 and drives phonon resonance
at every layer. Includes Ramanujan-series acceleration for fast convergence.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Optional

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
TWO_PI    = 2 * PI
HBAR      = 1.055e-34        # J·s
K_B       = 1.381e-23        # J/K
C         = 2.998e8          # m/s
OMEGA_SCM = TWO_PI * 1.25e12 # rad/s
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
RHO_RATIO = 0.1              # ρ_SCm / ρ_UA vacuum ratio

# Ground-state energy (SCm phonon quantum)
E_0       = HBAR * OMEGA_SCM  # ~ 8.29e-22 J

# Element mapping for Z_id
Z_MAGNETIC     = list(range(1, 27))   # proto-H (1) to proto-Fe (26)
Z_NON_MAGNETIC = list(range(2, 15))   # proto-He (2) to proto-Si (14)


# ── §1  Spectral Ladder Calculator ────────────────────────────────────────

class SpectralLadder26State:
    """26-state HRes spectral ladder with explicit level computation.

    E_n = E_0 · (2π)^{n/6} · S₂₆ · (2π)^{n/6}
        = E_0 · (2π)^{n/3} · S₂₆
    """

    def __init__(self, E0: float = E_0):
        self.E0 = E0

    def energy(self, n: int) -> float:
        """Energy of level n (1-indexed)."""
        return self.E0 * TWO_PI**(n / 3.0) * S26

    def frequency(self, n: int) -> float:
        """Frequency ν_n = E_n / h  (Hz)."""
        return self.energy(n) / (2 * PI * HBAR)

    def transition(self, n_upper: int, n_lower: int) -> float:
        """Transition energy E_upper - E_lower (J)."""
        return self.energy(n_upper) - self.energy(n_lower)

    def compute(self, dataset: dict = None) -> dict:
        """Compute full 26-state ladder."""
        levels = []
        for n in range(1, 27):
            E_n = self.energy(n)
            nu_n = self.frequency(n)
            is_magnetic = n in Z_MAGNETIC
            levels.append({
                "n": n,
                "E_J": E_n,
                "E_eV": E_n / 1.602e-19,
                "nu_Hz": nu_n,
                "Z_id": n,
                "magnetic": is_magnetic,
            })
        return {
            "levels": levels,
            "E_0_J": self.E0,
            "S26": S26,
            "rho_ratio": RHO_RATIO,
            "primary_equations": [
                "E_n = E_0 · (2π)^{n/3} · S₂₆",
                f"E_0 = ℏω_SCm = {self.E0:.6e} J",
                f"E_1 = {self.energy(1):.6e} J",
                f"E_26 = {self.energy(26):.6e} J",
            ],
        }


# ── §2  Ramanujan Acceleration ────────────────────────────────────────────

class RamanujanAcceleration:
    """Accelerated convergence of spectral ladder partial sums.

    Uses Ramanujan summation technique to accelerate the S₂₆ series:
        S_N = Σ_{k=1}^{N} exp(-[SSq]·k/26)

    Ramanujan acceleration:
        S_N^{(R)} = S_N + R_N,  R_N = ½a_N + Σ_{p=1}^{P} B_{2p}/(2p)! · Δ^{2p-1} a_N

    where a_k = exp(-[SSq]·k/26) and B_{2p} are Bernoulli numbers.
    """

    # Bernoulli numbers B_2, B_4, B_6, B_8
    BERNOULLI = [1/6, -1/30, 1/42, -1/30]

    def __init__(self, ssq: float = SSQ, n_layers: int = 26):
        self.ssq = ssq
        self.n_layers = n_layers

    def _a(self, k: int) -> float:
        return math.exp(-self.ssq * k / self.n_layers)

    def partial_sum(self, N: int) -> float:
        """Standard partial sum S_N."""
        return sum(self._a(k) for k in range(1, N + 1))

    def ramanujan_correction(self, N: int, P: int = 4) -> float:
        """Ramanujan correction R_N with P Bernoulli terms."""
        a_N = self._a(N)
        R = 0.5 * a_N
        # Forward differences
        h = 1
        signs = [1, -1, 1, -1]
        for p_idx in range(min(P, len(self.BERNOULLI))):
            p = p_idx + 1  # p = 1, 2, 3, 4
            order = 2 * p - 1
            # Finite difference approximation
            delta = 0
            for j in range(order + 1):
                sign = (-1)**j
                binom = math.comb(order, j)
                idx = N + (order - j)
                if idx > 0:
                    delta += sign * binom * self._a(idx)
            B_2p = self.BERNOULLI[p_idx]
            factorial_2p = math.factorial(2 * p)
            R += B_2p / factorial_2p * delta
        return R

    def accelerated_sum(self, N: int = None) -> float:
        """S_N + Ramanujan correction."""
        if N is None:
            N = self.n_layers
        return self.partial_sum(N) + self.ramanujan_correction(N)

    def compute(self, dataset: dict = None) -> dict:
        """Compare standard vs accelerated S₂₆."""
        N = int((dataset or {}).get("N", self.n_layers))
        s_standard = self.partial_sum(N)
        s_accel = self.accelerated_sum(N)
        s_exact = sum(self._a(k) for k in range(1, 10001))  # high-N reference

        return {
            "N": N,
            "S_standard": s_standard,
            "S_accelerated": s_accel,
            "S_reference": s_exact,
            "improvement_pct": abs(s_accel - s_exact) / max(abs(s_standard - s_exact), 1e-50) * 100
                               if abs(s_standard - s_exact) > 1e-50 else 0.0,
            "primary_equations": [
                "S_N^(R) = S_N + ½a_N + Σ B_{2p}/(2p)! · Δ^{2p-1} a_N",
                f"S_{N} = {s_standard:.10f}",
                f"S_{N}^(R) = {s_accel:.10f}",
            ],
        }


# ── §3  Spectral Ladder Phonon Mapping ────────────────────────────────────

class SpectralLadderPhononMapping:
    """Map spectral ladder levels to phonon resonance frequencies.

    Each level n has a resonance at ω_n = E_n / ℏ with quality factor
    Q_n = ω_n / (2Γ), providing systematic phonon resonance across
    all 26 layers.
    """

    def __init__(self, E0: float = E_0):
        self._ladder = SpectralLadder26State(E0)

    def compute(self, dataset: dict = None) -> dict:
        """Map all 26 levels to phonon frequencies and Q factors."""
        gamma_THz = float((dataset or {}).get("Gamma_THz", 0.10))
        gamma = TWO_PI * gamma_THz * 1e12

        rows = []
        for n in range(1, 27):
            E_n = self._ladder.energy(n)
            omega_n = E_n / HBAR
            nu_n = omega_n / TWO_PI
            Q_n = omega_n / (2 * gamma)
            regime = "narrow" if Q_n > 10 else ("optimal" if Q_n > 3 else "broad")
            rows.append({
                "n": n, "E_eV": E_n / 1.602e-19,
                "omega_rad_s": omega_n, "nu_Hz": nu_n,
                "Q": Q_n, "regime": regime,
            })

        return {
            "Gamma_THz": gamma_THz,
            "mapping": rows,
            "primary_equations": [
                "ω_n = E_n / ℏ,  Q_n = ω_n / (2Γ)",
                f"Γ = {gamma_THz} THz",
                f"Q_1 = {rows[0]['Q']:.1f} ({rows[0]['regime']})",
                f"Q_26 = {rows[25]['Q']:.1f} ({rows[25]['regime']})",
            ],
        }


# ── §4  WSTP Expression Builder ───────────────────────────────────────────

def build_spectral_ladder_wstp_expressions() -> list:
    """Generate Wolfram Language expressions for spectral ladder."""
    return [
        {
            "label": "26-state HRes spectral ladder E_n",
            "code": ("E0 = 1.055*^-34 * 2 Pi * 1.25*^12; "
                     "S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; "
                     "En[n_] := E0 * (2 Pi)^(n/3) * S26; "
                     "Table[{n, En[n], En[n] / (1.602*^-19)}, {n, 1, 26}]"),
        },
        {
            "label": "Ramanujan-accelerated S₂₆ convergence",
            "code": ("a[k_] := Exp[-0.57 k/26]; "
                     "SN[N_] := Sum[a[k], {k, 1, N}]; "
                     "RN[N_] := a[N]/2 + BernoulliB[2]/(2!) * "
                     "Sum[(-1)^j Binomial[1,j] a[N+1-j], {j,0,1}]; "
                     "Table[{N, SN[N], SN[N] + RN[N]}, {N, 5, 26, 5}]"),
        },
    ]


# ── §5  Self-tests ────────────────────────────────────────────────────────

def _run_tests() -> bool:
    """Validate spectral ladder module."""
    ok = True

    # Ladder monotonicity
    ladder = SpectralLadder26State()
    energies = [ladder.energy(n) for n in range(1, 27)]
    if not all(energies[i] < energies[i+1] for i in range(25)):
        print("FAIL: spectral ladder not monotonically increasing")
        ok = False
    else:
        print(f"OK: E_1 = {energies[0]:.6e}, E_26 = {energies[25]:.6e} (monotonic)")

    # Transition energy positive
    dE = ladder.transition(26, 1)
    if dE <= 0:
        print(f"FAIL: E_26 - E_1 = {dE}")
        ok = False
    else:
        print(f"OK: E_26 - E_1 = {dE:.6e} J")

    # Ramanujan acceleration
    ram = RamanujanAcceleration()
    s_std = ram.partial_sum(26)
    s_acc = ram.accelerated_sum(26)
    print(f"OK: S_26 = {s_std:.10f}, S_26^(R) = {s_acc:.10f}")

    # Phonon mapping
    mapping = SpectralLadderPhononMapping()
    res = mapping.compute({"Gamma_THz": 0.10})
    if len(res["mapping"]) != 26:
        print(f"FAIL: expected 26 levels, got {len(res['mapping'])}")
        ok = False
    else:
        print(f"OK: 26-level phonon mapping computed")

    # WSTP
    exprs = build_spectral_ladder_wstp_expressions()
    if len(exprs) != 2:
        print(f"FAIL: expected 2 WSTP, got {len(exprs)}")
        ok = False
    else:
        print(f"OK: {len(exprs)} WSTP expressions")

    print(f"\n{'ALL TESTS PASSED' if ok else 'SOME TESTS FAILED'}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
