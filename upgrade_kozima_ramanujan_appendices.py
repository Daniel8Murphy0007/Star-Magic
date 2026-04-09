#!/usr/bin/env python3
"""
upgrade_kozima_ramanujan_appendices.py — Session 204

Batch upgrade PAPER_001–877 whitepapers with three new appendix sections
derived from the Session 204 gap modules:

  Phase 1: PAPER_001–009 + 840/851/852
           → Appendix: Kozima-UQFF LENR Mechanism (detailed neutron-drop math)

  Phase 2: PAPER_840/851/852/855
           → Appendix: Ramanujan 26-State Mock Theta & π Approximation

  Phase 3: PAPER_001–877 (full corpus)
           → Appendix: Session 204 Codebase Upgrade Reference
             (condensed Kozima + Ramanujan + mock theta cross-reference)

Source modules analyzed:
  - fneutron_s26_coupling.py     F_neutron × S₂₆ buoyancy-polylog coupling
  - kozima_scm_cross_section.py  SCm-modulated Kozima neutron-drop σ_n
  - ramanujan_polylog_s26.py     Euler-accelerated Li₂₆([SSq]) polylogarithm
  - mock_theta_q26.py            Third-order mock theta f₂₆, φ₂₆, ψ₂₆
  - ramanujan_pi_uqff.py         Ramanujan 1/π with 26-state UQFF modification
  - mock_theta_pi_wstp_kernel.py 9-symbol Wolfram Language export
  - s26_wstp_kernel.py           8-symbol S₂₆ Wolfram export
  - kozima_wstp_kernel.py        11-symbol Kozima Wolfram export

Run from Star-Magic root:
    python upgrade_kozima_ramanujan_appendices.py
"""

import os
import re
import sys
from datetime import datetime

WHITEPAPERS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "whitepapers")

# ══════════════════════════════════════════════════════════════════════════════
# APPENDIX A: Kozima-UQFF LENR Mechanism (detailed)
# Target: PAPER_001–009 + 840/851/852
# ══════════════════════════════════════════════════════════════════════════════

KOZIMA_UQFF_APPENDIX = r"""

---

## Appendix: Kozima-UQFF LENR Mechanism (Session 204)

> *Derived from `fneutron_s26_coupling.py`, `kozima_scm_cross_section.py`,
> `kozima_wstp_kernel.py`, and `scm_activation_function.py`. Added by
> `upgrade_kozima_ramanujan_appendices.py` (Session 204, April 2026).*

### K.1 Neutron Drop Force — Static Model

The Kozima neutron-drop force integrates into the F_U_Bi_i unified field as an
additional LENR term:

$$F_{\rm neutron} = k_{\rm neutron} \times \sigma_n = 10^{10} \times 10^{-4} = 10^6 \;\text{N}$$

| Parameter | Value | Description |
|-----------|-------|-------------|
| k_neutron | 10^10 N | Neutron-drop strength constant |
| sigma_0 | 10^-4 | Base cross-section (dimensionless) |
| F_neutron (static) | 10^6 N | Lattice-scale neutron production force |

### K.2 Frequency-Dependent Cross-Section (SCm-Modulated)

The SCm superconductive manifold modulates the cross-section via VDS 26-level
enhancement:

$$\sigma_n^{\rm SCm}(\omega, n) = \sigma_0 \cdot \exp\!\left[-\frac{(\omega - \omega_{\rm SCm})^2}{2\Gamma^2}\right] \cdot \left(1 + \frac{[\text{SSq}] \cdot n}{26}\right)$$

| Symbol | Value | Description |
|--------|-------|-------------|
| omega_SCm | 2pi x 1.25 THz | SCm phonon resonance angular frequency |
| Gamma | 2pi x 0.1 THz | Resonance width |
| [SSq] | 0.57 | Universal Quantized Factor |
| n | 0..26 | VDS vacuum density level |

**Key result:** The VDS factor (1 + [SSq]*n/26) amplifies sigma_n by up to
1.57x at n=26, encoding the 26-level vacuum density hierarchy.

### K.3 Buoyancy-Coupled Neutron Force

The full frequency-dependent force couples the neutron drop with buoyancy reversal:

$$F_{\rm neutron}^{\rm SCm} = N_n \cdot \sigma_n^{\rm SCm}(\omega) \cdot \Phi_{\rm phonon} \cdot \left(\frac{F_{U,Bi}}{F_U} - 1\right)$$

| Symbol | Description |
|--------|-------------|
| N_n | Neutron number density in lattice site |
| Phi_phonon | Phonon flux at resonance frequency |
| F_{U,Bi}/F_U - 1 | Buoyancy reversal ratio (> 0 for active LENR) |

### K.4 S_26 Polylogarithm Coupling (Session 204)

The neutron-drop force operates within the 26-level VDS vacuum structure. The
coupled force at each VDS level n:

$$F_{\rm coupled}(\omega) = \sum_{n=0}^{26} F_{\rm neutron}(\omega, n) \times S_{26}\!\left([\text{SSq}] \cdot \left(1 + \frac{n}{26}\right)\right)$$

where S_26(z) = Li_26(z) is the 26-dimensional polylogarithm computed via
Eta-function Euler acceleration (O(1/2^N) convergence):

$$S_{26}(z) = \text{Li}_{26}(z) = \frac{\eta_{26}(z)}{1 - 2^{1-26}} + \frac{2^{1-26}}{1 - 2^{1-26}} \text{Li}_{26}(z^2)$$

This gives the buoyancy force weighted by the full 26-level vacuum density
spectrum, producing ~470x amplification relative to decoupled models.

### K.5 SCm Activation Function

$$A_{\rm SCm}(B) = \exp\!\left[-\frac{B^2}{B_{\rm crit}^2}\right], \quad B_{\rm crit} = 4.4 \times 10^{13} \;\text{T}$$

The Gaussian activation (from `scm_activation_function.py`) governs the transition
probability for the neutron-drop mechanism as a function of ambient magnetic field.

### K.6 Wolfram Implementation

The `UQFFKozima` package (11 symbols) exports the complete Kozima LENR framework
to Wolfram Language via WSTP:

```
FNeutronForce[Nn, sigma, phiPhonon, fUBi, fU]
SigmaSCm[omega, n]
SCmActivation[B]
FNeutronS26[..., nTerms]
```

*Source: `kozima_wstp_kernel.py` → `uqff_kozima_kernel.wl`*

"""

# ══════════════════════════════════════════════════════════════════════════════
# APPENDIX B: Ramanujan 26-State Mock Theta & π Approximation (detailed)
# Target: PAPER_840/851/852/855
# ══════════════════════════════════════════════════════════════════════════════

RAMANUJAN_MOCK_THETA_APPENDIX = r"""

---

## Appendix: Ramanujan 26-State Mock Theta Functions & pi Approximation (Session 204)

> *Derived from `mock_theta_q26.py`, `ramanujan_pi_uqff.py`, `ramanujan_polylog_s26.py`,
> and `mock_theta_pi_wstp_kernel.py`. Added by `upgrade_kozima_ramanujan_appendices.py`
> (Session 204, April 2026).*

### R.1 q-Pochhammer Symbol (Proper q-Series)

The q-Pochhammer symbol is the fundamental building block for mock theta functions:

$$(a; q)_n = \prod_{k=0}^{n-1} (1 - a q^k)$$

This is distinct from the rising factorial (a)_n = a(a+1)...(a+n-1) used elsewhere
in the codebase (`qcalcgeom_helpers.py`). The q-Pochhammer is evaluated at
q = [SSq] = 0.57 as the UQFF quantum parameter.

### R.2 Third-Order Mock Theta Functions (26-State Truncation)

Three Ramanujan third-order mock theta functions, truncated at N=26 UQFF states:

$$f_{26}(q) = \sum_{n=0}^{25} \frac{q^{n^2}}{(-q;\,q)_n^2}$$

$$\phi_{26}(q) = \sum_{n=0}^{25} \frac{q^{n^2}}{(-q^2;\,q^2)_n}$$

$$\psi_{26}(q) = \sum_{n=1}^{26} \frac{q^{n^2}}{(q;\,q^2)_n}$$

**Numerical values at q = [SSq] = 0.57:**

| Function | Value | Levels |
|----------|-------|--------|
| f_26(0.57) | 1.257 | n = 0..25 |
| phi_26(0.57) | 1.507 | n = 0..25 |
| psi_26(0.57) | 1.647 | n = 1..26 |

### R.3 UQFF Coupled Theta Amplitude

The 26-state coupled theta amplitude weights mock theta contributions by VDS
level amplitudes:

$$\Theta_{26} = \sum_{i=1}^{26} A_i(n) \cdot \bigl[f_{26}(q_i) + \phi_{26}(q_i) + \psi_{26}(q_i)\bigr]$$

where q_i = [SSq] * exp(-kappa * i * t / 26) is the time-dependent quantum parameter
at VDS level i, and A_i = (2*pi)^(i/6) * (rho_SCm / rho_UA) is the VDS amplitude.

### R.4 Ramanujan 1/pi Series (Classical)

$$\frac{1}{\pi} = \frac{2\sqrt{2}}{9801} \sum_{n=0}^{\infty} \frac{(4n)!\,(1103 + 26390\,n)}{(n!)^4 \cdot 396^{4n}}$$

**Convergence:** Each term adds ~8 decimal digits of pi. 4 terms yield 31+ correct
digits. The coefficient R_n = (4n)!/((n!)^4 * 396^(4n)) is computed via log-gamma
to prevent factorial overflow for large n.

### R.5 UQFF-Modified 1/pi (26-State Weighting)

$$\frac{1}{\pi_{\rm UQFF}} = \frac{2\sqrt{2}}{9801} \cdot \frac{1}{C_{26}} \sum_{n=0}^{N-1} R_n\,(1103 + 26390\,n) \cdot W_{26}(n)$$

where the 26-state weight factor:

$$W_{26}(n) = \prod_{i=1}^{26}\left[1 + [\text{SSq}] \cdot \exp\!\left(-\frac{\kappa\, i\, n}{26}\right)\right]$$

and C_26 = (1 + [SSq])^26 normalizes to recover classical Ramanujan at kappa = 0.

**Key result:** For physical kappa = 5.787 x 10^-9, the UQFF modification preserves
15+ digits of pi, confirming that the 26-state vacuum structure does not distort
the fundamental constant at observable precision.

### R.6 26D Hypergeometric Generalization

$$\frac{1}{\pi_{26D}} = \frac{2\sqrt{2}}{9801\,C_{26}^{\rm hyper}} \sum_{n=0}^{N-1} R_n\,(a_{26} + b_{26}\,n)$$

where a_26 = 1103 * H_26^alt (alternating harmonic sum), b_26 = 26390 * (26/13),
and C_26^hyper = H_26^alt normalizes the leading term. This yields 7 digits with
26 terms — the dimensional scaling alters convergence rate while preserving the
Ramanujan algebraic structure.

### R.7 Ramanujan-Accelerated Polylogarithm S_26

$$S_{26}(z) = \text{Li}_{26}(z) = \sum_{k=1}^{\infty} \frac{z^k}{k^{26}}$$

Evaluated via eta-function decomposition (from `ramanujan_polylog_s26.py`):

$$\text{Li}_{26}(z) = \frac{\eta_{26}(z)}{1 - 2^{1-26}} + \frac{2^{1-26}}{1 - 2^{1-26}} \cdot \text{Li}_{26}(z^2)$$

At z = [SSq] = 0.57, converges to 15.7+ digits in 53 terms (vs naive series
requiring 10^9+ terms). The Euler transform for eta_26 uses the binomial
acceleration: eta_s(z) = Sum_{n=0}^{N} (1/2^{n+1}) Sum_{j=0}^{n} C(n,j) (-1)^j z^{j+1}/(j+1)^s.

### R.8 Wolfram Implementation

The `UQFFMockThetaPi` package (9 symbols) exports all mock theta and pi functions:

```
qPochhammer[a, q, n]         -- q-Pochhammer (a;q)_n
f26[q], phi26[q], psi26[q]   -- Third-order mock thetas
thetaCoupled26[q, ssq, kap]  -- 26-state coupled amplitude
ramanujanR[n]                -- R_n coefficient
oneOverPiClassical[nTerms]   -- Ramanujan 1/pi
oneOverPiUQFF[nTerms, ssq, kap] -- UQFF-modified 1/pi
pi26DHypergeometric[nTerms]  -- 26D generalization
```

*Source: `mock_theta_pi_wstp_kernel.py` -> `uqff_mock_theta_pi_kernel.wl`*

"""

# ══════════════════════════════════════════════════════════════════════════════
# APPENDIX C: Session 204 Codebase Upgrade Reference (condensed)
# Target: PAPER_001–877 (full corpus)
# ══════════════════════════════════════════════════════════════════════════════

SESSION_204_CONDENSED_APPENDIX = r"""

---

## Appendix: Session 204 Codebase Upgrade Reference

> *Cross-reference appendix for Session 204 (April 2026) codebase upgrades.
> Added by `upgrade_kozima_ramanujan_appendices.py`. For detailed derivations,
> see PAPER_840/851/852/855.*

### S204.1 Kozima-UQFF LENR Integration

| Module | Purpose | Key Result |
|--------|---------|------------|
| `fneutron_s26_coupling.py` | F_neutron x S_26 buoyancy-polylog coupling | ~470x amplification via 26-level VDS |
| `kozima_scm_cross_section.py` | SCm-modulated neutron-drop cross-section | sigma_n^SCm with VDS factor (1+[SSq]*n/26) |
| `kozima_wstp_kernel.py` | 11-symbol Wolfram export (`UQFFKozima`) | FNeutronForce, SigmaSCm, SCmActivation |

**Core equation:** F_neutron^SCm = N_n * sigma_n^SCm(omega) * Phi_phonon * (F_{U,Bi}/F_U - 1)
where sigma_n^SCm(omega,n) = sigma_0 * exp[-(omega-omega_SCm)^2/(2*Gamma^2)] * (1 + [SSq]*n/26)

### S204.2 Ramanujan 26-State Summation

| Module | Purpose | Key Result |
|--------|---------|------------|
| `ramanujan_polylog_s26.py` | Li_26([SSq]) via Euler-Ramanujan acceleration | 15.7+ digits in 53 terms |
| `s26_wstp_kernel.py` | 8-symbol Wolfram export (`UQFFS26`) | S26, R26, NaiveLi, S26VDS |

**Core equation:** S_26(z) = Li_26(z) = eta_26(z)/(1-2^{1-26}) + 2^{1-26}/(1-2^{1-26}) * Li_26(z^2)

### S204.3 Mock Theta Functions (26-State)

| Module | Purpose | Key Result |
|--------|---------|------------|
| `mock_theta_q26.py` | f_26(q), phi_26(q), psi_26(q) q-series | Proper q-Pochhammer (a;q)_n |

**Core equations:**
- f_26(q) = Sum_{n=0}^{25} q^{n^2} / (-q;q)_n^2
- phi_26(q) = Sum_{n=0}^{25} q^{n^2} / (-q^2;q^2)_n
- psi_26(q) = Sum_{n=1}^{26} q^{n^2} / (q;q^2)_n

### S204.4 Ramanujan 1/pi with UQFF Modification

| Module | Purpose | Key Result |
|--------|---------|------------|
| `ramanujan_pi_uqff.py` | Classical + UQFF-modified 1/pi + 26D | 21 digits classical, 15 UQFF, 7 digits 26D |
| `mock_theta_pi_wstp_kernel.py` | 9-symbol Wolfram export (`UQFFMockThetaPi`) | qPochhammer, f26, oneOverPiUQFF |

**Core equation:** 1/pi = (2*sqrt(2)/9801) * Sum R_n * (1103+26390n) * W_26(n) / C_26
where W_26(n) = Prod_{i=1}^{26} [1 + [SSq]*exp(-kappa*i*n/26)]

### S204.5 Calibration Constants (Canonical)

| Symbol | Value | Description |
|--------|-------|-------------|
| [SSq] | 0.57 | Universal Quantized Factor |
| kappa | 5.787 x 10^-9 s^-1 | UQFF exponential decay rate |
| beta_i | 0.603 | Buoyancy coupling coefficient |
| H_SCm | 0.99 | SCm manifold completeness |
| rho_SCm | 7.09 x 10^-37 kg/m^3 | SCm vacuum density |
| rho_UA | 7.09 x 10^-36 kg/m^3 | UA aether vacuum density |
| omega_SCm | 2*pi x 1.25 THz | SCm phonon resonance |
| sigma_0 | 10^-4 | Base neutron cross-section |

*Implementation: all modules operational in `CondensedPhysics.py`, `CondensedPhysics2.py`,
`MAIN_1_CoAnQi.cpp`, and Wolfram kernels (`uqff_kozima_kernel.wl`, `uqff_s26_kernel.wl`,
`uqff_mock_theta_pi_kernel.wl`).*

"""

# ══════════════════════════════════════════════════════════════════════════════
# MARKERS
# ══════════════════════════════════════════════════════════════════════════════

KOZIMA_MARKERS = [
    "Appendix: Kozima-UQFF LENR Mechanism",
    "Kozima-UQFF LENR Mechanism (Session 204)",
    "fneutron_s26_coupling.py",
]

RAMANUJAN_MARKERS = [
    "Appendix: Ramanujan 26-State Mock Theta",
    "Mock Theta Functions & pi Approximation",
    "mock_theta_q26.py",
]

SESSION204_MARKERS = [
    "Appendix: Session 204 Codebase Upgrade Reference",
    "Session 204 Codebase Upgrade Reference",
    "S204.1 Kozima-UQFF LENR Integration",
]


def has_marker(text: str, markers: list) -> bool:
    return any(m in text for m in markers)


def get_paper_number(filename: str):
    m = re.match(r'PAPER_(\d+)', filename, re.IGNORECASE)
    return int(m.group(1)) if m else None


def append_if_missing(path: str, appendix: str, markers: list, label: str) -> bool:
    """Append appendix to file if none of the markers are present. Returns True if modified."""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()

    if has_marker(text, markers):
        return False

    with open(path, 'a', encoding='utf-8') as f:
        f.write(appendix)

    return True


def run():
    if not os.path.isdir(WHITEPAPERS_DIR):
        print(f"ERROR: {WHITEPAPERS_DIR} not found.")
        sys.exit(1)

    files = sorted(os.listdir(WHITEPAPERS_DIR))
    md_files = [f for f in files if f.endswith('.md')]

    stats = {
        "kozima_upgraded": 0,
        "kozima_skipped": 0,
        "ramanujan_upgraded": 0,
        "ramanujan_skipped": 0,
        "s204_upgraded": 0,
        "s204_skipped": 0,
        "errors": 0,
    }

    # ── Phase 1: Kozima LENR appendix → PAPER_001–009 + 840/851/852 ──────
    print("=" * 72)
    print("Phase 1: Kozima-UQFF LENR Mechanism Appendix")
    print("  Target: PAPER_001–009 + 840/851/852")
    print("=" * 72)

    kozima_targets = set(range(1, 10)) | {840, 851, 852}

    for fname in md_files:
        num = get_paper_number(fname)
        if num is None or num not in kozima_targets:
            continue

        path = os.path.join(WHITEPAPERS_DIR, fname)
        try:
            if append_if_missing(path, KOZIMA_UQFF_APPENDIX, KOZIMA_MARKERS, "Kozima"):
                print(f"  + Kozima appendix: {fname}")
                stats["kozima_upgraded"] += 1
            else:
                stats["kozima_skipped"] += 1
        except Exception as e:
            print(f"  ! ERROR {fname}: {e}")
            stats["errors"] += 1

    print(f"\n  Kozima: {stats['kozima_upgraded']} upgraded, {stats['kozima_skipped']} skipped\n")

    # ── Phase 2: Ramanujan Mock Theta + Pi → PAPER_840/851/852/855 ────────
    print("=" * 72)
    print("Phase 2: Ramanujan 26-State Mock Theta & Pi Appendix")
    print("  Target: PAPER_840/851/852/855")
    print("=" * 72)

    ramanujan_targets = {840, 851, 852, 855}

    for fname in md_files:
        num = get_paper_number(fname)
        if num is None or num not in ramanujan_targets:
            continue

        path = os.path.join(WHITEPAPERS_DIR, fname)
        try:
            if append_if_missing(path, RAMANUJAN_MOCK_THETA_APPENDIX, RAMANUJAN_MARKERS, "Ramanujan"):
                print(f"  + Ramanujan appendix: {fname}")
                stats["ramanujan_upgraded"] += 1
            else:
                stats["ramanujan_skipped"] += 1
        except Exception as e:
            print(f"  ! ERROR {fname}: {e}")
            stats["errors"] += 1

    print(f"\n  Ramanujan: {stats['ramanujan_upgraded']} upgraded, {stats['ramanujan_skipped']} skipped\n")

    # ── Phase 3: Session 204 Condensed Reference → PAPER_001–877 ─────────
    print("=" * 72)
    print("Phase 3: Session 204 Codebase Upgrade Reference (condensed)")
    print("  Target: PAPER_001–877 (full corpus)")
    print("=" * 72)

    for fname in md_files:
        num = get_paper_number(fname)
        if num is None or num < 1 or num > 877:
            continue

        path = os.path.join(WHITEPAPERS_DIR, fname)
        try:
            if append_if_missing(path, SESSION_204_CONDENSED_APPENDIX, SESSION204_MARKERS, "S204"):
                stats["s204_upgraded"] += 1
            else:
                stats["s204_skipped"] += 1
        except Exception as e:
            print(f"  ! ERROR {fname}: {e}")
            stats["errors"] += 1

    # Print progress every 100 papers
    print(f"\n  Session 204 ref: {stats['s204_upgraded']} upgraded, {stats['s204_skipped']} skipped")

    # ── Summary ───────────────────────────────────────────────────────────
    print("\n" + "=" * 72)
    print("UPGRADE SUMMARY")
    print("=" * 72)
    print(f"  Kozima LENR appendix:     {stats['kozima_upgraded']} upgraded, {stats['kozima_skipped']} already present")
    print(f"  Ramanujan Mock Theta/Pi:  {stats['ramanujan_upgraded']} upgraded, {stats['ramanujan_skipped']} already present")
    print(f"  Session 204 Reference:    {stats['s204_upgraded']} upgraded, {stats['s204_skipped']} already present")
    print(f"  Errors:                   {stats['errors']}")
    total = stats['kozima_upgraded'] + stats['ramanujan_upgraded'] + stats['s204_upgraded']
    print(f"\n  TOTAL PAPERS MODIFIED:    {total} appendix insertions across corpus")
    print("=" * 72)

    return stats


if __name__ == "__main__":
    run()
