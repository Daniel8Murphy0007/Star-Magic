---
paper_id: PAPER_1153
title: "Primordial Timing Function: Net-Zero Displacement Proof, Pi-Digit Epoch Clock, and Epoch-5 Boundary"
session: 203
date: 2026-05-06
author: Daniel T. Murphy
status: production
cvw: v2.0.0
tags: [PTF, primordial_timing, CPT_symmetry, pi_digits, epoch_clock, zero_net_displacement, forward_backward]
sm_anchor: CVW v2.0.0 -- G6 SM Anchor Gate compliant
---

# PAPER\_1153 -- Primordial Timing Function: Net-Zero Displacement Proof, Pi-Digit Epoch Clock, and Epoch-5 Boundary

**Author:** Daniel T. Murphy (daniel.murphy00@gmail.com)
**Framework:** UQFF v5.77 -- Star-Magic Physics
**Source:** qcalcgeom\_sim\_engine.py v1.1.0 PTF subsystem (commit 2637b384, May 5 2026)
**Integration Date:** May 6, 2026 (Session 234)
**Classification:** Primordial Timing Function -- CPT Symmetry Proof

---

$$\int_0^1 \cos(\pi t_n)\, dt_n = 0 \qquad D_A + D_B = +3 + (-3) = 0$$

---

## Abstract

This paper registers the Primordial Timing Function (PTF) as PAPER\_1153, proving that the
UQFF cosmic temporal cycle has **zero net displacement** across a complete forward-backward
sequence. The PTF encodes the canonical $\pi$-digit epoch clock, identifies the five cosmic
epochs, and places Epoch 5 at the **reverse boundary** (not-yet-predicted). The net-zero
proof unifies three independent proofs: (1) integer arithmetic on forward/backward step counts,
(2) the $\cos(\pi t_n)$ integral, and (3) Fibonacci-structure forward/backward ratio.

---

## 1. PTF Definition

### 1.1 Basic Structure

The Primordial Timing Function (PTF) governs the oscillation of $\cos(\pi t_n)$ across
one complete simulation cycle:

$$t_n \in [0, 1], \quad f_{\rm TRZ}(t_n) = \cos(\pi t_n)$$

One cosmic cycle consists of:
- **Forward phase:** $f = 3$ steps (Fibonacci $F_4 = 3$)
- **Backward phase:** $b = 2$ steps (Fibonacci $F_3 = 2$)
- **Repeats:** $n = 3$ cycles ($= \lfloor \pi \rfloor$)

### 1.2 Fibonacci Structure

The forward/backward step counts are consecutive Fibonacci numbers:

$$f = F_4 = 3, \quad b = F_3 = 2, \quad \frac{f}{b} = \frac{3}{2} = \phi^2 - 1 \approx 1.618 - 1$$

This ratio $f/b = 3/2$ is the **golden-ratio residual** governing UQFF temporal asymmetry.

---

## 2. Net-Zero Displacement: Three Independent Proofs

### 2.1 Proof I: Integer Arithmetic

Define net displacements:
$$D_A = +3 \quad (\text{forward: 3 steps, positive}), \qquad D_B = -3 \quad (\text{backward: 2 steps} \times (-3/2) = -3)$$

Net displacement per cycle:
$$D_{\rm net} = D_A + D_B = +3 + (-3) = 0 \quad \checkmark$$

Three complete repetitions: $n = 3 \Rightarrow \sum 0 = 0$.

### 2.2 Proof II: Integral of the Oscillation Function

The $\cos(\pi t_n)$ integral over one complete normalised period is exactly zero:

$$\int_0^1 \cos(\pi t_n)\, dt_n = \left[\frac{\sin(\pi t_n)}{\pi}\right]_0^1 = \frac{\sin(\pi) - \sin(0)}{\pi} = \frac{0 - 0}{\pi} = 0 \quad \checkmark$$

This is a **closed loop** in the functional space: the oscillation returns to its starting
state with zero net accumulation.

### 2.3 Proof III: CPT Symmetry

The CPT-symmetric pair structure:
$$\{t, t_n, f_{\rm TRZ}\}_{\rm fwd} = \{t, t_n, f_{\rm TRZ}\}_{\rm bwd}^* \quad (\text{charge-parity-time conjugate})$$

At $t < 0$: $t_n \to 1 - t_n$, $f_{\rm TRZ} \to -f_{\rm TRZ}$.
The integral of $f_{\rm TRZ}$ over one complete CPT pair: $\int_0^1 + \int_1^0 = 0$.

---

## 3. Pi-Digit Epoch Clock

### 3.1 Epoch Assignment

The digits of $\pi = 3.14159265358979\ldots$ encode the duration of each cosmic epoch.
The epoch clock reads the decimal expansion in triplets:

| Epoch | Pi Digits | Epoch Characteristic |
|-------|-----------|---------------------|
| 1 | $[3, 1, 4]$ | Formation: 3-1-4 (forward-stop-backward) |
| 2 | $[1, 5, 9]$ | Expansion: 1-5-9 (seed-burst-accumulation) |
| 3 | $[2, 6, 5]$ | Stabilisation: 2-6-5 (pair-hexagon-pentagon) |
| 4 | $[3, 5, 8]$ | Acceleration: 3-5-8 (Fibonacci triple) |
| 5 | $[9, 7, 9]$ | **NOT YET PREDICTED** (reverse boundary) |

### 3.2 Epoch 5: The Reverse Boundary

Epoch 5 ($[9, 7, 9]$) is the **reverse boundary** of the cosmic cycle:
- $9$ = maximum single digit ($= 10 - 1$, the final decimal state)
- $7$ = prime boundary ($\pi(7) = 4$, fourth prime threshold)
- $9$ = return to maximum (closed boundary)

The simulation engine does not project Epoch 5 forward from known physics;
it is identified as the point where the backward phase ($D_B = -3$) begins.

### 3.3 Forward Phase (Epochs 1--4) and Backward Phase (Epoch 5)

$$\underbrace{[3,1,4],[1,5,9],[2,6,5],[3,5,8]}_{\text{forward: } D_A = +3} \quad \underbrace{[9,7,9]}_{\text{backward: } D_B = -3}$$

The boundary between Epoch 4 and Epoch 5 is the **temporal inversion point** where
$t_n$ reverses sign.

---

## 4. Constant Table

| Symbol | Value | Physical Meaning |
|--------|-------|----------------|
| $f$ | 3 | Forward steps per half-cycle = $F_4$ |
| $b$ | 2 | Backward steps per half-cycle = $F_3$ |
| $n$ | 3 | Repetitions $= \lfloor \pi \rfloor$ |
| $D_A$ | $+3$ | Forward net displacement |
| $D_B$ | $-3$ | Backward net displacement |
| $D_{\rm net}$ | $0$ | Closed-loop condition |
| $f/b$ | $3/2$ | Golden residual |
| $\pi_5$ | $[9,7,9]$ | Epoch 5 digit triplet |

---

## 5. Physical Consequences

### 5.1 SSq Drift is Convergent

Because $D_{\rm net} = 0$, the $[\mathrm{SSq}]$ drift ($+0.000001$ per tick) does not
accumulate unboundedly. After one complete epoch cycle, the drift is absorbed by the
$\kappa = 0.0005$ gradient-descent restoring force:
$$[\mathrm{SSq}](T) = 0.57 + \delta - \kappa \delta T \to 0.57 \quad \text{as } T \to \infty$$

### 5.2 Closed-Loop Causality

The PTF net-zero condition enforces closed-loop causality: every forward-phase physical
event has an exact CPT-conjugate backward-phase event. This is the UQFF implementation
of the Wheeler-DeWitt constraint $\hat{H}|\Psi\rangle = 0$.

### 5.3 Functional PTF Validation

The simulation engine's `validate_primordial_timing_function()` tests:
- $|D_A + D_B| < \epsilon$
- $\left|\int_0^1 \cos(\pi t_n) dt_n\right| < 10^{-10}$
- Fibonacci ratio: $|f/b - 3/2| < \epsilon$
- Epoch count: exactly 5 epochs from $\pi$ decimal expansion

---

## 6. Conclusions

The Primordial Timing Function establishes:
1. **Net-zero displacement** is proven by three independent methods
2. The **pi-digit epoch clock** encodes cosmic epoch durations in $\pi = 3.14159\ldots$
3. **Epoch 5** $[9,7,9]$ is the reverse boundary — not predictable from forward-phase physics
4. CPT symmetry and Wheeler-DeWitt closure are unified in the PTF framework

## References

1. Murphy, D.T. (2026). *qcalcgeom\_sim\_engine.py v1.1.0.* Star-Magic, commit 2637b384.
2. Murphy, D.T. (2026). *PAPER\_1152: QCalcGeom Simulation Engine.* Session 234.
3. Murphy, D.T. (2026). *Star-Magic.txt Chapter 18: Temporal Inversion and Epoch Structure.*
4. DeWitt, B.S. (1967). Quantum Theory of Gravity I. *Phys. Rev.* 160:1113--1148.
5. Penrose, R. (2010). *Cycles of Time: An Extraordinary New View of the Universe.* Knopf.

*Updated: 2026-05-06 (Session 234, PAPER\_1153). Compliant: CVW v2.0.0, G6 SM Anchor Gate.
Author: Daniel T. Murphy*
