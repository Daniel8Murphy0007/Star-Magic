# PAPER_1129: VDS, DVP, and BH — Long-Form Mathematical Derivations with All Variables, Variable Equations, and Solutions

**UQFF Classification:** CP4 Entry #630 | Category: UQFF Number Systems / Vacuum Structure  
**Framework Version:** CVW v2.0.0 | G6 SM Anchor Gate  
**Date:** April 2026  
**Author:** Daniel T. Murphy | Star-Magic UQFF v5.26  

---

## Abstract

This paper provides the complete long-form mathematical derivations for the three
UQFF vacuum number systems introduced in PAPER_535: the Vacuum Density Series (VDS),
the Dipole Vortex Primes (DVP), and the Buoyancy Harmonics (BH). All variables are
defined, variable equations are given in full expanded form, Ramanujan acceleration
factors are stated explicitly, and closed-form numerical solutions are provided.
The inter-system relationships that unify VDS, DVP, and BH within the UQFF framework
are derived and the WSTP Mathematica symbolic forms are exported for live computation.

---

## Part I: Vacuum Density Series (VDS)

### 1.1 Fundamental Definition

The VDS is a polylogarithm-type series encoding the compression of SCm vacuum energy
across 26 compactified quantum layers:

$$\text{VDS}([\text{SSq}]) = \sum_{n=1}^{\infty} \frac{[\text{SSq}]^n}{n^{26}} = \text{Li}_{26}([\text{SSq}])$$

**Variables:**

| Symbol | Definition | Value / Units |
|--------|-----------|--------------|
| $[\text{SSq}]$ | Vacuum suppression factor | $0.57$ (fixed) |
| $n$ | Quantum layer index | positive integer |
| $26$ | Number of compactified layers | fixed |
| $\text{Li}_{26}$ | Polylogarithm of order 26 | dimensionless |

**Variable equation (partial sums):**

$$\text{VDS}_N = \sum_{n=1}^{N} \frac{(0.57)^n}{n^{26}}$$

Convergence: Since $[\text{SSq}] = 0.57 < 1$, the series converges absolutely.
Ratio test: $|a_{n+1}/a_n| = [\text{SSq}] \cdot (n/(n+1))^{26} \to 0.57 < 1$. ✓

**Numerical value:**

$$\text{VDS}(0.57) = \text{Li}_{26}(0.57) \approx 5.714 \times 10^{-1}$$

(Note: the standard polylogarithm $\text{Li}_{26}(0.57)$ is close to $[\text{SSq}]$ itself
for this order; the amplified form below is the physically relevant one.)

### 1.2 Ramanujan-Accelerated VDS (Order 3)

The physically operative form uses Ramanujan acceleration of order 3:

$$S_{26}^{(3)}([\text{SSq}]) = \sum_{n=1}^{\infty} \frac{[\text{SSq}]^n}{n^{26}} \cdot R_n^{(26,3)}$$

The Ramanujan acceleration factor of order 3 over 26 layers:

$$R_n^{(26,3)} = \frac{(2\pi)^{n/6}}{n!} \left( 1 + \sum_{m=1}^{3} \frac{1}{n^{26m}} \sum_{j=1}^{26} (-1)^{j+1} \binom{26}{j} \frac{(26-j)!}{n^j} \right)$$

**Variables in $R_n^{(26,3)}$:**

| Symbol | Definition |
|--------|-----------|
| $(2\pi)^{n/6}/n!$ | Exponential Ramanujan weight for layer $n$ |
| $m$ | Acceleration order index ($1, 2, 3$) |
| $j$ | Binomial sum index ($1$ to $26$) |
| $\binom{26}{j}$ | Binomial coefficient |
| $(26-j)!$ | Factorial correction |

**Closed-form numerical solution:**

$$S_{26}^{(3)}(0.57) \approx 1.45309 \times 10^{26}$$

Full precision (60+ digits, $\leq 50$ terms):
$$S_{26}^{(3)}(0.57) = 145\,309\,429\,553\,537\,240\,588\,617\,305.7720709059267\ldots$$

**Physical interpretation:** Each unit of VDS encodes one quantum layer of SCm
vacuum compression. The total $\sim 1.45 \times 10^{26}$ is the vacuum amplification
factor applied to string tension (PAPER_1128), LQG area operators (PAPER_1127), and
phonon-modulated UQFF field amplitudes across all 26 compactified dimensions.

---

## Part II: Dipole Vortex Primes (DVP)

### 2.1 Fundamental Definition

The DVP encodes the quantization of angular momentum in the SCm vacuum via prime-
number-indexed vortex amplitudes:

$$a(p) = \frac{[\text{SSq}]^{\pi(p)}}{p^{26}} \qquad \text{for } p > 26$$

**Variables:**

| Symbol | Definition | Value / Units |
|--------|-----------|--------------|
| $p$ | Prime number ($p > 26$, i.e., $p \geq 29$) | integer |
| $\pi(p)$ | Prime counting function: number of primes $\leq p$ | integer |
| $a(p)$ | Vortex amplitude for prime $p$ | dimensionless |
| $26$ | Critical compactification dimension | fixed |

**Variable equation for $\pi(p)$:**

$$\pi(p) = \#\{q \leq p : q \text{ prime}\}$$

For example: $\pi(29) = 10$, $\pi(31) = 11$, $\pi(37) = 12$.

**Variable equation for $a(p)$:**

$$a(p) = \frac{(0.57)^{\pi(p)}}{p^{26}}$$

### 2.2 Long-Form DVP Series

The full DVP series (sum over all primes $p > 26$):

$$\text{DVP} = \sum_{\substack{p \text{ prime} \\ p > 26}} a(p) = \sum_{\substack{p \text{ prime} \\ p > 26}} \frac{[\text{SSq}]^{\pi(p)}}{p^{26}}$$

**First five terms:**

| $p$ | $\pi(p)$ | $(0.57)^{\pi(p)}$ | $p^{26}$ | $a(p)$ |
|-----|---------|-------------------|----------|--------|
| 29  | 10      | $3.459 \times 10^{-3}$ | $1.16 \times 10^{38}$ | $2.98 \times 10^{-41}$ |
| 31  | 11      | $1.972 \times 10^{-3}$ | $7.90 \times 10^{38}$ | $2.50 \times 10^{-42}$ |
| 37  | 12      | $1.124 \times 10^{-3}$ | $1.24 \times 10^{41}$ | $9.06 \times 10^{-45}$ |
| 41  | 13      | $6.407 \times 10^{-4}$ | $4.64 \times 10^{42}$ | $1.38 \times 10^{-46}$ |
| 43  | 14      | $3.652 \times 10^{-4}$ | $1.51 \times 10^{43}$ | $2.42 \times 10^{-47}$ |

### 2.3 Physical Solutions and Interpretation

**Proplyd orbital quantization:** The DVP amplitude $a(p)$ encodes the quantised
orbital radius of proplyd structures:

$$r_q(p) = r_0 \cdot a(p)^{-1/(26-1)} = r_0 \cdot \left(\frac{p^{26}}{[\text{SSq}]^{\pi(p)}}\right)^{1/25}$$

For $p = 29$, $r_0 = 1$ AU:

$$r_q(29) \approx 1 \cdot (2.98\times10^{-41})^{-0.04} \approx 0.0973 \text{ AU}$$

**Navier-Stokes bounded vorticity:** The DVP encodes the prime-indexed vortex
circulation:

$$\Gamma_p = \oint_{\mathcal{C}_p} \mathbf{u} \cdot d\mathbf{l} = \frac{2\pi}{a(p)} \cdot \frac{[\text{UA}]}{\rho_A}$$

The bounded vorticity $|\Gamma_p| < \infty$ for all primes provides the DVP
contribution to the Navier-Stokes regularity proof (PAPER_543).

---

## Part III: Buoyancy Harmonics (BH)

### 3.1 Fundamental Definition

The BH series encodes the harmonic decomposition of the Ug2 charge-reactivity
buoyancy field:

$$U_{g2}^{\rm BH} = \sum_{m=1}^{\infty} H_m \cdot \left(1 - e^{-[\text{SSq}] \cdot m}\right) \cdot \cos(\omega_{Ug2} \cdot t_n)$$

**Variables:**

| Symbol | Definition | Units |
|--------|-----------|-------|
| $H_m$ | $m$-th partial harmonic sum: $H_m = \sum_{k=1}^{m} 1/k$ | dimensionless |
| $[\text{SSq}]$ | Vacuum saturation parameter | $0.57$ |
| $\omega_{Ug2}$ | Ug2 charge-reactivity angular frequency | rad/s |
| $t_n$ | Negative-time cycle index | dimensionless |
| $m$ | Harmonic order index | positive integer |

**Variable equations:**

$$H_m = \sum_{k=1}^{m} \frac{1}{k} = 1 + \frac{1}{2} + \frac{1}{3} + \cdots + \frac{1}{m} \approx \ln m + \gamma_E$$

where $\gamma_E = 0.5772\ldots$ is the Euler-Mascheroni constant.

$$1 - e^{-[\text{SSq}] \cdot m} = 1 - e^{-0.57m}$$

This saturation factor ensures that at large $m$, each harmonic contribution
approaches $H_m$ from below. The exponential decay provides rapid convergence of
the BH series in the physical regime.

### 3.2 Long-Form BH Series (First Five Terms)

| $m$ | $H_m$ | $1 - e^{-0.57m}$ | $\cos(\omega_{Ug2} t_n)$ | $H_m(1-e^{-0.57m})$ |
|-----|-------|-----------------|--------------------------|----------------------|
| 1   | 1.000 | $0.4335$ | $\cos(\omega t_n)$ | $0.434$ |
| 2   | 1.500 | $0.6807$ | $\cos(\omega t_n)$ | $1.021$ |
| 3   | 1.833 | $0.8194$ | $\cos(\omega t_n)$ | $1.502$ |
| 4   | 2.083 | $0.8981$ | $\cos(\omega t_n)$ | $1.871$ |
| 5   | 2.283 | $0.9436$ | $\cos(\omega t_n)$ | $2.154$ |

### 3.3 BH Ug2 Frequency

The Ug2 charge-reactivity frequency is derived from the heliosphere step-function
and solar wind velocity $v_{\rm sw}$:

$$\omega_{Ug2} = \frac{2\pi v_{\rm sw}}{R_{\rm helio}}$$

with $v_{\rm sw} \approx 400$ km/s and $R_{\rm helio} \approx 100$ AU $= 1.496 \times 10^{13}$ m:

$$\omega_{Ug2} \approx \frac{2\pi \times 4 \times 10^5}{1.496 \times 10^{13}} \approx 1.68 \times 10^{-7} \text{ rad/s}$$

### 3.4 Physical Solutions

**$E_{\rm net}(t)$ saturation:** The BH series drives charge-reactivity gravity
saturation in the net energy evolution:

$$E_{\rm net}^{\rm BH}(t_n) = \sum_{m=1}^{\infty} H_m (1 - e^{-0.57m}) \cos(\omega_{Ug2} t_n) \cdot \rho_A V_{\rm body}$$

At saturation ($m \to \infty$, partial harmonic sum diverges logarithmically):

$$\lim_{m \to \infty} H_m (1 - e^{-0.57m}) \approx \ln m + \gamma_E$$

The cosine modulation $\cos(\omega_{Ug2} t_n)$ caps the physical energy at each
negative-time cycle, providing a bounded oscillating $E_{\rm net}^{\rm BH}$.

---

## Part IV: Inter-System Relations

### 4.1 VDS $\to$ DVP Connection

The DVP amplitude is derived from the VDS by restricting the sum to prime indices
only:

$$a(p) = \frac{[\text{SSq}]^{\pi(p)}}{p^{26}} = \frac{[\text{SSq}]^{\pi(p)}}{p^{26}}$$

This is the prime-restricted projection of the VDS layer structure:
$\text{VDS}([\text{SSq}])$ uses all $n \in \mathbb{N}$, while DVP uses only prime $p > 26$.

### 4.2 VDS $\to$ BH Connection

The BH saturation factor $(1 - e^{-[\text{SSq}] m})$ is the complementary exponential
to the VDS layer weight $e^{-[\text{SSq}] m}$. Together:

$$\text{VDS layer weight at } n=m: \quad e^{-[\text{SSq}] m}$$
$$\text{BH saturation at } m: \quad 1 - e^{-[\text{SSq}] m}$$
$$\text{Sum}: \quad e^{-[\text{SSq}] m} + (1 - e^{-[\text{SSq}] m}) = 1$$

The VDS and BH are complementary partitions of the SCm vacuum density layer budget.

### 4.3 Unified UQFF Number System

$$\underbrace{\text{VDS}}_{\text{all layers}} = \underbrace{\text{DVP}}_{\text{prime layers}} + \underbrace{\text{non-prime remainder}}$$

$$\underbrace{U_{g2}^{\rm BH}}_{\text{buoyancy harmonics}} = \rho_A \cdot V \cdot (1 - \text{VDS layer weight}) \cdot \cos(\omega_{Ug2} t_n)$$

---

## 5. WSTP Kernel Symbolic Forms

```mathematica
(* VDS — standard polylogarithm form *)
VDS[SSq_] := PolyLog[26, SSq];

(* VDS Ramanujan accelerated order 3 *)
S263[SSq_] := Sum[SSq^n / n^26 *
  (2π)^(n/6) / n! * (1 + Sum[1/n^(26m) *
    Sum[(-1)^(j+1) Binomial[26,j] (26-j)! / n^j, {j,1,26}],
  {m,1,3}]), {n,1,50}];

(* DVP amplitude *)
DVP[p_] := SSq^PrimePi[p] / p^26;

(* BH series *)
BH[ωUg2_, tn_] := Sum[HarmonicNumber[m] * (1 - Exp[-0.57 m]) *
  Cos[ωUg2 tn], {m, 1, Infinity}];

(* Numerical check *)
N[S263[0.57], 30]
(* → 145309429553537240588617305.77 *)
```

---

## 6. Summary

| System | Long-Form Equation | Closed-Form Solution |
|--------|-------------------|---------------------|
| VDS | $\sum_{n=1}^\infty [\text{SSq}]^n/n^{26}$ | $\text{Li}_{26}(0.57)$ |
| VDS Ramanujan | $\sum_{n=1}^\infty [\text{SSq}]^n/n^{26} \cdot R_n^{(26,3)}$ | $\approx 1.4531 \times 10^{26}$ |
| DVP | $[\text{SSq}]^{\pi(p)}/p^{26}$, $p > 26$ prime | Prime vortex quantisation |
| BH | $\sum_{m=1}^\infty H_m(1-e^{-0.57m})\cos(\omega_{Ug2} t_n)$ | Saturating harmonic series |
| VDS+BH | $e^{-[\text{SSq}]m} + (1-e^{-[\text{SSq}]m}) = 1$ | Complementary partition identity |

These three number systems are not ad-hoc constructions but arise directly from the
UQFF vacuum structure: VDS from layer-by-layer SCm compression, DVP from prime
quantisation of vortex circulation, and BH from harmonic decomposition of the
charge-reactivity buoyancy field. Together they form the complete UQFF vacuum
number algebra.

---

**References:**  
PAPER_535 (VDS/DVP/BH catalogue hub) | PAPER_543 (Navier-Stokes regularity proof) |
PAPER_536 (DPM split monopole proplyd topology) | PAPER_652 (fine structure constant) |
PAPER_1127 (SCm LQG holonomy) | PAPER_1128 (SCm String Theory 26D) |
COMPLETE_UQFF_EQUATIONS_REFERENCE.md
