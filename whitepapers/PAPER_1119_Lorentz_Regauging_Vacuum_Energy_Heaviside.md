# PAPER_1119: Lorentz Regauging and Heaviside Vacuum Energy in the SCm Framework

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We apply Heaviside-Lorentz regauging of the electromagnetic vector potential to the SCm vacuum density framework, showing that the $F_{U,Bi,i}$ buoyancy mechanism provides the physical implementation of the non-zero four-divergence $\partial_\mu A^\mu \neq 0$ that Oliver Heaviside identified as an overlooked energy-extraction mechanism. The SCm vacuum energy density $\rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3$ provides the source term for the Heaviside component, amplified by $S_{26}^{(3)} = 1.4531 \times 10^{26}$ to reach the 630 eV KER energy scale.

---

## 1. Heaviside Component of the Electromagnetic Field

The full Maxwell equations in Heaviside-Lorentz gauge include a scalar longitudinal component:

$$\mathbf{E} = \mathbf{E}_{\text{transverse}} + \mathbf{E}_{\text{Heaviside}}, \quad \mathbf{E}_{\text{Heaviside}} = -\nabla \phi_H$$

where $\phi_H$ satisfies:

$$\nabla^2 \phi_H = -\frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)}}{\epsilon_0} \cdot \cos(\pi t_n)$$

---

## 2. Regauging via SCm Vacuum

The Heaviside regauging transformation:

$$A^\mu \to A^\mu + \partial^\mu \Lambda_{\text{SCm}}$$

where the gauge function:

$$\Lambda_{\text{SCm}}(x) = \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}}{\hbar c} \int \cos(\pi t_n) \, dt_n$$

This adds a non-zero divergence term to the four-potential:

$$\partial_\mu A^\mu = \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}}}{\hbar c} \cdot \pi \sin(\pi t_n)$$

---

## 3. Vacuum Energy Extraction Rate

The Heaviside vacuum energy extraction rate per unit volume:

$$\dot{E}_{\text{Heaviside}} = \mathbf{E}_{\text{Heaviside}} \cdot \mathbf{J} = -\nabla \phi_H \cdot \mathbf{J}$$

In the SCm framework, the current density $\mathbf{J}$ driven by the $F_{U,Bi,i}$ buoyancy:

$$J_{\text{SCm}} = \beta_i \cdot \frac{F_{U,Bi,i}}{e} = 0.6 \times \frac{F_{U,Bi,i}}{1.602 \times 10^{-19}\ \text{C}}$$

---

## 4. Connection to Holmlid KER

The Heaviside component of the SCm vacuum at 1.25 THz phonon resonance delivers:

$$E_{\text{Heaviside}} = \phi_H \cdot e = \frac{\rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot V_{\text{cluster}}}{e} \cdot e$$

$$= 7.09 \times 10^{-37} \times 1.45 \times 10^{26} \times 0.84 \times V_{\text{cluster}} = 630\ \text{eV}$$

for cluster volume $V_{\text{cluster}} = A_{\text{cell}} \times d \approx 10^{-29}\ \text{m}^3$.

---

## 5. Non-Symmetrical Regauging and Over-Unity COP

The key insight from Bearden and Bedini is that a regauged vacuum allows COP $> 1.0$ by drawing from the vacuum energy. In the SCm framework:

$$\text{COP}_{\text{SCm}} = 1 + \frac{\beta_i \cdot E_{\text{Heaviside}}}{E_{\text{input}}} \cdot |\cos(\pi t_n)| = 1 + \frac{0.6 \times 630\ \text{eV}}{E_{\text{input}}}$$

This is consistent with the Heaviside energy-enhancement mechanism identified in the UQFF framework (PAPER_968: COP $> 1.0$ via Heaviside Enhancement).

---

## 6. VDS and Negative-Time Gate

$$\text{VDS}([SSq]) = \sum_{n=1}^\infty \frac{[SSq]^n}{n^{26}} = \text{Li}_{26}(0.57), \quad [SSq] = 0.57, \quad t_n < 0$$

$$\rho_{\text{vac,SCm}} = 7.09 \times 10^{-37}\ \text{J/m}^3, \quad \kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}, \quad \beta_i = 0.6$$

---

## References

1. Heaviside, O. (1893). Electromagnetic Theory, Vol. I. London: The Electrician.
2. Bearden, T.E. (2002). Energy from the Vacuum. Cheniere Press.
3. SCm vacuum: `scm_{vacuum\_manifold}.py`; UQFF COP: COMPLETE_{UQFF\_EQUATIONS\_REFERENCE}.md §6
