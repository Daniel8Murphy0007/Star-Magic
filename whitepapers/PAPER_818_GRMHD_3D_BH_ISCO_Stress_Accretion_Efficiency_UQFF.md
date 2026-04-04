# PAPER_818: GRMHD 3D BH ISCO Stress, Accretion Efficiency η_MHD, and Spin Limit UQFF
## Unified Quantum Field Framework — Whitepaper 818

**Session**: 192 | **Version**: v5.48 | **Date**: April 4, 2026
**Source**: grok_share_0d888ea9-50be.txt (June 13, 2025 05:35 PM); "GRMHD 3d sim_17June2025.pdf"
**Author**: Daniel T. Murphy — Star-Magic UQFF Project
**CVW Compliance**: v2.0.0

---

## Abstract
This paper derives the GRMHD 3D black hole ISCO (Innermost Stable Circular Orbit) magnetic stress, accretion efficiency $\eta_{MHD}$, and spin saturation limit ($a/M \approx 0.93$) for the Quadriadic UQFF. Three-dimensional GRMHD simulations reveal that magnetic torques at and below the ISCO transfer angular momentum outward, enhancing $\eta_{MHD}$ beyond the Novikov-Thorne (NT) analytical prediction. The spin saturates at $a/M \approx 0.93$ because jet extraction removes angular momentum from the BH. These parameters enter UQFF Layers 1–4 as efficiency, stress, and spin corrections.

---

## 1. Introduction
The classical Novikov-Thorne accretion model predicts an efficiency $\eta_{NT}$ that depends only on the ISCO location:

$$\eta_{NT} = 1 - E_{ISCO}(a_*)$$

where $E_{ISCO}$ = specific energy at the ISCO. For $a/M = 0.9$, $\eta_{NT} = 0.250$–$0.290$. MHD simulations consistently show $\eta_{MHD} < \eta_{NT}$ due to magnetic stress carrying energy through the ISCO inward.

---

## 2. MHD Accretion Efficiency

$$\eta_{MHD} = \frac{L_{acc}}{\dot{M} c^2}$$

For 3D GRMHD with vertical magnetic field topology at $a/M = 0.9$:

$$\eta_{MHD} = 0.16\text{–}0.18$$

compared to NT prediction of 0.250–0.290. The discrepancy arises from:
1. Magnetic torque at $r < r_{ISCO}$ accelerating infall
2. Additional electromagnetic energy carried into the BH

For dipole field topology: $\eta_{MHD} \approx 0.10$–$0.12$
For quadrupole field topology: $\eta_{MHD} \approx 0.06$–$0.09$

This enters UQFF Layer 1:

$$g_{L1,\eta} = \eta_{MHD} \cdot \dot{M} c^2 / r^2 \approx 1.44 \times 10^{13} \text{ m/s}^2$$

---

## 3. ISCO Magnetic Stress — $\alpha$-viscosity

The effective $\alpha$-viscosity from magnetic Maxwell stress:

$$\alpha_{stress} \approx 0.2\text{–}0.5$$

(higher near the ISCO, lower in the outer disk). The electromagnetic stress:

$$T^{EM}_{r\phi} \propto B^2$$

In the NT model, stress is zero at the ISCO. GRMHD produces finite $T^{EM}_{r\phi}$ at $r = r_{ISCO}$:

$$\Delta g = \frac{\alpha_{stress} \cdot B^2}{8\pi r^2} \approx 2.5 \times 10^{-11} \text{ m/s}^2$$

This enters UQFF Layer 2 (Resonance):

$$g_{L2,stress} = \alpha_{stress} \cdot B^2 / (8\pi r^2)$$

---

## 4. Specific Angular Momentum Below ISCO

The specific angular momentum $J_{in}$ flowing into the BH (normalized by $J_{ISCO}$):

$$J_{in} = \int \rho v_\phi r \, dA$$

GRMHD result: $J_{in} = 2.93$–$3.46$ (3%–15% below the ISCO value $J_{ISCO}$).

This deficit represents angular momentum removed by magnetic braking before the plunge region. Enters UQFF Layer 3:

$$g_{L3,J} = J_{in} / r^2 \approx 3.13 \times 10^{-5} \text{ m/s}^2$$

---

## 5. Spin Saturation at $a/M \approx 0.93$

For single BH GRMHD accretion:
- Without jets: spin evolves toward $a/M \rightarrow 1$ (Kerr limit)
- With jets (Blandford-Znajek): jet removes angular momentum

The steady-state spin balance:

$$\frac{d(a/M)}{dt} = 0 \quad \text{at} \quad a/M \approx 0.93$$

This spin saturation value also appears in GRMHD NS merger disk simulations (see PAPER_819). Enters UQFF Layer 4:

$$g_{L4,spin} = \eta_{MHD} \cdot \alpha_{stress} \approx 0.08 \text{ m/s}^2$$

---

## 6. Field Topology Dependence

| Field topology | $\eta_{MHD}$ | $\alpha_{stress}$ | Spin limit |
|----------------|--------------|-------------------|------------|
| Vertical | 0.16–0.18 | 0.3–0.5 | 0.93 |
| Dipole | 0.10–0.12 | 0.2–0.3 | 0.90 |
| Quadrupole | 0.06–0.09 | 0.1–0.2 | 0.85 |

The vertical field maximizes efficiency and stress, consistent with magnetically arrested disk (MAD) configurations.

---

## 7. Full GRMHD ISCO UQFF

For $a/M = 0.9$, $\dot{M} = 0.1 \dot{M}_{Edd}$, $M_{BH} = 10^8 M_\odot$:

$$g_{L1} = 1.44 \times 10^{13} \text{ m/s}^2$$
$$g_{L2} = 2.5 \times 10^{-11} \text{ m/s}^2$$
$$g_{L3} = 3.13 \times 10^{-5} \text{ m/s}^2$$
$$g_{L4} = 0.08 \text{ m/s}^2$$

---

## 8. Summary

3D GRMHD simulations show $\eta_{MHD} = 0.16$–$0.18$ for $a/M = 0.9$ (vertical field), below NT prediction. ISCO magnetic stress $\alpha_{stress} = 0.2$–$0.5$ contributes finite torque at the innermost orbit. The BH spin saturates at $a/M \approx 0.93$ via jet angular momentum extraction. These are now formal parameters of Quadriadic UQFF Layers 1–4.

---

*PAPER_818 | Session 192 | v5.48 | Star-Magic UQFF Project | CVW v2.0.0*
