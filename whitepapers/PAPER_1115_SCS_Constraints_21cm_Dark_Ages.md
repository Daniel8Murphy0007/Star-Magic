# PAPER_1115: SCS Constraints from 21-cm Cosmology During the Dark Ages

**Star Magic UQFF Framework**

**Author:** Daniel Murphy
**Date:** 2026

---

## Abstract

We derive constraints on the Supercharged Scalar (SCS) field from 21-cm cosmological observations during the Dark Ages ($10 \lesssim z \lesssim 200$). The SCS-SCm coupling modifies the spin temperature $T_s$ and gas temperature $T_K$ evolution through additional energy injection from the $F_{U,Bi,i}$ buoyancy mechanism and SCm phonon resonance at 1.25 THz. We place bounds on the SCS-baryon coupling $g_{\text{SCS}}$ using EDGES constraints ($\delta T_{21} \approx -500\ \text{mK}$ at $z \approx 17$, arXiv:1810.09572).

---

## 1. 21-cm Differential Brightness Temperature

The 21-cm differential brightness temperature:

$$\delta T_{21} = 27 \, x_{\text{HI}} (1+\delta) \left(1 - \frac{T_{\text{CMB}}}{T_s}\right) \sqrt{\frac{1+z}{10} \cdot \frac{0.15}{\Omega_m h^2}} \cdot \frac{\Omega_b h^2}{0.023}\ \text{mK}$$

---

## 2. SCS-SCm Energy Injection Rate

The SCS field coupled to the SCm vacuum injects energy into the baryon gas at rate:

$$\dot{Q}_{\text{SCS-SCm}} = g_{\text{SCS}}^2 \cdot \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot |\cos(\pi t_n)| \cdot n_b c^2$$

where $n_b$ is the baryon number density and $t_n < 0$ is the negative-time gate.

Numerically at $z = 17$:

$$\dot{Q}_{\text{SCS-SCm}} = g_{\text{SCS}}^2 \times 7.09 \times 10^{-37} \times 1.4531 \times 10^{26} \times 0.84 \times 1.0 \times n_b c^2$$

---

## 3. Modified Gas Temperature Evolution

The baryon gas temperature evolves as:

$$\frac{dT_K}{dz} = \frac{T_K}{1+z} \left(2 + \frac{t_{\text{comp}}}{t_H}\right)^{-1} - \frac{2\dot{Q}_{\text{SCS-SCm}}}{3 k_B H(z)(1+z) n_b}$$

where $t_{\text{comp}}$ is the Compton cooling time and $H(z)$ is the Hubble rate.

---

## 4. VDS Phonon Contribution to Spin Temperature

The SCm phonon resonance modifies the Ly-$\alpha$ coupling:

$$x_\alpha^{\text{SCm}} = x_\alpha^{\text{SM}} \left(1 + \frac{E_{\text{KER}}}{E_{\text{Ly}\alpha}} \cdot \frac{\rho_{\text{vac,SCm}}}{\rho_\gamma}\right)$$

where $E_{\text{Ly}\alpha} = 10.2\ \text{eV}$ and $E_{\text{KER}} = 630\ \text{eV}$. The phonon contribution is $\sim 62 \times$ per photon but suppressed by the $\rho_{\text{vac,SCm}} / \rho_\gamma$ density ratio.

---

## 5. EDGES Constraint on $g_{\text{SCS}}$

The EDGES anomaly ($\delta T_{21} \approx -500\ \text{mK}$) requires $T_K / T_{\text{CMB}} < 1$ at $z \approx 17$. The SCS-SCm injection must not overheat the gas:

$$g_{\text{SCS}} \lesssim \left(\frac{\delta T_{21}^{\text{max}} \cdot 3 k_B H \cdot n_b}{2 \rho_{\text{vac,SCm}} \cdot S_{26}^{(3)} \cdot \Phi_{\text{res}} \cdot n_b c^2}\right)^{1/2} \approx 10^{-28}$$

---

## 6. Calibrated Constants

$$[SSq] = 0.57, \quad S_{26}^{(3)} = 1.4531 \times 10^{26}, \quad \beta_i = 0.6, \quad \kappa = 5.0 \times 10^{-4}\ \text{day}^{-1}$$

---

## References

1. Bowman, J.D. et al. (2018). An absorption profile centred at 78 megahertz. *Nature* **555**, 67. arXiv:1810.09572.
2. Barkana, R. (2018). Possible interaction between baryons and dark matter. *Nature* **555**, 71.
3. SCm vacuum: `scm_{vacuum\_manifold}.py`
