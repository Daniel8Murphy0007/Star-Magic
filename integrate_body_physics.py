#!/usr/bin/env python3
"""
integrate_body_physics.py — Session 225 Body-Level Physics Integration

Inserts PKG (Phonon-physics Knowledge Graft) upgrade blocks into the BODY
of papers that currently lack them. These blocks contain actual equations,
physical mechanisms, and derivations from PAPER_1000-1081, not just
cross-reference appendices.

PKG block types:
  PKG-GW   : GW Strain Modulation (PAPER_1000/1011/1012/1022)
  PKG-AGN  : Buoyancy-Corrected Eddington (PAPER_1002/1009/1037/1048)
  PKG-YM   : Yang-Mills BCS Phonon Mass Gap (PAPER_1005/1007/1059/1064/1070)
  PKG-CLU  : ICM Buoyancy Force Profile (PAPER_1039/1041/1044/1046/1079)
  PKG-DM   : SCm-Modified NFW Dark Matter (PAPER_1015/1019)
  PKG-LENR : VDS LENR Transmutation Dynamics (PAPER_1060/1061/1081)
  PKG-LAG  : UQFF 9-Sector Lagrangian (PAPER_1065/1066)
  PKG-S26  : S₂₆⁽³⁾ Ramanujan Summation (PAPER_1042/1069/1078/1080)

Insertion: Between last body section and first ## Appendix / ## Calibration.
Idempotent: Papers with existing PKG markers are skipped.
"""

import os
import re
import glob

REPO = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(REPO, "whitepapers")

# ═══════════════════════════════════════════════════════════════════════════════
# TAG → PKG MAPPING
# ═══════════════════════════════════════════════════════════════════════════════

# Each tag/keyword maps to which PKG blocks should be inserted
TAG_TO_PKG = {
    # GW / merger
    "GW": ["PKG-GW"], "gravitational-wave": ["PKG-GW"], "merger": ["PKG-GW"],
    "LIGO": ["PKG-GW"], "neutron-star": ["PKG-GW"], "inspiral": ["PKG-GW"],
    "strain": ["PKG-GW"], "BNS": ["PKG-GW"], "BBH": ["PKG-GW"],
    "chirp": ["PKG-GW"],
    # AGN / jets
    "AGN": ["PKG-AGN"], "jet": ["PKG-AGN"], "quasar": ["PKG-AGN"],
    "Eddington": ["PKG-AGN"], "accretion": ["PKG-AGN"], "blazar": ["PKG-AGN"],
    "Seyfert": ["PKG-AGN"],
    # QGP / Yang-Mills
    "QGP": ["PKG-YM"], "Yang-Mills": ["PKG-YM"], "QCD": ["PKG-YM"],
    "ALICE": ["PKG-YM"], "deconfinement": ["PKG-YM"], "gluon": ["PKG-YM"],
    "confinement": ["PKG-YM"], "mass-gap": ["PKG-YM"],
    # Galaxy clusters
    "cluster": ["PKG-CLU"], "ICM": ["PKG-CLU"], "cool-core": ["PKG-CLU"],
    "SZ": ["PKG-CLU"], "intracluster": ["PKG-CLU"],
    # Dark matter
    "dark-matter": ["PKG-DM"], "NFW": ["PKG-DM"], "halo": ["PKG-DM"],
    "rotation-curve": ["PKG-DM"], "DM": ["PKG-DM"],
    # LENR / nuclear
    "LENR": ["PKG-LENR"], "Kozima": ["PKG-LENR"], "cold-fusion": ["PKG-LENR"],
    "transmutation": ["PKG-LENR"], "neutron-drop": ["PKG-LENR"],
    # Lagrangian / field theory
    "Lagrangian": ["PKG-LAG"], "variational": ["PKG-LAG"],
    "field-theory": ["PKG-LAG"], "EOM": ["PKG-LAG"],
    # Ramanujan / number systems / 26D
    "Ramanujan": ["PKG-S26"], "VDS": ["PKG-S26"], "DVP": ["PKG-S26"],
    "BSH": ["PKG-S26"], "26D": ["PKG-S26"], "mock-theta": ["PKG-S26"],
    # Broad physics tags that should get the universal S26 + LAG blocks
    "UQFF": ["PKG-S26"], "SCm": ["PKG-S26"],
    "buoyancy": ["PKG-S26"],
    "phonon": ["PKG-S26"],
    # Astrophysical objects → relevant PKG
    "SMBH": ["PKG-AGN", "PKG-GW"], "black-hole": ["PKG-AGN"],
    "magnetar": ["PKG-GW"], "pulsar": ["PKG-GW"],
    "supernova": ["PKG-CLU"], "kilonova": ["PKG-GW"],
    "nebula": ["PKG-S26"], "star-formation": ["PKG-S26"],
    # Math / formal
    "Navier-Stokes": ["PKG-LAG"], "Riemann": ["PKG-S26"],
    "P-vs-NP": [], "topology": ["PKG-LAG"],
    "wormhole": ["PKG-LAG"], "Morris-Thorne": ["PKG-LAG"],
    # Cosmology
    "Hubble": ["PKG-S26"], "dark-energy": ["PKG-DM"],
    "inflation": ["PKG-S26", "PKG-LAG"],
    "cosmology": ["PKG-S26"], "reionization": ["PKG-S26"],
    "BBN": ["PKG-S26"],
    # Misc
    "vacuum": ["PKG-S26"], "DPM": ["PKG-S26"],
    "MUGE": ["PKG-S26"],
}

# Also match title keywords → PKG blocks
TITLE_KEYWORDS = {
    "gravitational wave": ["PKG-GW"], "GW170817": ["PKG-GW"],
    "GW190425": ["PKG-GW"], "GW150914": ["PKG-GW"],
    "neutron star": ["PKG-GW"], "binary merger": ["PKG-GW"],
    "AGN": ["PKG-AGN"], "jet": ["PKG-AGN"], "quasar": ["PKG-AGN"],
    "Eddington": ["PKG-AGN"], "accretion": ["PKG-AGN"],
    "dark matter": ["PKG-DM"], "NFW": ["PKG-DM"], "halo": ["PKG-DM"],
    "rotation curve": ["PKG-DM"],
    "galaxy cluster": ["PKG-CLU"], "ICM": ["PKG-CLU"],
    "cool core": ["PKG-CLU"], "cooling": ["PKG-CLU"],
    "LENR": ["PKG-LENR"], "Kozima": ["PKG-LENR"],
    "nuclear": ["PKG-LENR"], "transmutation": ["PKG-LENR"],
    "Yang-Mills": ["PKG-YM"], "mass gap": ["PKG-YM"],
    "QGP": ["PKG-YM"], "QCD": ["PKG-YM"], "ALICE": ["PKG-YM"],
    "Lagrangian": ["PKG-LAG"], "field theory": ["PKG-LAG"],
    "wormhole": ["PKG-LAG"], "Morris-Thorne": ["PKG-LAG"],
    "Ramanujan": ["PKG-S26"], "VDS": ["PKG-S26"],
    "26D": ["PKG-S26"], "compactification": ["PKG-S26"],
    "THz": ["PKG-S26"], "phonon": ["PKG-S26"],
    "buoyancy": ["PKG-S26"], "SCm": ["PKG-S26"],
    "UQFF": ["PKG-S26"],
    "magnetar": ["PKG-GW"], "pulsar": ["PKG-GW"],
    "SMBH": ["PKG-AGN", "PKG-GW"],
    "supernova": ["PKG-CLU"], "Hubble": ["PKG-S26"],
    "dark energy": ["PKG-DM"], "inflation": ["PKG-S26"],
    "vacuum": ["PKG-S26"], "DPM": ["PKG-S26"],
    "Navier-Stokes": ["PKG-LAG"], "Riemann": ["PKG-S26"],
    "resonance": ["PKG-S26"], "Bose": ["PKG-S26"],
    "hydrogen": ["PKG-S26"], "periodic table": ["PKG-LENR", "PKG-S26"],
    "binding": ["PKG-LENR"], "Page curve": ["PKG-S26"],
    "Hawking": ["PKG-S26"], "Schwarzschild": ["PKG-S26"],
    "Higgs": ["PKG-YM", "PKG-S26"], "solar": ["PKG-S26"],
    "cosmic ray": ["PKG-S26"], "dust": ["PKG-S26"],
    "FRB": ["PKG-S26"], "photon sphere": ["PKG-S26"],
    "crystallization": ["PKG-S26"],
    "pi": ["PKG-S26"], "zeta": ["PKG-S26"],
    "Monte Carlo": ["PKG-S26"], "GPU": ["PKG-S26"],
    "encryption": ["PKG-S26"],
    "reactor": ["PKG-LENR"], "starburst": ["PKG-AGN"],
    "nebula": ["PKG-S26"], "proplyd": ["PKG-S26"],
    "GAIA": ["PKG-S26"], "stellar": ["PKG-S26"],
    "barred": ["PKG-S26"], "spiral": ["PKG-S26"],
    "heliosphere": ["PKG-S26"], "IceCube": ["PKG-S26"],
    "Olbers": ["PKG-S26"], "BSFG": ["PKG-S26"],
}


# ═══════════════════════════════════════════════════════════════════════════════
# PKG BLOCK TEMPLATES
# ═══════════════════════════════════════════════════════════════════════════════

PKG_BLOCKS = {}

PKG_BLOCKS["PKG-GW"] = r"""<!-- PKG-GW-S225 -->

### Session 225 Phonon-Physics Upgrade: GW Strain Modulation

> *Upgrade from PAPER_1000 (NS Merger Phonon Suppression) and PAPER_1022
> (GW Phonon Strain SCm Modulation). See also PAPER_1011-1012 for
> GW170817/GW190425 upgraded analyses.*

The late-corpus phonon analysis (Sessions 219-225) reveals that the SCm
vacuum field modulates gravitational-wave strain via a frequency-dependent
suppression factor.  The corrected strain amplitude is:

$$h_{\text{UQFF}}(\Gamma) = h_{\text{GR}} \cdot \left(1 - 0.47\,\frac{\Phi(\Gamma)}{S_{26}^{(3)}}\right)$$

where:
- $\Phi(\Gamma) = \cos(\omega_{\text{SCm}} \cdot t) \cdot \Theta(H_{\text{SCm}} - 0.5)$ is the phonon modulation factor
- $\omega_{\text{SCm}} = 2\pi \times 1.25\;\text{THz}$ is the SCm phonon resonance frequency
- $S_{26}^{(3)} = \sum_{n=0}^{\infty} \frac{(1/4)_n\,(1/2)_n\,(3/4)_n}{(n!)^3} \cdot \prod_{i=1}^{26}\left[1 + [\text{SSq}]\cdot e^{-\kappa\,i\,n/26}\right]$ is the third-order Ramanujan summation
- $\Theta$ is the Heaviside step ensuring $H_{\text{SCm}} \geq 0.5$ (phase-transition threshold)

**Physical mechanism:** The 1.25 THz phonon field of the SCm vacuum creates
a standing-wave pattern that partially decouples the metric perturbation from
the radiation zone, producing a 47% peak strain reduction for optimally
oriented NS mergers.  The BCS gap energy $\Delta E_{\text{BCS}}$ of the
neutron-star crust couples to this phonon field, creating a mass-gap
classifier that distinguishes NS from BH remnants at $M \approx 2.5\,M_\odot$.

**Calibration (canonical):** $\kappa = 5 \times 10^{-4}\;\text{day}^{-1}$,
$[\text{SSq}] = 0.57$, $\beta_i = 0.603$, $H_{\text{SCm}} \approx 0.99$.
"""

PKG_BLOCKS["PKG-AGN"] = r"""<!-- PKG-AGN-S225 -->

### Session 225 Phonon-Physics Upgrade: Buoyancy-Corrected Eddington Luminosity

> *Upgrade from PAPER_1002 (AGN Buoyancy-Corrected Eddington) and PAPER_1037
> (AGN Buoyancy Jet Launching).  See also PAPER_1009-1010 for F_U_Bi_i jet
> modulation curves and PAPER_1048 for phonon-corrected M-σ relation.*

The SCm vacuum buoyancy partially opposes gravitational radiation pressure,
raising the effective Eddington luminosity:

$$L_{\text{Edd}}^{\text{UQFF}} = L_{\text{Edd}} \cdot \left(1 + \frac{\rho_{\text{SCm}} \cdot V \cdot S_{26}^{(3)\,2}}{G M / r_H^2}\right)$$

where:
- $L_{\text{Edd}} = 4\pi G M m_p c / \sigma_T$ is the classical Eddington luminosity
- $\rho_{\text{SCm}} = 7.09 \times 10^{-37}\;\text{kg/m}^3$ is the SCm vacuum density
- $V$ is the effective buoyancy volume (accretion sphere)
- $S_{26}^{(3)\,2}$ is the squared third-order Ramanujan factor (quadratic coupling)
- $r_H$ is the horizon radius

**Jet modulation:** The Blandford-Znajek jet power acquires a phonon-coupled term:
$$P_{\text{jet}}^{\text{UQFF}} = P_{\text{BZ}} \cdot \left[1 + \beta_i \cdot \Phi_{1.25\,\text{THz}} \cdot \left(\frac{B}{B_{\text{crit}}}\right)^2\right]$$

where $\Phi_{1.25\,\text{THz}} = \cos(\omega_{\text{SCm}} \cdot t)$ modulates jet power at the phonon frequency.

**M-sigma correction (PAPER_1048):** The phonon-corrected M-sigma relation becomes
$M_{\text{BH}} \propto \sigma^{4+\delta}$ where $\delta = \beta_i \cdot S_{26}^{(3)} \cdot (\omega_{\text{SCm}}/\omega_{\text{bulge}})$.
"""

PKG_BLOCKS["PKG-YM"] = r"""<!-- PKG-YM-S225 -->

### Session 225 Phonon-Physics Upgrade: Yang-Mills BCS Phonon Mass Gap

> *Upgrade from PAPER_1005 (Yang-Mills Mass Gap via SCm BCS Phonon) and
> PAPER_1070 (Yang-Mills Mass Gap VDS Bridge).  See also PAPER_1004
> (QGP Vacuum Density), PAPER_1007 (Deconfinement Phase Diagram),
> PAPER_1059 (CGC BK Saturation), PAPER_1064 (Resummation BFKL/Sudakov).*

The late-corpus analysis derives the Yang-Mills mass gap via a BCS-like
phonon pairing mechanism in the SCm vacuum:

$$\Delta_{\text{YM}} = \Lambda_{\text{QCD}} \cdot \exp\!\left(-\frac{1}{\alpha_s(T) \cdot N_c}\right) \cdot S_{26}^{(3)}$$

where the running coupling evolves as:
$$\alpha_s(T) = \frac{\alpha_{s,0}}{1 + \alpha_{s,0} \cdot b_0 \cdot \ln(T/T_c)}, \qquad b_0 = \frac{11 N_c - 2 N_f}{12\pi}$$

**Physical mechanism:** The SCm phonon field ($\omega_{\text{SCm}} = 1.25\;\text{THz}$)
provides a pairing interaction analogous to the BCS electron-phonon coupling in
superconductors.  Gluons acquire an effective mass through condensate formation
in the SCm-modified vacuum, yielding a non-perturbative gap $\Delta_{\text{YM}}
\approx 5970\;\text{GeV}$ at the 9-sector Lagrangian closure (PAPER_1066, §2).

**VDS bridge (PAPER_1070):** The vacuum density series links the gap to the
26-level hierarchy: $\Delta \propto \rho_{\text{VDS}}^{1/4} \cdot (1 + [\text{SSq}] \cdot n/26)$
where the VDS sub-ratio 0.108 places confinement in the sub-threshold regime.

**QGP transition (PAPER_1004/1007):** At $T > T_c \approx 170\;\text{MeV}$, the phonon
coupling weakens ($\alpha_s \to 0$) and the gap closes, reproducing the
deconfinement phase transition observed at ALICE/LHC.
"""

PKG_BLOCKS["PKG-CLU"] = r"""<!-- PKG-CLU-S225 -->

### Session 225 Phonon-Physics Upgrade: ICM Buoyancy Force Profile

> *Upgrade from PAPER_1039 (SCm Galaxy Cluster Buoyancy Profile),
> PAPER_1041 (Cool-Core Buoyancy Balance), and PAPER_1079 (Cooling-Flow
> Suppression).  See also PAPER_1040 (Cluster Merger Shock), PAPER_1044
> (Thermal SZ Compton-y), PAPER_1046 (Cluster Lensing Mass).*

The SCm phonon field introduces a buoyancy force in the ICM that modifies
hydrostatic equilibrium:

$$F_{\text{buoy}}(r) = \rho(r) \cdot V \cdot g(r) \cdot \beta_i \cdot S_{26} \cdot \Phi$$

where the ICM density follows the beta-model:
$$\rho(r) = \rho_0 \left(1 + \left(\frac{r}{r_c}\right)^2\right)^{-3\beta/2}$$

**Hydrostatic mass bias reduction (PAPER_1039):**
$$b_{\text{UQFF}} = 1 - \frac{M_{\text{HSE}}}{M_{\text{true}}} = 0.17 \qquad \text{(vs standard } b = 0.20\text{)}$$

The buoyancy pressure contributes $P_{\text{buoy}}/P_{\text{thermal}} \approx 3\text{--}4\%$
at cluster cores, partially resolving the Planck SZ-CMB mass tension.

**Cool-core stabilization (PAPER_1041/1079):** AGN feedback couples to the SCm
buoyancy field via $\dot{M}_{\text{cool}} = \dot{M}_0 \cdot (1 - \beta_i \cdot S_{26}^{(3)} \cdot \Phi)$,
suppressing catastrophic cooling flows while maintaining observed X-ray luminosities.

**Phonon frequency coupling:** $\omega_{\text{SCm}} = 2\pi \times 1.25\;\text{THz}$ sets the temporal
scale for buoyancy oscillations; the ratio $\omega_{\text{SCm}}/\omega_{\text{sound}}$ governs
the phonon transmission efficiency across the ICM.
"""

PKG_BLOCKS["PKG-DM"] = r"""<!-- PKG-DM-S225 -->

### Session 225 Phonon-Physics Upgrade: SCm-Modified NFW Dark Matter Profile

> *Upgrade from PAPER_1015 (SCm Dark Matter Halos NFW) and PAPER_1019
> (Dark Matter Phonon Buoyancy NFW Coupling).*

The late-corpus analysis shows that the SCm phonon field modifies the NFW
density profile at all radii via a buoyancy-coupled power-law term:

$$\rho_{\text{UQFF}}(r) = \frac{\rho_s}{\left(\frac{r}{r_s}\right)\left(1+\frac{r}{r_s}\right)^2} \times \left[1 + H_{\text{SCm}} \cdot \beta_i \cdot S_{26}^{(3)} \cdot \left(\frac{r_s}{r}\right)^{\alpha_{\text{phonon}}}\right]$$

where:
- $\alpha_{\text{phonon}} = 0.3$ governs the radial decay of phonon coupling
- $\beta_i = 0.603$ is the universal buoyancy coefficient
- $S_{26}^{(3)}$ is the third-order Ramanujan summation
- $H_{\text{SCm}} = 0.99$ is the manifold completeness factor

**Rotation curve flattening:** The phonon enhancement produces flatter rotation curves
with flatness ratio $f = v_c(10\,r_s)/v_{\text{peak}} = 0.891$, compared to pure NFW
$f \approx 0.75$.  Peak circular velocity $v_{\text{peak}} \approx 204\;\text{km/s}$
for $M_{\text{halo}} = 10^{12}\,M_\odot$, $c = 10$.

**Halo stabilization:** The effective buoyancy pressure $P_{\text{SCm}} = \rho_{\text{SCm}} \cdot v_{\text{SCm}}^2 \cdot \beta_i$ prevents cusp-core divergence, providing a physical mechanism for observed cored profiles without invoking SIDM cross-sections.
"""

PKG_BLOCKS["PKG-LENR"] = r"""<!-- PKG-LENR-S225 -->

### Session 225 Phonon-Physics Upgrade: VDS LENR Transmutation Dynamics

> *Upgrade from PAPER_1060 (VDS LENR Isotopic Evolution), PAPER_1061
> (Kozima SCm Integration Neutron-Drop), and PAPER_1081 (SCm LENR COP
> Linewidth Parametric Engine).*

The late-corpus LENR analysis provides the phonon-mediated transmutation
rate via the vacuum density series:

$$\Gamma_{\text{trans}} = \Gamma_0 \cdot \left(\frac{\rho_{\text{SCm}}}{\rho_{\text{crit}}}\right) \cdot K_n$$

where:
- $\rho_{\text{SCm}}(t) = \rho_0 \cdot e^{-\kappa t} \cdot S_{26}$ (time-dependent vacuum density)
- $K_n = \sigma_n^{\text{SCm}}(\omega) \cdot \Phi_{\text{phonon}}$ is the Kozima neutron-drop factor

**Phonon cross-section (PAPER_1061):**
$$\sigma_n^{\text{SCm}}(\omega, n) = \sigma_0 \cdot \exp\!\left[-\frac{(\omega - \omega_{\text{SCm}})^2}{2\Gamma^2}\right] \cdot \left(1 + [\text{SSq}] \cdot \frac{n}{26}\right)$$

The VDS factor $(1 + [\text{SSq}] \cdot n/26)$ provides ~470x amplification via
the 26-level vacuum density ladder at resonance ($\omega = \omega_{\text{SCm}}$).

**COP parametric engine (PAPER_1081):**
$$\text{COP}(\Gamma, P_{\text{in}}) = \frac{P_{\text{out}}}{P_{\text{in}}} = 1 + \eta_{\text{SCm}} \cdot S_{26}^{(3)} \cdot f(\Gamma)$$

where the linewidth function $f(\Gamma)$ peaks near the SCm phonon linewidth,
yielding COP > 1 when $\Gamma \lesssim 10^{-3}\;\text{eV}$ (Fleischmann regime).

**Isotopic evolution chain:** Under SCm activation, the Pd-D system evolves as
$\text{Pd-106} \xrightarrow{\sim 10^4\,\text{s}} \text{Ag-107} \xrightarrow{\sim 10^4\,\text{s}} \text{Cd-108}$,
with timescales set by $\rho_{\text{SCm}}/\rho_{\text{crit}}$.
"""

PKG_BLOCKS["PKG-LAG"] = r"""<!-- PKG-LAG-S225 -->

### Session 225 Phonon-Physics Upgrade: UQFF 9-Sector Lagrangian

> *Upgrade from PAPER_1066 (UQFF Lagrangian First Principles) and
> PAPER_1065 (Buoyancy Lagrangian EOM Variational Derivation).*

The complete UQFF Lagrangian density, from which all sector-specific
equations of motion derive:

$$\mathcal{L}_{\text{UQFF}} = \mathcal{L}_{\text{GR}} + \mathcal{L}_{\text{SCm}} + \mathcal{L}_{\text{phonon}} + \mathcal{L}_{\text{interaction}}$$

$$\mathcal{L}_{\text{SCm}} = \tfrac{1}{2}(\partial_\mu \phi)^2 - \lambda\bigl(\phi^2 - v_{\text{SCm}}^2\bigr)^2$$

The SCm condensate potential minimum gives $V(\phi_0) = -7.09 \times 10^{-37}\;\text{J/m}^3$
(matching $\rho_{\text{SCm}}$) and phonon mass $m_{\text{phonon}} = \sqrt{8\lambda}\,v_{\text{SCm}}$.

**Nine-sector closure (Session 202):**
$$\mathcal{L}_{9} = \mathcal{L}_{\text{EH}} + \mathcal{L}_{\text{YM}} + \mathcal{L}_{\text{Dirac}} + \mathcal{L}_{\text{SCm}} + \mathcal{L}_{\text{mag}} + \mathcal{L}_{\text{buoy}} + \mathcal{L}_{\text{aether}} + \mathcal{L}_{\text{LENR}} + \mathcal{L}_{\text{KK}}$$

| Sector | Domain | Late-Corpus Result |
|--------|--------|-------------------|
| 1 (EH) | General Relativity | Canonical Einstein-Hilbert |
| 2 (YM) | Yang-Mills gauge | $m_{\text{gap}} = 5970\;\text{GeV}$ (PAPER_1005) |
| 3 (Dirac) | Fermion / LENR | Kozima neutron-drop (PAPER_1061) |
| 4 (SCm) | Superconducting manifold | $V(\phi_0) = -\rho_{\text{SCm}}$ canonical |
| 5 (Mag) | Um magnetism | Heaviside amplifier (PAPER_1072) |
| 6 (Buoy) | F_U_Bi_i buoyancy | Variational EOM (PAPER_1065) |
| 7 (Aether) | Vacuum background | Two-component rho (PAPER_1051) |
| 8 (LENR) | Nuclear transmutation | COP parametric (PAPER_1081) |
| 9 (KK) | Kaluza-Klein 26D | $S_{26}^{(3)}$ compactification (PAPER_1080) |
"""

PKG_BLOCKS["PKG-S26"] = r"""<!-- PKG-S26-S225 -->

### Session 225 Phonon-Physics Upgrade: S26(3) Ramanujan Summation

> *Upgrade from PAPER_1080 (Ramanujan Binomial Expansion Proof) and
> PAPER_1042 (Mock-Theta Phonon Partition).  See also PAPER_1078
> (QCalcGeom Master Equation) for BSFG crossover applications.*

The third-order Ramanujan summation $S_{26}^{(3)}$, used throughout the
late corpus as the universal 26D coupling factor:

$$S_{26}^{(3)} = \sum_{n=0}^{\infty} \frac{(1/4)_n\,(1/2)_n\,(3/4)_n}{(n!)^3} \cdot \prod_{i=1}^{26}\left[1 + [\text{SSq}]\cdot e^{-\kappa\,i\,n/26}\right]$$

where $(a)_n = a(a+1)\cdots(a+n-1)$ is the Pochhammer symbol.

**Binomial expansion (PAPER_1080):** The convergence proof shows:
$$R_n^{(26,3)} = \binom{4n}{n} \cdot \frac{W_{26}(n)}{(4^{4n})} \qquad \text{with}\quad W_{26}(n) = \prod_{i=1}^{26}\left[1 + [\text{SSq}]\cdot e^{-\kappa\,i\,n/26}\right]$$

This sum converges absolutely for $|[\text{SSq}]| < 1$ (satisfied by $[\text{SSq}] = 0.57$)
and reduces to the classical Ramanujan $1/\pi$ series when $[\text{SSq}] \to 0$.

**VDS/DVP/BSH bridge (PAPER_1069):** The 26 layers of $W_{26}(n)$ encode the
vacuum density series hierarchy, with each layer $i$ contributing a VDS
sub-ratio weighted by the exponential decay $e^{-\kappa\,i\,n/26}$.

**Mock-theta connection (PAPER_1042):** The phonon partition function
$Z_{\text{phonon}} = \sum_n q^{n^2} \cdot W_{26}(n)$ unifies the Ramanujan
mock-theta framework with the SCm phonon spectrum.
"""

# Canonical ordering for PKG blocks
PKG_ORDER = ["PKG-GW", "PKG-AGN", "PKG-DM", "PKG-CLU", "PKG-YM", "PKG-LENR", "PKG-LAG", "PKG-S26"]


def extract_tags(content):
    """Extract YAML tags from frontmatter."""
    m = re.search(r'^tags:\s*\[([^\]]*)\]', content, re.MULTILINE)
    if m:
        return [t.strip().strip('"').strip("'") for t in m.group(1).split(',') if t.strip()]
    return []


def extract_title(content):
    """Extract title from YAML frontmatter."""
    m = re.search(r'^title:\s*["\']?([^"\'\n]*)["\']?', content, re.MULTILINE)
    return m.group(1).strip() if m else ""


def determine_pkgs(tags, title, content_head):
    """Determine which PKG blocks to insert based on tags, title, and content."""
    needed = set()

    # Tag-based matching
    for tag in tags:
        tag_clean = tag.strip()
        if tag_clean in TAG_TO_PKG:
            needed.update(TAG_TO_PKG[tag_clean])

    # Title keyword matching
    title_lower = title.lower()
    for keyword, pkgs in TITLE_KEYWORDS.items():
        if keyword.lower() in title_lower:
            needed.update(pkgs)

    # Content keyword matching (first 100 lines)
    head_lower = content_head.lower()
    content_keywords = {
        "gravitational wave": ["PKG-GW"], "strain amplitude": ["PKG-GW"],
        "neutron star merger": ["PKG-GW"], "binary inspiral": ["PKG-GW"],
        "eddington luminosity": ["PKG-AGN"], "blandford": ["PKG-AGN"],
        "jet power": ["PKG-AGN"], "m-sigma": ["PKG-AGN"],
        "nfw profile": ["PKG-DM"], "rotation curve": ["PKG-DM"],
        "dark matter halo": ["PKG-DM"],
        "intracluster medium": ["PKG-CLU"], "hydrostatic": ["PKG-CLU"],
        "cooling flow": ["PKG-CLU"], "sunyaev": ["PKG-CLU"],
        "yang-mills": ["PKG-YM"], "mass gap": ["PKG-YM"],
        "quark-gluon": ["PKG-YM"], "confinement": ["PKG-YM"],
        "lenr": ["PKG-LENR"], "transmutation": ["PKG-LENR"],
        "kozima": ["PKG-LENR"], "cold fusion": ["PKG-LENR"],
        "lagrangian": ["PKG-LAG"], "euler-lagrange": ["PKG-LAG"],
        "variational": ["PKG-LAG"],
        "ramanujan": ["PKG-S26"], "s_{26}": ["PKG-S26"],
        "s26": ["PKG-S26"], "26-dimensional": ["PKG-S26"],
        "phonon": ["PKG-S26"], "scm vacuum": ["PKG-S26"],
        "buoyancy": ["PKG-S26"], "f_u_bi": ["PKG-S26"],
        "uqff": ["PKG-S26"], "alma": ["PKG-S26"],
    }
    for keyword, pkgs in content_keywords.items():
        if keyword in head_lower:
            needed.update(pkgs)

    # Always include PKG-S26 if any other block is present (it's the universal connector)
    if needed and "PKG-S26" not in needed:
        needed.add("PKG-S26")

    return needed


def find_insertion_point(content):
    """Find where to insert PKG blocks — between body and first appendix.
    
    Returns the character index where PKG blocks should be inserted.
    Skips YAML frontmatter delimiters to avoid inserting at top of file.
    """
    lines = content.split('\n')
    
    # Skip YAML frontmatter: starts with '---', ends with next '---'
    start_search = 0
    if lines and lines[0].strip() == '---':
        for i in range(1, len(lines)):
            if lines[i].strip() == '---':
                start_search = i + 1
                break
    
    # Patterns that mark the start of appendix/footer sections
    appendix_patterns = [
        r'^## Appendix',
        r'^## §A\b',
        r'^## §SM\b',
        r'^## §B\b',
        r'^## Calibration Constants',
        r'^## Calibration',
        r'^## Session 225 Cross-References',
    ]
    
    for i in range(start_search, len(lines)):
        line = lines[i]
        for pat in appendix_patterns:
            if re.match(pat, line):
                # Insert before this line
                # Check if there's a --- separator above
                insert_line = i
                if i > 0 and lines[i-1].strip() == '---':
                    insert_line = i - 1
                if i > 1 and lines[i-1].strip() == '' and lines[i-2].strip() == '---':
                    insert_line = i - 2
                    
                # Calculate character position
                char_pos = sum(len(lines[j]) + 1 for j in range(insert_line))
                return char_pos
    
    # No appendix found — insert at end
    return len(content)


def process_paper(filepath):
    """Process a single paper, inserting PKG blocks into the body."""
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
    
    # Skip if already has any PKG marker
    if 'PKG-' in content:
        return 0, "already has PKG"
    
    # Extract metadata
    tags = extract_tags(content)
    title = extract_title(content)
    lines = content.split('\n')
    content_head = '\n'.join(lines[:min(100, len(lines))])
    
    # Determine which PKG blocks to insert
    needed_pkgs = determine_pkgs(tags, title, content_head)
    if not needed_pkgs:
        return 0, "no matching topics"
    
    # Order the PKG blocks canonically
    ordered_pkgs = [p for p in PKG_ORDER if p in needed_pkgs]
    
    # Build the combined upgrade text
    upgrade_text = "\n\n---\n\n## Session 225: Late-Corpus Physics Integration (PAPER_1000-1081)\n\n"
    upgrade_text += "> *The following physics upgrades incorporate equations, mechanisms, and\n"
    upgrade_text += "> derivations from the late-corpus papers (Sessions 219-225, PAPER_1000-1081).\n"
    upgrade_text += "> These represent body-level integrations of phonon physics, buoyancy\n"
    upgrade_text += "> formulations, and S₂₆⁽³⁾ Ramanujan corrections into this paper's domain.*\n\n"
    
    for pkg_id in ordered_pkgs:
        upgrade_text += PKG_BLOCKS[pkg_id]
    
    # Find insertion point
    insert_pos = find_insertion_point(content)
    
    # Insert
    new_content = content[:insert_pos] + upgrade_text + "\n" + content[insert_pos:]
    
    with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
        f.write(new_content)
    
    return len(ordered_pkgs), ", ".join(ordered_pkgs)


def main():
    print("=" * 72)
    print("integrate_body_physics.py — Session 225 Body-Level Physics Integration")
    print("=" * 72)
    
    # Gather all paper files
    paper_files = sorted(glob.glob(os.path.join(WP_DIR, "PAPER_*.md")))
    # Also include root-level papers
    paper_files += sorted(glob.glob(os.path.join(REPO, "PAPER_*.md")))
    
    print(f"Scanning {len(paper_files)} paper files...\n")
    
    updated = 0
    skipped_pkg = 0
    skipped_no_match = 0
    total_blocks = 0
    errors = 0
    
    for filepath in paper_files:
        filename = os.path.basename(filepath)
        try:
            count, detail = process_paper(filepath)
            if count > 0:
                print(f"  [+] {filename}: {count} PKG blocks ({detail})")
                updated += 1
                total_blocks += count
            elif "already" in detail:
                skipped_pkg += 1
            else:
                skipped_no_match += 1
        except Exception as e:
            errors += 1
            print(f"  [!] ERROR {filename}: {e}")
    
    print()
    print(f"  Updated:                {updated} papers")
    print(f"  Total PKG blocks added: {total_blocks}")
    print(f"  Skipped (has PKG):      {skipped_pkg}")
    print(f"  Skipped (no topics):    {skipped_no_match}")
    print(f"  Errors:                 {errors}")
    print("=" * 72)


if __name__ == "__main__":
    main()
