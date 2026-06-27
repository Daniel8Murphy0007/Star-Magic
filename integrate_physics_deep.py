#!/usr/bin/env python3
"""
integrate_physics_deep.py — Session 225 Deep Physics Integration
════════════════════════════════════════════════════════════════════════════════

Inserts actual physics equations and derivations from PAPER_1000-1081 into the
BODY of earlier papers (001-999) where topically relevant.

This is NOT an appendix-only update.  Each upgrade package inserts a clearly
marked "### Session 225 Phonon-Physics Upgrade" subsection into the paper body,
placed BEFORE the first appendix/calibration/SM-Anchor heading.

Upgrade Packages
───────────────────
  PKG-GW   : GW strain phonon modulation        (from PAPER_1000/1022)
  PKG-AGN  : Buoyancy-corrected Eddington        (from PAPER_1002/1037)
  PKG-DM   : Phonon-modified NFW density          (from PAPER_1015/1019)
  PKG-CLU  : ICM buoyancy force profile           (from PAPER_1039/1041)
  PKG-YM   : Yang-Mills BCS phonon mass gap       (from PAPER_1318/1070)
  PKG-LENR : VDS LENR transmutation dynamics      (from PAPER_1060/1061/1081)
  PKG-LAG  : 9-Sector UQFF Lagrangian foundation  (from PAPER_1066)
  PKG-S26  : S₂₆⁽³⁾ Ramanujan definition          (from PAPER_1042/1080)

Idempotent: each package has a unique marker; will not re-insert.
"""

import os
import re
import glob
from typing import Optional, List, Tuple, Set

REPO = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(REPO, "whitepapers")

# ═══════════════════════════════════════════════════════════════════════════════
# UPGRADE PACKAGES
# ═══════════════════════════════════════════════════════════════════════════════

PACKAGES = {}

# ── PKG-GW: Gravitational-Wave Phonon Strain Modulation ──────────────────────
PACKAGES["PKG-GW"] = {
    "marker": "<!-- PKG-GW-S225 -->",
    "match_tags": {"GW", "gravitational-wave", "gravitational wave", "strain",
                   "LIGO", "BNS", "GW170817", "GW190425", "neutron-star merger",
                   "NS merger", "inspiral", "merger", "ringdown"},
    "match_title_kw": [r"GW\d{6}", r"gravitational.wave", r"strain",
                       r"neutron.star.merger", r"NS.merger", r"inspiral",
                       r"coalescence", r"LIGO", r"Virgo", r"merger.*phonon",
                       r"BNS", r"binary.*neutron"],
    "exclude_papers": set(range(1000, 1082)),  # don't upgrade late-corpus
    "content": r"""
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

""",
}

# ── PKG-AGN: Buoyancy-Corrected Eddington Luminosity ─────────────────────────
PACKAGES["PKG-AGN"] = {
    "marker": "<!-- PKG-AGN-S225 -->",
    "match_tags": {"AGN", "quasar", "Seyfert", "accretion", "Eddington",
                   "jet", "blazar", "active galactic", "3C273", "M87",
                   "SgrA", "Sgr A", "Centaurus A", "TON618", "SMBH",
                   "active galaxy", "radio galaxy"},
    "match_title_kw": [r"AGN", r"quasar", r"Seyfert", r"accretion",
                       r"Eddington", r"blazar", r"jet", r"active.galact",
                       r"radio.galax", r"3C273", r"M87", r"TON618",
                       r"Centaurus.A", r"SMBH.*accret", r"Blandford",
                       r"Sgr.*A"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
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

**Jet modulation:** The Blandford–Znajek jet power acquires a phonon-coupled term:
$$P_{\text{jet}}^{\text{UQFF}} = P_{\text{BZ}} \cdot \left[1 + \beta_i \cdot \Phi_{1.25\,\text{THz}} \cdot \left(\frac{B}{B_{\text{crit}}}\right)^2\right]$$

where $\Phi_{1.25\,\text{THz}} = \cos(\omega_{\text{SCm}} \cdot t)$ modulates jet power at the phonon frequency.

**M–σ correction (PAPER_1048):** The phonon-corrected M-σ relation becomes
$M_{\text{BH}} \propto \sigma^{4+\delta}$ where $\delta = \beta_i \cdot S_{26}^{(3)} \cdot (\omega_{\text{SCm}}/\omega_{\text{bulge}})$.

""",
}

# ── PKG-DM: Phonon-Modified NFW Dark Matter Profile ──────────────────────────
PACKAGES["PKG-DM"] = {
    "marker": "<!-- PKG-DM-S225 -->",
    "match_tags": {"dark matter", "dark-matter", "NFW", "halo",
                   "rotation curve", "DM", "SIDM", "virial",
                   "CDM", "WIMP", "dark-matter-halo"},
    "match_title_kw": [r"dark.matter", r"NFW", r"halo", r"rotation.curve",
                       r"SIDM", r"virial.*dark", r"CDM", r"WIMP",
                       r"DM.halo", r"dark.halo"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
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

""",
}

# ── PKG-CLU: Galaxy Cluster ICM Buoyancy Force ───────────────────────────────
PACKAGES["PKG-CLU"] = {
    "marker": "<!-- PKG-CLU-S225 -->",
    "match_tags": {"cluster", "ICM", "galaxy cluster", "galaxy-cluster",
                   "X-ray cluster", "cool core", "cool-core",
                   "Perseus", "Coma", "Virgo", "Abell",
                   "intracluster", "AGN feedback", "SZ",
                   "Sunyaev-Zeldovich", "radio relic", "cluster merger"},
    "match_title_kw": [r"cluster", r"ICM", r"intracluster", r"cool.core",
                       r"Perseus", r"Coma", r"Virgo", r"Abell",
                       r"X.ray.*cluster", r"AGN.*feedback", r"SZ.*effect",
                       r"Sunyaev", r"radio.*relic", r"hydrostatic.*bias"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
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

The buoyancy pressure contributes $P_{\text{buoy}}/P_{\text{thermal}} \approx 3\text{–}4\%$
at cluster cores, partially resolving the Planck SZ–CMB mass tension.

**Cool-core stabilization (PAPER_1041/1079):** AGN feedback couples to the SCm
buoyancy field via $\dot{M}_{\text{cool}} = \dot{M}_0 \cdot (1 - \beta_i \cdot S_{26}^{(3)} \cdot \Phi)$,
suppressing catastrophic cooling flows while maintaining observed X-ray luminosities.

**Phonon frequency coupling:** $\omega_{\text{SCm}} = 2\pi \times 1.25\;\text{THz}$ sets the temporal
scale for buoyancy oscillations; the ratio $\omega_{\text{SCm}}/\omega_{\text{sound}}$ governs
the phonon transmission efficiency across the ICM.

""",
}

# ── PKG-YM: Yang-Mills BCS Phonon Mass Gap ──────────────────────────────────
PACKAGES["PKG-YM"] = {
    "marker": "<!-- PKG-YM-S225 -->",
    "match_tags": {"Yang-Mills", "Yang Mills", "mass gap", "mass-gap",
                   "QCD", "confinement", "gluon", "QGP",
                   "quark-gluon", "quark gluon", "deconfinement",
                   "running coupling", "strong force", "color"},
    "match_title_kw": [r"Yang.Mills", r"mass.gap", r"QCD", r"confine",
                       r"gluon", r"QGP", r"quark.gluon", r"deconfin",
                       r"running.coupling", r"strong.force", r"color.force",
                       r"ALICE", r"Pb.Pb", r"centrality"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
### Session 225 Phonon-Physics Upgrade: Yang-Mills BCS Phonon Mass Gap

> *Upgrade from PAPER_1318 (Yang-Mills Mass Gap via SCm BCS Phonon) and
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
\approx 1.736\;\text{GeV}$ at the 9-sector Lagrangian closure (PAPER_1066, §2).

**VDS bridge (PAPER_1070):** The vacuum density series links the gap to the
26-level hierarchy: $\Delta \propto \rho_{\text{VDS}}^{1/4} \cdot (1 + [\text{SSq}] \cdot n/26)$
where the VDS sub-ratio 0.108 places confinement in the sub-threshold regime.

**QGP transition (PAPER_1004/1007):** At $T > T_c \approx 170\;\text{MeV}$, the phonon
coupling weakens ($\alpha_s \to 0$) and the gap closes, reproducing the
deconfinement phase transition observed at ALICE/LHC.

""",
}

# ── PKG-LENR: VDS LENR Transmutation Dynamics ────────────────────────────────
PACKAGES["PKG-LENR"] = {
    "marker": "<!-- PKG-LENR-S225 -->",
    "match_tags": {"LENR", "cold fusion", "Kozima", "Widom-Larsen",
                   "Widom Larsen", "transmutation", "excess heat",
                   "Fleischmann", "palladium", "Pd-D", "neutron-drop",
                   "heavy electron", "nuclear reaction", "lattice confinement"},
    "match_title_kw": [r"LENR", r"cold.fusion", r"Kozima", r"Widom.Larsen",
                       r"transmut", r"excess.heat", r"Fleischmann",
                       r"palladium", r"Pd.D", r"neutron.drop",
                       r"heavy.electron", r"lattice.confine"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
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

The VDS factor $(1 + [\text{SSq}] \cdot n/26)$ provides ~470× amplification via
the 26-level vacuum density ladder at resonance ($\omega = \omega_{\text{SCm}}$).

**COP parametric engine (PAPER_1081):**
$$\text{COP}(\Gamma, P_{\text{in}}) = \frac{P_{\text{out}}}{P_{\text{in}}} = 1 + \eta_{\text{SCm}} \cdot S_{26}^{(3)} \cdot f(\Gamma)$$

where the linewidth function $f(\Gamma)$ peaks near the SCm phonon linewidth,
yielding COP > 1 when $\Gamma \lesssim 10^{-3}\;\text{eV}$ (Fleischmann regime).

**Isotopic evolution chain:** Under SCm activation, the Pd-D system evolves as
$\text{Pd-106} \xrightarrow{\sim 10^4\,\text{s}} \text{Ag-107} \xrightarrow{\sim 10^4\,\text{s}} \text{Cd-108}$,
with timescales set by $\rho_{\text{SCm}}/\rho_{\text{crit}}$.

""",
}

# ── PKG-LAG: UQFF 9-Sector Lagrangian Foundation ────────────────────────────
PACKAGES["PKG-LAG"] = {
    "marker": "<!-- PKG-LAG-S225 -->",
    "match_tags": {"Lagrangian", "field theory", "variational",
                   "Euler-Lagrange", "action principle", "EOM",
                   "first principles", "unified field", "unification",
                   "9-sector", "nine-sector"},
    "match_title_kw": [r"Lagrangian", r"field.theory", r"variation",
                       r"Euler.Lagrange", r"action.princ", r"EOM",
                       r"first.princ", r"unified.field", r"unifi",
                       r"9.sector", r"nine.sector", r"F_U.*master"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
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
| 2 (YM) | Yang-Mills gauge | $m_{\text{gap}} = 1.736\;\text{GeV}$ (PAPER_1318) |
| 3 (Dirac) | Fermion / LENR | Kozima neutron-drop (PAPER_1061) |
| 4 (SCm) | Superconducting manifold | $V(\phi_0) = -\rho_{\text{SCm}}$ canonical |
| 5 (Mag) | Um magnetism | Heaviside amplifier (PAPER_1072) |
| 6 (Buoy) | F_U_Bi_i buoyancy | Variational EOM (PAPER_1065) |
| 7 (Aether) | Vacuum background | Two-component $\rho$ (PAPER_1051) |
| 8 (LENR) | Nuclear transmutation | COP parametric (PAPER_1081) |
| 9 (KK) | Kaluza-Klein 26D | $S_{26}^{(3)}$ compactification (PAPER_1080) |

""",
}

# ── PKG-S26: S₂₆⁽³⁾ Ramanujan Definition ────────────────────────────────────
PACKAGES["PKG-S26"] = {
    "marker": "<!-- PKG-S26-S225 -->",
    "match_tags": {"26D", "Ramanujan", "S26", "compactification",
                   "Kaluza-Klein", "extra dimensions", "string theory",
                   "26 dimensions", "polylogarithm", "mock-theta",
                   "q-series"},
    "match_title_kw": [r"26.?D", r"Ramanujan", r"S.?26", r"compactif",
                       r"Kaluza.Klein", r"extra.dimen", r"string.theory",
                       r"26.*dimen", r"polylog", r"mock.theta",
                       r"q.series", r"Bosonic.*string"],
    "exclude_papers": set(range(1000, 1082)),
    "content": r"""
### Session 225 Phonon-Physics Upgrade: S₂₆⁽³⁾ Ramanujan Summation

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

""",
}


# ═══════════════════════════════════════════════════════════════════════════════
# PAPER PROCESSING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

def extract_paper_number(filename: str) -> Optional[int]:
    m = re.search(r'PAPER_(\d+)', filename)
    return int(m.group(1)) if m else None


def extract_tags(content: str) -> Set[str]:
    m = re.search(r'^tags:\s*\[([^\]]*)\]', content, re.MULTILINE)
    if m:
        return {t.strip().strip("'\"").lower()
                for t in m.group(1).split(',') if t.strip()}
    return set()


def extract_title(content: str) -> str:
    m = re.search(r'^title:\s*["\']?([^"\'\n]*)["\']?', content, re.MULTILINE)
    return m.group(1).strip() if m else ""


def find_insertion_point(content: str) -> int:
    """Find the byte offset to insert upgrade sections.

    Strategy: Insert BEFORE the first of these headings:
      - ## Calibration
      - ## Appendix
      - ## §SM
      - ## §A
      - ## §B
      - ## Session 225 Cross-References
      - ## ---   (horizontal rule before appendix)

    If none found, insert before the last --- separator, or at end.
    """
    # Try to find anchor headings (in order of preference)
    patterns = [
        r'^##\s+Calibration',
        r'^##\s+(?:Appendix|§SM|§A\b|§B\b)',
        r'^##\s+Session\s+225\s+Cross-References',
        r'^##\s+Appendix:\s+Session\s+225',
        r'^##\s+Appendix:\s+Session\s+204',
        r'^---\s*$',  # Last resort: any --- separator in final quarter
    ]

    lines = content.split('\n')
    total = len(lines)

    for pattern in patterns[:-1]:
        for i, line in enumerate(lines):
            if re.match(pattern, line, re.IGNORECASE):
                # Found it; insert before this line
                # Back up over any preceding blank lines
                insert_line = i
                while insert_line > 0 and lines[insert_line - 1].strip() == '':
                    insert_line -= 1
                return sum(len(l) + 1 for l in lines[:insert_line])

    # Last resort: find --- separator in last 25% of file
    cutoff = int(total * 0.75)
    for i in range(total - 1, cutoff, -1):
        if re.match(r'^---\s*$', lines[i]):
            insert_line = i
            while insert_line > 0 and lines[insert_line - 1].strip() == '':
                insert_line -= 1
            return sum(len(l) + 1 for l in lines[:insert_line])

    # Absolute last resort: end of file
    return len(content)


def matches_package(pkg: dict, tags: Set[str], title: str,
                    content_head: str) -> bool:
    """Check if a paper matches a package's criteria."""
    # Tag matching (case-insensitive)
    pkg_tags_lower = {t.lower() for t in pkg["match_tags"]}
    if tags & pkg_tags_lower:
        return True

    # Title / content keyword matching
    search_text = (title + "\n" + content_head).lower()
    for kw_pattern in pkg["match_title_kw"]:
        if re.search(kw_pattern, search_text, re.IGNORECASE):
            return True

    return False


def process_paper(filepath: str) -> Tuple[int, List[str]]:
    """Process one paper, inserting all matching upgrade packages.
    Returns (number_of_packages_inserted, list_of_package_names).
    """
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()

    filename = os.path.basename(filepath)
    paper_num = extract_paper_number(filename)
    if paper_num is None:
        return 0, []

    tags = extract_tags(content)
    title = extract_title(content)
    # first ~80 lines for keyword matching
    content_head = '\n'.join(content.split('\n')[:80])

    inserted = []
    for pkg_name, pkg in PACKAGES.items():
        # Skip if paper is in exclude range
        if paper_num in pkg["exclude_papers"]:
            continue
        # Skip if already has this package's marker
        if pkg["marker"] in content:
            continue
        # Check if paper matches this package
        if not matches_package(pkg, tags, title, content_head):
            continue

        # Insert the upgrade content
        insert_pos = find_insertion_point(content)
        upgrade_text = pkg["marker"] + "\n" + pkg["content"]
        content = content[:insert_pos] + "\n" + upgrade_text + "\n" + content[insert_pos:]
        inserted.append(pkg_name)

    if inserted:
        with open(filepath, 'w', encoding='utf-8', newline='\n') as f:
            f.write(content)

    return len(inserted), inserted


def main():
    print("=" * 72)
    print("integrate_physics_deep.py — Session 225 Deep Physics Integration")
    print("=" * 72)

    # Gather all paper files
    paper_files = []
    for f in glob.glob(os.path.join(WP_DIR, "PAPER_*.md")):
        paper_files.append(f)
    for f in glob.glob(os.path.join(REPO, "PAPER_*.md")):
        paper_files.append(f)

    paper_files.sort(key=lambda x: extract_paper_number(os.path.basename(x)) or 0)
    print(f"Found {len(paper_files)} paper files\n")

    updated = 0
    skipped = 0
    total_pkgs = 0
    pkg_counts = {name: 0 for name in PACKAGES}

    for filepath in paper_files:
        count, names = process_paper(filepath)
        if count > 0:
            updated += 1
            total_pkgs += count
            for name in names:
                pkg_counts[name] += 1
            fname = os.path.basename(filepath)
            print(f"  [+] {fname}: {', '.join(names)}")
        else:
            skipped += 1

    print(f"\n{'─' * 72}")
    print(f"  Papers updated:  {updated}")
    print(f"  Papers skipped:  {skipped}")
    print(f"  Total packages:  {total_pkgs}")
    print(f"\n  Package breakdown:")
    for name, count in sorted(pkg_counts.items()):
        print(f"    {name:12s} : {count:4d} papers")
    print(f"{'─' * 72}")


if __name__ == "__main__":
    main()
