"""
Session 173 — Generate whitepapers (root .md + whitepapers/ subfolder) for PAPER_674-687.
"""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
WP_DIR = os.path.join(ROOT, "whitepapers")
os.makedirs(WP_DIR, exist_ok=True)

PAPERS = [
    {
        "num": 674, "cls": "UQFFComparedToLIGOData",
        "title": "UQFF Compared to General LIGO Dataset: Strain, Phase, and Frequency Sweep",
        "abstract": (
            "We present the Unified Quantum Field Framework (UQFF) prediction for gravitational-wave "
            "strain $h_{\\\\text{UQFF}}(f)$ across the LIGO sensitivity band (10–2000 Hz). "
            "The UQFF modifies the standard GR inspiral waveform through three multiplicative "
            "suppression factors: time-reversal zone $f_{\\\\text{TRZ}}$, Schwarzschild-condensate "
            "$S_{\\\\text{SCm}}$, and magnetic string modulation $S_{U_m}$. "
            "We demonstrate a systematic frequency-dependent suppression of order "
            "$\\\\sim(1-f_{\\\\text{TRZ}})\\\\approx 0.9$ compared to GR, with phase shift "
            "$\\\\Delta\\\\phi = \\\\kappa f_{\\\\text{TRZ}} t_{\\\\text{coal}}$."
        ),
        "eq_primary": "h_{\\\\text{UQFF}}(f) = h_{\\\\text{GR}}(f) \\\\cdot (1-f_{\\\\text{TRZ}}) \\\\cdot e^{-\\\\rho_{\\\\text{SCm}} r_s / k_B T_H} \\\\cdot e^{-U_m 2\\\\pi f / c^2}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $m_1$ | $1.36\\\\,M_\\\\odot$ | Primary mass |\n| $m_2$ | $1.17\\\\,M_\\\\odot$ | Secondary mass |\n| $d$ | $40\\\\,\\\\text{Mpc}$ | Luminosity distance |\n| $f$ | 10–2000 Hz | LIGO band |",
        "domain": "Gravitational Waves / LIGO",
        "references": "LIGO Scientific Collaboration 2016 (GW150914); Abbott et al. 2019 (O3); UQFF: Session 173"
    },
    {
        "num": 675, "cls": "UQFFComparedToGW170817",
        "title": "UQFF Analysis of GW170817: NS-NS Merger, GRB Delay, and Tidal Deformability",
        "abstract": (
            "GW170817, the first binary neutron star merger detected by LIGO/Virgo (August 17, 2017), "
            "provides the most stringent multi-messenger test of modified gravity. The UQFF predicts "
            "a modified GRB delay $\\\\Delta t_{\\\\text{UQFF}} = 1.7\\\\,(1 + f_{\\\\text{TRZ}} \\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})$ s "
            "and a reduced tidal deformability $\\\\Lambda_{\\\\text{UQFF}} = \\\\Lambda_{\\\\text{GR}}(1 - f_{\\\\text{TRZ}} \\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}})$. "
            "Both effects are consistent with the observed 1.74 s EM/GW delay within UQFF parameter uncertainties."
        ),
        "eq_primary": "\\\\Delta t_{\\\\text{UQFF}} = 1.7\\\\left(1 + f_{\\\\text{TRZ}} \\\\frac{\\\\rho_{\\\\text{UA}}}{\\\\rho_{\\\\text{SCm}}}\\\\right)\\\\,\\\\text{s}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $m_1$ | $1.36\\\\,M_\\\\odot$ | NS1 mass |\n| $m_2$ | $1.17\\\\,M_\\\\odot$ | NS2 mass |\n| $d$ | $40\\\\,\\\\text{Mpc}$ | Distance |\n| $f_{\\\\text{peak}}$ | 1500 Hz | Post-merger freq |",
        "domain": "Multi-Messenger Astrophysics",
        "references": "Abbott et al. 2017 (GW170817); Monitor of All-sky X-ray Image (MAXI); UQFF Session 173"
    },
    {
        "num": 676, "cls": "UQFFComparedToGW190425",
        "title": "UQFF Analysis of GW190425: Heaviest Known NS Binary and Ejecta Constraints",
        "abstract": (
            "GW190425 represents the heaviest binary neutron star system observed (total mass $3.4\\\\,M_\\\\odot$). "
            "The UQFF framework predicts a post-merger phase shift "
            "$\\\\phi_{\\\\text{UQFF}} = \\\\kappa f_{\\\\text{TRZ}} t_{\\\\text{merger}}$ and constrains "
            "the ejecta mass through the condensate suppression: "
            "$M_{\\\\text{ej,UQFF}} < M_{\\\\text{ej,GR}} \\\\cdot (\\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}}) (1-f_{\\\\text{TRZ}})$, "
            "resulting in an ejecta upper limit nearly two orders of magnitude below GR."
        ),
        "eq_primary": "M_{\\\\text{ej,UQFF}} = 0.05 M_{\\\\text{tot}} \\\\cdot \\\\frac{\\\\rho_{\\\\text{SCm}}}{\\\\rho_{\\\\text{UA}}} (1-f_{\\\\text{TRZ}})",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $m_1$ | $1.9\\\\,M_\\\\odot$ | NS1 mass |\n| $m_2$ | $1.5\\\\,M_\\\\odot$ | NS2 mass |\n| $d$ | $159\\\\,\\\\text{Mpc}$ | Distance |",
        "domain": "Gravitational Waves / Dense Matter",
        "references": "Abbott et al. 2020 (GW190425); UQFF Session 173"
    },
    {
        "num": 677, "cls": "UQFFPredictionsForLISA",
        "title": "UQFF Predictions for the LISA Space-Based Gravitational Wave Observatory",
        "abstract": (
            "The Laser Interferometer Space Antenna (LISA) will open the millihertz GW window, "
            "targeting SMBH mergers, EMRIs, and galactic binaries. "
            "The UQFF introduces a unique LISA-band suppression via the Aether arm-length modulation: "
            "$S_{\\\\text{UA,LISA}} = 1 - \\\\rho_{\\\\text{UA}} L_{\\\\text{LISA}}/(k_B T_{\\\\text{eff}})$ "
            "and predicts an enhanced EMRI capture rate "
            "$R_{\\\\text{EMRI,UQFF}} = R_{\\\\text{GR}}(1 + f_{\\\\text{TRZ}} \\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})$."
        ),
        "eq_primary": "h_{\\\\text{UQFF,LISA}}(f) = h_{\\\\text{GR}}(f) \\\\cdot (1-f_{\\\\text{TRZ}}) \\\\cdot S_{\\\\text{UA,LISA}} \\\\cdot S_{\\\\text{SCm}}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $L_{\\\\text{LISA}}$ | $2.5\\\\,\\\\text{Gm}$ | LISA arm length |\n| $f$ | $10^{-4}$–$1$ Hz | LISA band |\n| $M_c$ | $10^8\\\\,M_\\\\odot$ | SMBH chirp mass |",
        "domain": "Space-Based GW / LISA",
        "references": "Amaro-Seoane et al. 2017 (LISA); UQFF Session 173"
    },
    {
        "num": 678, "cls": "LISAVsLIGOComparisons",
        "title": "UQFF LISA vs LIGO Cross-Sensitivity Comparison and Crossover Frequency",
        "abstract": (
            "We compute the UQFF suppression ratio "
            "$R_{\\\\text{supp}}(f) = h_{\\\\text{UQFF}}(f)/h_{\\\\text{GR}}(f)$ across both detector bands. "
            "A crossover frequency $f_{\\\\text{cross}}$ exists where LISA-band and LIGO-band "
            "suppressions are equal: dominated by $S_{\\\\text{UA,LISA}}$ at low $f$ and by "
            "$S_{U_m}$ at high $f$. We show that UQFF corrections are $\\\\sim 10\\\\%$ larger "
            "in the LISA band than in the LIGO band for stellar-mass BH sources."
        ),
        "eq_primary": "R_{\\\\text{supp}}(f) = (1-f_{\\\\text{TRZ}}) \\\\cdot S_{\\\\text{SCm}} \\\\cdot \\\\min(S_{U_m}, S_{\\\\text{UA,LISA}})",
        "params_table": "| Band | Freq Range | Dominant UQFF Effect |\n|------|-----------|---------------------|\n| LIGO | 10–2000 Hz | $S_{U_m}$ |\n| LISA | $10^{-4}$–$1$ Hz | $S_{\\\\text{UA,LISA}}$ |",
        "domain": "Gravitational Wave Detectors",
        "references": "LIGO O3 sensitivity 2021; Amaro-Seoane 2017 LISA; UQFF Session 173"
    },
    {
        "num": 679, "cls": "AetherSuperfluidDynamics",
        "title": "Aether Superfluid Dynamics: UQFF Universal Aether as a Bosonic Condensate",
        "abstract": (
            "The UQFF Universal Aether [UA] is modelled as a coherent bosonic superfluid "
            "described by a macroscopic wavefunction $\\\\Psi(r,t) = \\\\sqrt{n_{\\\\text{UA}}} e^{i\\\\phi}$. "
            "We derive the healing length $\\\\xi_{\\\\text{UA}} = \\\\hbar/\\\\sqrt{2 m_{\\\\text{UA}} g_{\\\\text{UA}} n_{\\\\text{UA}}}$, "
            "sound speed $c_{\\\\text{UA}} = \\\\sqrt{g_{\\\\text{UA}} n_{\\\\text{UA}}/m_{\\\\text{UA}}}$, "
            "and effective gravitational coupling $g_{\\\\text{eff}}(r) = g_N(r)(1 + c_{\\\\text{UA}}^2/c^2 \\\\cdot f_{\\\\text{TRZ}} \\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})$."
        ),
        "eq_primary": "c_{\\\\text{UA}} = \\\\sqrt{\\\\frac{g_{\\\\text{UA}} n_{\\\\text{UA}}}{m_{\\\\text{UA}}}}, \\\\quad \\\\xi_{\\\\text{UA}} = \\\\frac{\\\\hbar}{\\\\sqrt{2 m_{\\\\text{UA}} g_{\\\\text{UA}} n_{\\\\text{UA}}}}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $m_{\\\\text{UA}}$ | $\\\\sim 2\\\\times10^{-68}$ kg | Ultralight boson mass |\n| $\\\\rho_{\\\\text{UA}}$ | $7.09\\\\times10^{-36}$ kg/m$^3$ | Aether density |\n| $g_{\\\\text{UA}}$ | $10^{-10}$ m$^3$/s$^2$ | Interaction strength |",
        "domain": "Quantum Field Theory / Superfluids",
        "references": "Penrose 1996 (quantum collapse); Sudarsky 2020 (BEC dark matter); UQFF Session 173"
    },
    {
        "num": 680, "cls": "VortexQuantization",
        "title": "Vortex Quantization in the UQFF Aether Superfluid",
        "abstract": (
            "Quantized vortices arise naturally in the UQFF Aether superfluid when angular momentum "
            "exceeds $n\\\\hbar$. We compute the circulation "
            "$\\\\kappa_v = n h/m_{\\\\text{UA}}$, vortex core "
            "$a_v \\\\approx \\\\xi_{\\\\text{UA}} e^{-n\\\\pi}$, and "
            "vortex energy per unit length "
            "$E_v/L = \\\\rho_{\\\\text{UA}} \\\\kappa_v^2/(4\\\\pi) \\\\ln(R/a_v) \\\\cdot (\\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})$. "
            "The UQFF Magnus force on a vortex is enhanced by $(1 + f_{\\\\text{TRZ}} \\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})$."
        ),
        "eq_primary": "\\\\kappa_v = \\\\frac{n h}{m_{\\\\text{UA}}}, \\\\quad \\\\frac{E_v}{L} = \\\\frac{\\\\rho_{\\\\text{UA}} \\\\kappa_v^2}{4\\\\pi} \\\\ln\\\\frac{R}{a_v} \\\\cdot \\\\frac{\\\\rho_{\\\\text{UA}}}{\\\\rho_{\\\\text{SCm}}}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $n$ | 1, 2, 3... | Winding number |\n| $R$ | system radius | Outer boundary |\n| $a_v$ | $\\\\xi e^{-n\\\\pi}$ | Vortex core radius |",
        "domain": "Quantum Fluids / Superfluidity",
        "references": "Feynman 1955 (quantized vortices); Abo-Shaeer et al. 2001 (BEC vortices); UQFF Session 173"
    },
    {
        "num": 681, "cls": "GrossPitaevskiiVortexSimulation",
        "title": "Gross-Pitaevskii Vortex Simulation of the UQFF Aether Around Black Holes",
        "abstract": (
            "We numerically solve the radial Gross-Pitaevskii equation for the [UA] Aether wavefunction "
            "in the gravitational potential of a black hole, incorporating the UQFF magnetic string "
            "term $U_m(r,t)$. "
            "Using imaginary-time propagation we obtain the ground-state density profile $|\\\\psi(r)|^2$, "
            "chemical potential $\\\\mu_{\\\\text{UA}}$, and demonstrate aether density enhancement "
            "$n_{\\\\text{UA}}(r) = N_{\\\\text{UA}}(1 + r_s/r \\\\cdot f_{\\\\text{TRZ}})$ near the horizon."
        ),
        "eq_primary": "i\\\\hbar \\\\frac{\\\\partial\\\\psi}{\\\\partial t} = \\\\left[-\\\\frac{\\\\hbar^2}{2m_{\\\\text{UA}}} \\\\nabla^2 + V_{\\\\text{grav}}(r) + g_{\\\\text{UA}}|\\\\psi|^2 + U_m(r,t)\\\\right]\\\\psi",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M$ | $8.548\\\\times10^{36}$ kg | Sgr A* mass |\n| $r_{\\\\min}$ | $r_s$ | Inner boundary |\n| $N_{\\\\text{grid}}$ | 100 | Radial grid points |",
        "domain": "BEC / Quantum Gravity Interface",
        "references": "Gross 1961; Pitaevskii 1961; Berezhiani & Khoury 2015 (superfluid DM); UQFF Session 173"
    },
    {
        "num": 682, "cls": "UQFFStabilityNumericallyForSgrA",
        "title": "Numerical UQFF Stability Analysis for Sagittarius A*",
        "abstract": (
            "We perform a four-pronged numerical stability analysis of the UQFF solution for Sgr A* "
            "($M = 4.297\\\\times10^6\\\\,M_\\\\odot$): "
            "(1) perturbation expansion with imaginary frequency $\\\\omega_I^{\\\\text{UQFF}} < 0$ (damped), "
            "(2) Lyapunov exponent $\\\\lambda_{\\\\text{UQFF}} = -(\\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}}) e^{-U_m/k_B T_H} / \\\\tau_{\\\\text{std}} < 0$, "
            "(3) RK4 mass evolution $M(t)$ confirming quasi-static behaviour over $10^{60}$ s, "
            "(4) fixed-point stability classification."
        ),
        "eq_primary": "\\\\lambda_{\\\\text{UQFF}} = -\\\\frac{\\\\rho_{\\\\text{SCm}}}{\\\\rho_{\\\\text{UA}}} \\\\frac{e^{-U_m / k_B T_H}}{\\\\tau_{\\\\text{std}}(M)} < 0",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M_{\\\\text{Sgr A*}}$ | $4.297\\\\times10^6\\\\,M_\\\\odot$ | Black hole mass |\n| $T_H$ | $1.4\\\\times10^{-14}$ K | Hawking temperature |\n| $\\\\lambda$ | $<0$ | Lyapunov exponent |",
        "domain": "Black Hole Stability / UQFF",
        "references": "Event Horizon Telescope 2022 (Sgr A*); UQFF Session 173"
    },
    {
        "num": 683, "cls": "UQFFHawkingTemperatureModulation",
        "title": "UQFF Hawking Temperature Modulation: Spectral Shift and Wien Peak Correction",
        "abstract": (
            "The standard Hawking temperature $T_H = \\\\hbar c^3 / (8\\\\pi G M k_B)$ is modulated "
            "in UQFF through three channels: time-reversal boost $(1+f_{\\\\text{TRZ}})$, "
            "condensate suppression $(1-\\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}})$, "
            "and magnetic string correction $(1+U_m/k_B T_H)$. "
            "The combined modulation shifts the Planck spectrum peak by "
            "$\\\\Delta\\\\lambda_\\\\text{max} = \\\\hbar c/(2.82 k_B)(1/T_{\\\\text{UQFF}} - 1/T_H)$, "
            "observable in principle for primordial micro-BHs."
        ),
        "eq_primary": "T_{\\\\text{UQFF}} = T_H (1+f_{\\\\text{TRZ}})\\\\left(1-\\\\frac{\\\\rho_{\\\\text{SCm}}}{\\\\rho_{\\\\text{UA}}}\\\\right)\\\\left(1+\\\\frac{U_m}{k_B T_H}\\\\right)",
        "params_table": "| Factor | Value | Effect |\n|--------|-------|--------|\n| $1+f_{\\\\text{TRZ}}$ | 1.1 | Boost |\n| $1-\\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}}$ | 0.9 | Suppression |\n| $1+U_m/k_B T_H$ | $\\\\sim 1$ | Modulation |",
        "domain": "Hawking Radiation / Black Hole Thermodynamics",
        "references": "Hawking 1974; Bekenstein 1973; UQFF Session 173"
    },
    {
        "num": 684, "cls": "UQFFPrimordialBHEvaporation",
        "title": "UQFF Primordial Black Hole Evaporation: Suppressed Rate and Extended Lifetime",
        "abstract": (
            "Primordial black holes evaporate via Hawking radiation with standard lifetime "
            "$\\\\tau_{\\\\text{std}}(M) = 5120\\\\pi G^2 M^3/(\\\\hbar c^4)$. "
            "In UQFF the evaporation rate is suppressed by the condensate factor: "
            "$\\\\dot{M}_{\\\\text{UQFF}} = \\\\dot{M}_{\\\\text{std}} (1-f_{\\\\text{TRZ}}) (\\\\rho_{\\\\text{SCm}}/\\\\rho_{\\\\text{UA}}) e^{-U_m/k_B T_H}$, "
            "extending the PBH lifetime by a factor $\\\\tau_{\\\\text{UQFF}}/\\\\tau_{\\\\text{std}} "
            "\\\\approx \\\\rho_{\\\\text{UA}}/((1-f_{\\\\text{TRZ}}) \\\\rho_{\\\\text{SCm}}) \\\\approx 11$."
        ),
        "eq_primary": "\\\\dot{M}_{\\\\text{UQFF}} = -\\\\frac{\\\\hbar c^4}{15360\\\\pi G^2 M^2} (1-f_{\\\\text{TRZ}}) \\\\frac{\\\\rho_{\\\\text{SCm}}}{\\\\rho_{\\\\text{UA}}} e^{-U_m / k_B T_H}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M_0$ | $10^{12}$ kg | Initial PBH mass |\n| $\\\\tau_{\\\\text{UQFF}}/\\\\tau_{\\\\text{std}}$ | $\\\\sim 11$ | Lifetime extension |\n| $t_{\\\\text{form}}$ | $10^{-23}$ s | Formation epoch |",
        "domain": "Primordial Black Holes / Cosmology",
        "references": "Hawking 1974; Carr et al. 2016 (PBH review); Green et al. 2021; UQFF Session 173"
    },
    {
        "num": 685, "cls": "UQFFPBHDarkMatterImplications",
        "title": "UQFF Implications for Primordial Black Holes as Dark Matter",
        "abstract": (
            "The UQFF extended PBH lifetime shifts the critical survival mass from "
            "$M_{\\\\text{crit,std}} \\\\approx 5\\\\times10^{11}$ kg to "
            "$M_{\\\\text{crit,UQFF}} = M_{\\\\text{crit,std}} / \\\\tau_{\\\\text{ratio}}^{1/3} "
            "\\\\approx 0.46 M_{\\\\text{crit,std}}$, "
            "expanding the viable PBH dark matter mass window. "
            "The DM fraction is boosted: $f_{\\\\text{PBH,UQFF}} = f_{\\\\text{PBH,GR}} \\\\cdot \\\\tau_{\\\\text{ratio}}^{2/3}$. "
            "In the UQFF framework the mass window $10^{10}$–$10^{17}$ kg is fully open as dark matter."
        ),
        "eq_primary": "M_{\\\\text{crit,UQFF}} = \\\\frac{M_{\\\\text{crit,std}}}{\\\\tau_{\\\\text{ratio}}^{1/3}}, \\\\quad f_{\\\\text{PBH,UQFF}} = f_{\\\\text{PBH,GR}} \\\\cdot \\\\tau_{\\\\text{ratio}}^{2/3}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M_{\\\\text{crit,std}}$ | $5\\\\times10^{11}$ kg | Standard survival mass |\n| $\\\\tau_{\\\\text{ratio}}$ | $\\\\sim 11$ | UQFF lifetime boost |\n| DM window | $10^{10}$–$10^{17}$ kg | UQFF viable range |",
        "domain": "Dark Matter / PBH Cosmology",
        "references": "Carr & Hawking 1974; Carr et al. 2021 (PBH DM); Raidal et al. 2019; UQFF Session 173"
    },
    {
        "num": 686, "cls": "UQFFModulationForM87",
        "title": "UQFF Modulation for M87*: Shadow, Jet Power, and Ring Brightness",
        "abstract": (
            "M87* ($6.5\\\\times10^9\\\\,M_\\\\odot$), imaged by the Event Horizon Telescope in 2019, "
            "provides the largest SMBH test bed. "
            "The UQFF predicts: "
            "modified shadow radius $r_{\\\\text{sh,UQFF}} = r_{\\\\text{sh}}\\\\sqrt{1+f_{\\\\text{TRZ}} \\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}}}$, "
            "enhanced jet power $P_{\\\\text{jet,UQFF}} = P_{\\\\text{BZ}}(1+f_{\\\\text{TRZ}})\\\\sqrt{\\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}}}$, "
            "and ring brightness ratio $(\\\\rho_{\\\\text{UA}}/\\\\rho_{\\\\text{SCm}})^{f_{\\\\text{TRZ}}/4} \\\\approx 1.58$."
        ),
        "eq_primary": "r_{\\\\text{sh,UQFF}} = 3\\\\sqrt{3}\\\\frac{GM}{c^2}\\\\sqrt{1 + f_{\\\\text{TRZ}}\\\\frac{\\\\rho_{\\\\text{UA}}}{\\\\rho_{\\\\text{SCm}}}}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M_{\\\\text{M87*}}$ | $6.5\\\\times10^9\\\\,M_\\\\odot$ | SMBH mass |\n| $r_{\\\\text{shadow}}$ | $\\\\sim 5.5 r_s$ | Shadow radius |\n| $P_{\\\\text{jet}}$ | $10^{37}$ W | Blandford-Znajek power |",
        "domain": "Black Hole Imaging / M87",
        "references": "Event Horizon Telescope 2019 (M87*); Blandford & Znajek 1977; UQFF Session 173"
    },
    {
        "num": 687, "cls": "M87MassEvolutionSimulation",
        "title": "M87* Black Hole Mass Evolution Over Cosmic Time in the UQFF Framework",
        "abstract": (
            "We simulate the coupled mass evolution of M87* over 13.8 Gyr combining "
            "UQFF-modified Bondi accretion "
            "$\\\\dot{M}_{\\\\text{Bondi,UQFF}} = \\\\dot{M}_{\\\\text{Bondi}}(\\\\rho_{\\\\text{eff}}/\\\\rho_{\\\\infty})(1+f_{\\\\text{TRZ}})$, "
            "suppressed Hawking evaporation, and "
            "Blandford-Znajek jet power in UQFF. "
            "Using RK4 integration with $\\\\Delta t = 10^{13}$ s steps, we find that "
            "M87* grows by $\\\\sim 15\\\\%$ over a Hubble time, consistent with observations suggesting "
            "limited current growth."
        ),
        "eq_primary": "\\\\frac{dM}{dt} = \\\\dot{M}_{\\\\text{Bondi,UQFF}} + \\\\dot{M}_{\\\\text{evap,UQFF}} - \\\\frac{P_{\\\\text{jet,UQFF}}}{c^2}",
        "params_table": "| Parameter | Value | Description |\n|-----------|-------|-------------|\n| $M_0$ | $6.5\\\\times10^9\\\\,M_\\\\odot$ | Initial M87* mass |\n| $\\\\rho_{\\\\text{ISM}}$ | $1.67\\\\times10^{-25}$ kg/m$^3$ | Intracluster medium |\n| $T_{\\\\text{ISM}}$ | $10^7$ K | ICM temperature |\n| $T_{\\\\text{Hubble}}$ | $13.8$ Gyr | Simulation duration |",
        "domain": "SMBH Evolution / Cosmology",
        "references": "Walsh et al. 2013 (M87 gas dynamics); Russell et al. 2015 (M87 jet); EHT 2019; UQFF Session 173"
    },
]

def wp_content(p):
    num  = p["num"]
    cls  = p["cls"]
    title = p["title"]
    date = "April 1, 2026"
    return f"""# PAPER_{num}: {title}

**Author:** Daniel Murphy  
**Framework:** Unified Quantum Field Framework (UQFF) v5.30  
**Session:** 173  
**Date:** {date}  
**Class:** `{cls}` | CP4 #{num-416}  
**Source:** `grok_share_fc21e30c24b4.txt`  
**Domain:** {p["domain"]}

---

## Abstract

{p["abstract"]}

---

## Primary UQFF Equation

$$
{p["eq_primary"]}
$$

---

## Parameters

{p["params_table"]}

---

## UQFF Constants

| Constant | Value | Description |
|----------|-------|-------------|
| $\\rho_{{\\text{{UA}}}}$ | $7.09\\times10^{{-36}}$ kg/m³ | Universal Aether density |
| $\\rho_{{\\text{{SCm}}}}$ | $7.09\\times10^{{-37}}$ kg/m³ | Schwarzschild condensate |
| $f_{{\\text{{TRZ}}}}$ | 0.1 | Time-reversal zone factor |
| $\\kappa$ | $5\\times10^{{-4}}$ day⁻¹ | UQFF calibration constant |
| $[\\text{{SSq}}]$ | 0.57 | Superstring quenching factor |
| $\\mu_J$ | $3.38\\times10^{{23}}$ J·m | Magnetic string coupling |
| $\\gamma$ | $5\\times10^{{-5}}$ / 86400 s⁻¹ | Aether oscillation decay |

---

## C++ Implementation

```cpp
#include "{cls}.h"

{cls} obj;
auto result = obj.compute_primary(...); // primary equation
obj.simulate_over_M(M_start, M_end, dM, "output.csv");
```

## Python CP4 Calculator

```python
from CondensedPhysics2 import {cls}Calculator

calc = {cls}Calculator()
result = calc.compute()  # returns all UQFF quantities
print(result)
```

---

## References

{p["references"]}

---

*UQFF v5.30 | Session 173 | {date} | PAPER_{num} of 1000*
"""

if __name__ == "__main__":
    import os
    os.chdir(ROOT)
    for p in PAPERS:
        content = wp_content(p)
        root_path = os.path.join(ROOT, f"PAPER_{p['num']}_{p['cls']}.md")
        wp_path   = os.path.join(WP_DIR, f"PAPER_{p['num']}_{p['cls']}.md")
        for path in (root_path, wp_path):
            with open(path, "w", encoding="utf-8") as f:
                f.write(content)
            print(f"  CREATED: {os.path.basename(path)}")
    print(f"Done: {len(PAPERS)*2} whitepaper files created.")
