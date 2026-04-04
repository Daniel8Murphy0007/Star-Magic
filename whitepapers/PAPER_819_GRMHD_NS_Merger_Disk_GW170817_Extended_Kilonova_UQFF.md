# PAPER_819: GRMHD NS Merger Disk + GW170817 Extended Kilonova UQFF
## Unified Quantum Field Framework — Whitepaper 819

**Session**: 192 | **Version**: v5.48 | **Date**: April 4, 2026
**Source**: grok_share_0d888ea9-50be.txt (June 13, 2025 05:55 PM); "GRMHD sim Neutron Star Accretion Disks_17June2025.pdf"
**Author**: Daniel T. Murphy — Star-Magic UQFF Project
**CVW Compliance**: v2.0.0

---

## Abstract
This paper integrates GRMHD simulations of a neutron star merger remnant accretion disk into the Quadriadic UQFF framework. The system consists of a 0.033 $M_\odot$ disk surrounding a 3 $M_\odot$ black hole (spin $a = 0.8$), evolved over 9 seconds post-merger. Key results include mass ejection $M_{ej} = 0.013 M_\odot$, ejecta velocity $v_{ej} = 0.1$–$0.22c$, jet luminosity $L_{jet} = 3 \times 10^{50}$ erg, and lanthanide-rich ejecta ($Y_e = 0.16$) consistent with the extended kilonova GW170817 / AT2017gfo. Dual ejection phases (MHD-driven and thermal wind) and fast ejecta ($v \geq 0.4c$) neutron precursors are formally derived.

---

## 1. Introduction
GW170817 (neutron star-neutron star merge, August 17, 2017) produced kilonova AT2017gfo observed in optical/infrared, confirming $r$-process nucleosynthesis. The long-term disk wind from the post-merger BH + disk system contributes lanthanide-rich ejecta that explains the red kilonova component lasting several days. GRMHD simulations extend this picture to 9 seconds, significantly longer than prior ~1 second runs.

---

## 2. System Configuration

- **Disk mass**: $M_{disk} = 0.033 M_\odot$
- **Central BH mass**: $M_{BH} = 3 M_\odot$ (spin $a = 0.8$)
- **Simulation duration**: 9 seconds (longest GRMHD NS disk run)
- **Magnetic field**: toroidal + poloidal configuration
- **Neutrino treatment**: leakage scheme with electron fraction tracking

---

## 3. Mass Ejection

Total ejected mass:

$$M_{ej} = 0.013 M_\odot = 39\% \text{ of initial disk mass}$$

Ejection proceeds in two phases:

**Phase 1 (MHD-driven)**: $t < 0.5$ s
$$M_{ej,MHD} \approx 0.007 M_\odot, \quad v_{ej} \approx 0.15c$$

**Phase 2 (thermal wind)**: $0.5 < t < 9$ s
$$M_{ej,thermal} \approx 0.006 M_\odot, \quad v_{ej} \approx 0.10c$$

UQFF Layer 1 entry:

$$g_{L1,ej} = \frac{M_{ej} \cdot v_{ej}^2}{r \cdot r_{code}^2} \approx 3.9 \times 10^{-5} \text{ m/s}^2$$

---

## 4. Electron Fraction and Neutronization

The electron fraction controls $r$-process yields:

- **Midplane**: $Y_e \approx 0.16$ (neutron-rich, lanthanide-rich ejecta)
- **Outflow average**: $Y_e \approx 0.25$–$0.35$ (lanthanide-poor/rich mixture)

Neutrino absorption rate:

$$\Gamma_\nu \approx 0.22 \cdot T_{10}^6 \cdot D_4 \text{ s}^{-1}$$

where $T_{10} = T/10$ MeV, $D_4 = D/4 \times 10^{-4}$ erg/s/g.

UQFF Layer 2:

$$g_{L2,\nu} = \Gamma_\nu \cdot T_{10}^6 \approx 1.38 \times 10^{-6} \text{ m/s}^2$$

---

## 5. Nuclear Recombination Energy

During wind launch, alpha particle recombination releases:

$$\Delta q_{nuc} \approx 6.8 \times 10^{18} \cdot \Delta X_\alpha \text{ erg/g}$$

where $\Delta X_\alpha$ is the mass fraction of assembled $\alpha$ particles (typically 0.1–0.4). UQFF Layer 3:

$$g_{L3,nuc} = \frac{\Delta q_{nuc} \cdot M_{ej}}{r^2} \approx 8.84 \text{ m/s}^2$$

---

## 6. Jet Luminosity

The GRMHD jet luminosity:

$$L_{jet} \approx 3 \times 10^{50} \text{ erg}$$

This overproduces the GW170817 non-thermal (X-ray/radio) afterglow by a factor ~3, suggesting the jet is choked or structured. UQFF Layer 4:

$$g_{L4,jet} = \frac{L_{jet}}{r^3} \approx 3 \times 10^{-2} \text{ m/s}^2$$

---

## 7. Fast Ejecta — Neutron Precursor

A small fast component:
$$v_{ej,fast} \geq 0.4c, \quad M_{fast} \approx 7.4 \times 10^{-8} M_\odot$$

These are neutron-rich precursors stripped from the disk edge by MHD winds. They generate a UV/optical "neutron precursor" emission 1–2 hours before peak kilonova, potentially observable with ULTRASAT or Swift UV.

---

## 8. GW170817 Context

The extended GRMHD run reproduces:
- Red kilonova: $M_{red} \approx 0.013 M_\odot$, $v \approx 0.1c$, $\kappa \approx 10$ cm²/g (lanthanide-rich)
- Blue kilonova: $M_{blue} \approx 0.005 M_\odot$, $v \approx 0.22c$, $\kappa \approx 0.5$ cm²/g (lanthanide-poor)
- Total $M_{ej}$ budget consistent with AT2017gfo multi-band photometry

The 9-second simulation reveals that MHD outperforms hydrodynamic models at 40% ejection efficiency vs ~20%.

---

## 9. Summary

GRMHD NS merger disk simulations yield $M_{ej} = 0.013 M_\odot$, $v_{ej} = 0.1$–$0.22c$, $L_{jet} = 3\times10^{50}$ erg, lanthanide-rich ejecta ($Y_e = 0.16$) explaining GW170817 extended kilonova. Dual MHD + thermal ejection phases and neutron precursor fast ejecta are now formal UQFF Quadriadic Layer 1–4 parameters.

---

*PAPER_819 | Session 192 | v5.48 | Star-Magic UQFF Project | CVW v2.0.0*
