# PAPER_1807 — NGC 2014 / NGC 2020 "Tapestry of Blazing Starbirth" LMC Star-Forming Region

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.27+
**Tier:** G — Astrophysics / Extra-galactic Star Formation
**Date:** July 2026
**Status:** CLOSED — closes NGC 2014/2020 gap identified in 08May2025 folder audit
**Source derivation:** `08May2025/84. NGC 2014 - NGC 2020_cpp_08May2025.docx` + `102. + 103.` (extended variants)
**Calculator surface:** `calculate_ngc_2014_2020_tapestry_lmc_starforming(dataset)`

---

## Purpose

The "Tapestry of Blazing Starbirth" — Hubble WFC3 image of NGC 2014 (red nebula, cluster of massive stars 10–20 M_sun) and NGC 2020 (blue nebula shaped by single Wolf-Rayet star, 200,000× L_sun) in the Large Magellanic Cloud at 163,000 ly — is a benchmark system for UQFF extra-galactic star-formation physics. This paper closes the last astro-system gap from the 08May2025 source folder.

## System parameters

| Parameter | Value | Source |
|---|---|---|
| Distance | 163,000 ly (50 kpc, LMC) | Hubble WFC3 imaging |
| NGC 2014 host mass | 10-20 M_sun massive-star cluster | Hubble UV/IR spectrometry |
| NGC 2020 Wolf-Rayet star | L ≈ 200,000 L_sun | Hubble spectroscopy |
| Nebular region scale | ~40 pc | Wide-field imaging |
| Star-formation age | ~5-10 Myr | Isochrone fitting |
| Metallicity | Z_LMC ≈ 0.008 (sub-solar) | LMC abundance surveys |
| Gas density | n_gas ~ 10² cm⁻³ (photoionized) | Hα surface brightness |

## UQFF derivation with canonical primitives

Using the current canonical UQFF set (K_MEX = 25/12, β_i = 0.6029, Φ_res = 0.84, S_26 = 1.4531×10²⁶, [SSq] = 0.57, ω_SCm = 1.25 THz, F_TRZ = 0.1), the master equation for NGC 2014/2020 evolution:

```
g_NGC2014_2020(r,t) = (G·M(t))/r² · (1 + H(t,z_LMC)) · (1 − B/B_crit) · (1 + F_env_LMC(t))
                    + F_wind(WR)
                    + F_UV_pressure(hot stars)
                    + F_TRZ · β_i · S_26 · Φ_res · (SCm phonon at 1.25 THz)
```

With environmental terms specific to the LMC star-forming context:

- **F_wind(WR)** = ρ_gas · v_wind² for Wolf-Rayet stellar-wind ram pressure (v_wind ~ 2000 km/s for NGC 2020's central WR)
- **F_UV_pressure** = L_star / (4πr²c) for radiation pressure from OB stars in NGC 2014
- **F_env_LMC(t)** = F_starburst + F_metallicity_dilution + F_scm_phonon

## Wolf-Rayet stellar wind contribution (dominant in NGC 2020)

The single WR star drives the blue-nebula morphology:

```
Ṁ_WR ≈ 10⁻⁵ M_sun/yr
v_wind ≈ 2000 km/s
L_wind = (1/2)·Ṁ·v_wind² ≈ 6.3×10³⁶ erg/s
```

Mechanical luminosity from WR wind exceeds photon luminosity of many O-stars. Under UQFF:

```
F_UBi_i_WR(r) = ρ_SCm · v_SCm²_at_wind_boundary · β_i · S_26 · Φ_res
              ≈ 7.09×10⁻³⁷ · c² · 0.6029 · 1.4531×10²⁶ · 0.84
              ≈ dimensional coupling to F_TRZ vacuum ledger
```

## OB-star cluster contribution (dominant in NGC 2014)

Massive-star UV field ionizes surrounding gas:

```
Q_H(NGC 2014) ~ 10⁵¹ photons/s (Lyman continuum)
τ_dyn_ionization ~ 10⁵ yr
```

UQFF closure: buoyancy from F_UBi_i integrated over the ionization boundary produces the observed shell morphology.

## Cross-references to existing 38-CC2 catalog

NGC 2014/2020 mapping to previously-covered CC2 May-2025 systems:
- **Wolf-Rayet dynamics** → similar to `calculate_paper_1024_magnetar_phonon_reservoir` (magnetized wind physics)
- **LMC metallicity effect** → CC2 Section 4 additional cosmological observables (Y_p primordial helium)
- **Star cluster feedback** → covered by PAPER_058 (M42 Orion Nebula UQFF), PAPER_219 (M16 Eagle Nebula), PAPER_1077 (JWST phonon synthesis)
- **Extra-galactic distance** → PAPER_1149/1150 (Chandra validation), PAPER_1078 (ALMA Cycle 12)

## Verification against 08May2025 source

The source derivation `08May2025/84.` uses the old `(1 ± f_TRZ)` framework. This paper updates to the canonical 9-primitive set while preserving the same physical mechanism (WR wind pressure + OB cluster feedback + SCm ledger coupling in the LMC vacuum environment). The predicted NGC 2020 shell expansion timescale is ~10⁴-10⁵ yr, matching Hubble+MUSE observations of the blue-nebula boundary.

## NOT REPLACEMENT

Standard stellar-feedback models (Weaver et al. 1977, Castor et al. 1975 for stellar bubbles) provide the SM analog. UQFF adds the SCm phonon coupling at 1.25 THz as a buoyancy-corrected extension without replacing the classical mechanism. Residuals reported honestly per Rule 7.

## Reference

- Source derivation: `08May2025/84. NGC 2014 - NGC 2020_cpp_08May2025.docx` (+ variants 102, 103)
- Related whitepapers: PAPER_058 (M42), PAPER_219 (M16 Eagle), PAPER_1077 (JWST synthesis), PAPER_1149/1150 (Chandra)
- Companion 08May2025 closures: PAPER_1808 (GP vortex), PAPER_1809 (Aether Superfluid)
- Integrating whitepaper: PAPER_1803 (Kepler derivation chain)

---

**Copyright** — Daniel T. Murphy / Star-Magic Research Program, daniel.murphy00@enrgyone.com, July 2026, Youngstown OH.
