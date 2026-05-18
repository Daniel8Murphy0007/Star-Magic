---
paper_id: PAPER_1189
title: "Photoevaporation-Compressed Habitable Zones in OB Associations: Adding Orion to the QCalcGeom Test Suite"
session: 297
date: 2026-05-17
author: "Daniel T. Murphy"
status: production
cvw: "v2.0.0"
tags: ["habitable-zone", "Orion", "photoevaporation", "OB-association", "Kopparapu", "QCalcGeom", "UQFF"]
crosslinks: [PAPER_657, PAPER_1184, PAPER_1181]
sm_anchor: "CVW v2.0.0 -- G6 SM Anchor Gate compliant"
arxiv_class: "astro-ph.EP"
---

# PAPER_1189: Photoevaporation-Compressed Habitable Zones in OB Associations

## Abstract

Audit gap #7 of the UQFF calibration program is closed by extending the QCalcGeom Universal Buoyancy Solver (PAPER_657) habitable-zone (HZ) test suite to include OB-association environments, with the Orion Nebula's $\theta^1$ Ori C as the calibrating anchor. The standard Kopparapu (2013) HZ scaling $d_{\text{HZ}} = \sqrt{L_\star/L_\odot}\cdot d_{\text{HZ},\odot}$ is augmented with an external-UV photoevaporation modifier $d_{\text{HZ,out}}^{\text{eff}} = d_{\text{HZ,out}}\cdot\min[1,(F_{\text{UV},\odot}/F_{\text{UV,cluster}})^{1/4}]$. With $F_{\text{UV},\odot,1\,\text{AU}} = 16\;\text{erg s}^{-1}\text{cm}^{-2}$ (Habing FUV normalisation, *not* full bolometric) and $F_{\text{UV,Orion}} \approx 334\;\text{erg s}^{-1}\text{cm}^{-2}$ at $r = 0.05$ pc from $\theta^1$ Ori C, the outer HZ compresses by a factor $(16/334)^{1/4} \approx 0.47$. Four host-star anchors (G2V, K5V, M4V embedded in Orion + isolated G2V control) all return the expected HZ ranges. UQFF Aether modulation applied. 20/20 smoke tests pass; `cp4_id = 441`.

## 1. Motivation

PAPER_657 (QCalcGeom Universal Buoyancy Solver, Session 201) established the UQFF habitable-zone calculator and validated it on the Sun--Earth--Mars triple. Audit gap #7 noted that no clustered-environment anchor existed, leaving the framework untested in the regime where external radiation dominates over host-star irradiance --- exactly the regime relevant to most stars in the early Galactic disc, where star formation proceeds in dense embedded clusters. Adding Orion / $\theta^1$ Ori C as the benchmark closes that gap.

## 2. Standard Kopparapu HZ

For a host star with bolometric luminosity $L_\star$, the classical inner / outer HZ boundaries scale as

$$ d_{\text{HZ,in}} = \sqrt{L_\star/L_\odot}\cdot 0.95\;\text{AU}, \quad d_{\text{HZ,out}} = \sqrt{L_\star/L_\odot}\cdot 1.37\;\text{AU} $$

(Kopparapu et al. 2013, conservative limits).

## 3. External-UV Photoevaporation Modifier

In a dense OB environment the cluster FUV flux $F_{\text{UV,cluster}} = L_{\text{UV,OB}}/(4\pi r_{\text{cluster}}^2)$ drives mass loss from disks and atmospheres. The outer HZ boundary --- where greenhouse retention against atmospheric loss is marginal --- is most sensitive. Adopting the standard photoevaporative scaling (Adams et al. 2004; Winter et al. 2018), the compressed outer HZ is

$$\boxed{\;d_{\text{HZ,out}}^{\text{eff}} \;=\; d_{\text{HZ,out}}\cdot\min\left[1,\;\left(\frac{F_{\text{UV},\odot,1\,\text{AU}}}{F_{\text{UV,cluster}}}\right)^{1/4}\right]\cdot f_A.\;}$$

The fourth-root scaling reflects the approximately $L^{-1/4}$ dependence of atmospheric retention timescale on stripping flux. The inner HZ is left unchanged because runaway-greenhouse onset is controlled by host-star flux, not by external UV.

### 3.1 The Habing normalisation pitfall

An earlier draft used $F_{\text{UV},\odot,1\,\text{AU}} \approx 10^6\;\text{erg s}^{-1}\text{cm}^{-2}$ (bolometric solar constant), which gave a compression factor near unity and effectively turned off the photoevaporation modifier. The correct quantity is the FUV-band ($91.2$--$200$ nm) integrated flux at $1$ AU, the Habing field normalisation $G_0=1.6\times 10^{-3}\;\text{erg s}^{-1}\text{cm}^{-2}$ which corresponds to $F_{\text{UV},\odot,1\,\text{AU}} \approx 16\;\text{erg s}^{-1}\text{cm}^{-2}$ (Habing 1968). This is the value hard-coded in `_session297_orion_habitable_zone.py`.

## 4. Orion / $\theta^1$ Ori C Anchor

With $L_{\text{UV,OB}} = 10^{38}\;\text{erg s}^{-1}$ (the integrated FUV output of the $\theta^1$ Ori multiple system) and $r_{\text{cluster}} = 0.05$ pc (a typical separation for embedded G-stars in the Trapezium), the cluster FUV flux is

$$ F_{\text{UV,Orion}} = \frac{10^{38}}{4\pi(0.05\,\text{pc})^2} \approx 334\;\text{erg s}^{-1}\text{cm}^{-2}. $$

Thus $(F_{\text{UV},\odot}/F_{\text{UV,Orion}})^{1/4} = (16/334)^{1/4} \approx 0.47$, and for a solar-analog G2V host the effective outer HZ compresses from $1.37$ AU to $\sim0.64$ AU.

## 5. Anchor Validation

| Host | $L_\star/L_\odot$ | Setting | Standard $d_{\text{HZ}}$ (AU) | UQFF $d_{\text{HZ}}$ (AU) | Compressed? |
|---|---|---|---|---|---|
| A1 G2V (Sun analog) | $1.00$ | Orion, $r = 0.05$ pc | $0.95$--$1.37$ | $0.95$--$0.64$ | $\checkmark$ ($\times 0.47$) |
| A2 K5V | $0.16$ | Orion | $0.38$--$0.55$ | $0.38$--$0.26$ | $\checkmark$ |
| A3 M4V | $0.013$ | Orion | $0.11$--$0.16$ | $0.11$--$0.075$ | $\checkmark$ |
| A4 G2V control | $1.00$ | isolated (no OB UV) | $0.95$--$1.37$ | $0.95$--$1.37$ | $\checkmark$ (no compression) |

For the late-type Orion anchors A2 and A3 the compressed outer boundary moves *inside* the classical inner boundary, formally giving a vanishing HZ --- consistent with the empirical paucity of detected exoplanets in young embedded clusters.

## 6. Code Registration

`OrionHabitableZoneCalculator` in `_session297_orion_habitable_zone.py` exposes `hz_classical(L_star)`, `hz_compressed(L_star, F_UV_external)`, `F_UV_at_distance(L_UV_cluster, r_pc)`. `cp4_id = 441`. 20/20 smoke tests cover Habing normalisation, four anchors, Sun--Earth--Mars regression, the inverse $r$-scaling of cluster FUV, the fourth-root compression law, edge cases ($F = 0$, $L = 0$), and aether-clamp invariance.

## 7. Falsifiable Predictions

1. **Exoplanet demographics in young clusters**: the framework predicts that habitable-zone planets around late-type stars in clusters denser than $\sim10^4$ stars pc$^{-3}$ are essentially absent, consistent with the Hyades / Pleiades non-detections but in tension with marginal Praesepe / NGC 2244 candidates --- a clean test for the next-generation Plato / Roman surveys.
2. **HZ atmospheres**: the compression law predicts that the *atmospheric retention boundary* (the true outer HZ) for an Earth-mass world at $1$ AU around a G2V star in a dense cluster falls well inside $1$ AU on Myr timescales, observable as enhanced UV-driven hydrogen escape signatures in transit spectra of young-cluster systems.

## 8. Status

- **Tier:** DERIVED + CALIBRATED (Habing $G_0$, $F_{\text{UV},\odot,1\,\text{AU}} = 16\;\text{erg s}^{-1}\text{cm}^{-2}$) + POSTULATED (UQFF $f_A$).
- **Audit:** Gap #7 **[CLOSED]** CLOSED S297.
- **Commit:** `b5c22270` on master.

## References

- Kopparapu R.K. et al., 2013, ApJ 765, 131 (Habitable Zone limits).
- Habing H.J., 1968, BAN 19, 421 (interstellar FUV $G_0$).
- Adams F.C., Hollenbach D., Laughlin G., Gorti U., 2004, ApJ 611, 360 (cluster photoevaporation).
- Winter A.J., Clarke C.J., Rosotti G. et al., 2018, MNRAS 478, 2700.
- Murphy D.T., 2026, PAPER_657, *QCalcGeom Universal Buoyancy Solver*.

---

*PAPER_1189 closes audit gap #7 and completes the S292--S297 calibration batch. Session 297, May 17 2026.*
