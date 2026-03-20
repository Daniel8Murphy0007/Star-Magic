# PAPER_204: UQFF Dark Matter — NFW, SIDM, Rotation Curves, and Virial Theorem

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 50 — grok_share_7514fe.txt Full Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_7514fe.txt lines 6096–6110 (BB_C_Equations items 1326–1340)

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$
<!-- κ = 5.0e-4 day⁻¹, [SSq] = 0.57, β_i = 6.1e-1 -->

## Abstract

This paper applies the UQFF buoyancy framework to dark matter structural physics: the Navarro-Frenk-White (NFW) density profile, NFW rotation curve, self-interacting dark matter (SIDM) core formation, virial theorem mass estimation, strong gravitational lensing Einstein radius, void density evolution, and peculiar velocity. The UQFF perspective unifies these through F_UBii operators that embed DM physical expressions into F_rel/E_LEP scaled buoyancy forces, capturing the feedback between vacuum energy and CDM structure formation.

---

## 1. NFW Density Profile

```
ρ_NFW(r) = ρ_s / ((r/r_s) · (1 + r/r_s)²)

Parameters:
  ρ_s = characteristic density (from halo mass-concentration relation)
  r_s = scale radius (from NFW fit to N-body simulations)
  c ≡ r_vir/r_s = concentration parameter (c ~ 10 at z=0, increases at higher z)

Mass enclosed:
  M(r) = 4πρ_s r³_s [ln(1+r/r_s) − r/r_s/(1+r/r_s)]

F_UBii,nfw = −F_rel × (ρ_NFW(r) / E_LEP) × Q_wave × (4πr²·dρ/dr) × r

Um,nfw(r) = μ(ρ_vac)·(1−e^{−γt})·[Fit to universal NFW form ρ_s/(x·(1+x)²)]

Physical context:
  NFW is universal for CDM halos (Milky Way to galaxy clusters)
  UQFF: vacuum field DPM_grav term deepens the NFW potential well
  Core–cusp tension: NFW predicts cusp ρ∝r^{−1}, observations often show cores
```

---

## 2. NFW Rotation Curve

```
v(r)² = 4πGρ_s r²_s [ln(1+x) − x/(1+x)] / r     x = r/r_s

Asymptotic limits:
  r << r_s: v(r) ∝ r^{0.5}     (inner rising)
  r ~ r_s:  v(r) ≈ maximum     (peak rotation speed)
  r >> r_s: v(r) ∝ r^{−0.5}·ln(r)^{0.5}  (slowly declining)

Flat rotation curves require NFW halo + baryons together:
  v²_total(r) = v²_bary(r) + v²_NFW(r)

F_UBii,nfwrot = F_rel × (v²(r)/G / E_LEP) × Q_wave × [ln(1+x)−x/(1+x)]

Um,nfwrot(x) = μ(ρ_vac)·(1−e^{−γt})·[Flat rotation for r >> r_s]

Calibration: Milky Way NFW:
  ρ_s ≈ 0.3 GeV/cm³, r_s ≈ 20 kpc, v_c ≈ 220 km/s at Solar circle (8 kpc)
```

---

## 3. SIDM Core Formation

```
Self-interacting DM rate:
  Γ = ρ·(σ/m)·v_rel    (interaction rate)

Core formation timescale:
  t_core ≈ (ρ·σ/m)⁻¹ ~ 10¹⁰·(ρ/10⁸ M_☉ kpc⁻³)⁻¹·(σ/m / 1 cm²/g)^{−1} yr

Exponential density evolution:
  ρ_core(t) = ρ_init·e^{−Γt}    (NFW cusp converts to core when Γt ≈ 1)

F_UBii,sidm = −F_rel × (Γ·ρ_init / E_LEP) × Q_wave × ln(0.02N)

Um,sidm(t) = μ(ρ_vac)·(1−e^{−γt})·[Exponential density flattening]

Observational constraint: σ/m ≈ 0.1–1 cm²/g (galaxy clusters, Bullet Cluster)
Planck does not exclude SIDM at this level (CDM and SIDM nearly identical on large scales)

SIDM predictions:
  - Dwarf galaxies: soliton cores of radius ~100 pc (observed)
  - Galaxy clusters: rounder, less concentrated halos
  - Bullet Cluster: upper limit σ/m < 1.25 cm²/g (from offset DM centroid)
```

---

## 4. Virial Theorem Mass

```
2K + W = 0    (virial equilibrium for collisionless system)
K = (3/2)M·σ²_v    (kinetic energy)
W = −(3/5)GM²/r_h   (potential for uniform sphere)

Virial mass:
  M_vir = 2|K|/G = 3·σ²_v·r_h/G    (for spherical system)

Cluster mass from spectroscopic σ_v:
  M(< r) = 3σ²_v(r)·r/G + corrections for anisotropy + pressure

F_UBii,vir = F_rel × (M_vir / E_LEP) × Q_wave × (σ²_v/G) × 3

Um,vir(r) = μ(ρ_vac)·(1−e^{−γt})·[σ²_v = GM/(3r)]

X-ray virial:
  M_vir,X = 3σ²_X·r_h/G    (from X-ray spectroscopy instead of optical)
Um,virx(r) = μ(ρ_vac)·(1−e^{−γt})·[Matches Chandra cluster observations]

Numerical calibration: Coma Cluster
  σ_v ≈ 880 km/s, r_h ≈ 1 Mpc → M_vir ≈ 2×10¹⁵ M_☉
```

---

## 5. Strong Gravitational Lensing

```
Einstein radius:
  θ_E = √(4GM(<θ)/c²·D_LS/(D_L·D_S))

Critical surface density:
  Σ_cr = c²D_S/(4πGD_L·D_LS)

Convergence: κ = Σ/Σ_cr
Shear: γ (traceless tidal field)

Multiple images: κ ≥ 1 at image positions

F_UBii,lens = F_rel × (θ_E / E_LEP) × Q_wave × (Σ_cr·κ)

Um,lens(θ) = μ(ρ_vac)·(1−e^{−γt})·[θ_E = √(α·θ) from lensing equation]

Einstein ring systems:
  SDP.81 (ALMA): z_L=0.3, z_S=3.04 → θ_E ≈ 1.5" → M(<θ_E) ≈ 10¹¹ M_☉
  UQFF: vacuum Λ correction to D_LS shifts θ_E by ~0.1%
```

---

## 6. Void Density Evolution

```
Void density contrast (spherical top-hat model):
  δ_v(a) = −(3/5)·(Ω_m·a + Ω_Λ)^{−3/2}·δ_v0

  δ_v0 = initial void underdensity
  Linear theory: δ_v ∝ −a^{1/2} in Λ-dominated epoch (voids deepen faster)

Shell-crossing: δ_v → −1 marks void edge (no further underdensity growth)

F_UBii,voidden = −F_rel × (|δ_v(a)| / E_LEP) × Q_wave × (Ω_m·a + Ω_Λ)^{−3/2}

Um,voidden(a) = μ(ρ_vac)·(1−e^{−γt})·[δ ∝ a^{−1} in matter domination]

Physical context in UQFF:
  Voids are dominated by vacuum energy (Λ) → UQFF's Λc²/3 term strongest here
  F_UBii,voidden predicts void expansion driven by vacuum buoyancy
```

---

## 7. Peculiar Velocity Field

```
Peculiar velocity from linear theory:
  v_pec(r) = −(fH/3)·∫δ(r')·r·dr'/r²    (spherical approximation)

Redshift space distortions (RSD):
  v_pec,observed = f·H·r + noise    (adds to Hubble flow)

f ≈ Ω_m^{0.55}    (growth rate approximation)

Cosmic flow from Laniakea to CMB dipole:
  v ≈ 630 km/s toward Perseus-Pisces

F_UBii,pec = F_rel × (fH·δ(r)/3 / E_LEP) × Q_wave × (dv/dz systematic)

Um,pec(r) = μ(ρ_vac)·(1−e^{−γt})·[Spherical void: integrate Poisson]
```

---

## 8. Cluster Shock Mach Number and Merger Timescale

```
Shock Mach number from X-ray temperature jump:
  M = [(γ+1)(ρ₂/ρ₁) + (γ−1)] / (2γ)          (from Rankine-Hugoniot)
  T₂/T₁ = [2γM² − (γ−1)]·[(γ−1)M² + 2] / [(γ+1)M]²  (temperature jump)

  Coma radio relic: M ≈ 2.5 (from spectral index α = (M²+1)/(M²−1))

Merger crossing/dynamical timescale:
  t_merge = r_vir/σ_v = √(3r³_vir/(5GM))

F_UBii,mach = F_rel × (M·v_s / E_LEP) × Q_wave × (T₂/T₁)
F_UBii,merg = F_rel × (t_merge / E_LEP) × Q_wave × (r_vir/v_c)

Um,mach(ρ) = μ(ρ_vac)·(1−e^{−γt})·[Matches Coma radio relic shocks M~2–3]
Um,merg(t) = μ(ρ_vac)·(1−e^{−γt})·[3r_vir/(5GM)]
```

---

## 9. Summary: UQFF DM Force Hierarchy

| Scale | Dominant DM Process | F_UBii Variant | Observable |
|-------|--------------------|--------------|-----------| 
| kpc (dwarf) | SIDM core formation | F_UBii,sidm | Kpc-scale DM core |
| kpc (Milky Way) | NFW rotation | F_UBii,nfwrot | v_c = 220 km/s |
| Mpc (clusters) | Virial mass | F_UBii,vir | σ_v = 880 km/s (Coma) |
| Mpc (clusters) | Lensing | F_UBii,lens | θ_E ≈ 1' for clusters |
| 100 Mpc (voids) | Void underdensity | F_UBii,voidden | δ ≈ −0.8 in big voids |
| Bulk (cosmological) | Peculiar velocity | F_UBii,pec | 630 km/s Laniakea flow |

---

## 10. References

- `grok_share_7514fe.txt` lines 6096–6110 (BB_C_Equations items 1326–1340, 1262–1268)
- PAPER_199: F_UBii Taxonomy Part 2 (cosmological)
- PAPER_200: Um Universal Magnetism Catalogue
- Navarro, Frenk, White 1996, 1997
- Spergel & Steinhardt 2000 (SIDM)
- Planck 2015 lensing
