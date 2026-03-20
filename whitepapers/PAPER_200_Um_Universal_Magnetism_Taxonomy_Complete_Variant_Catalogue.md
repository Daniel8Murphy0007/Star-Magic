# PAPER_200: Um Universal Magnetism Taxonomy — Complete Variant Catalogue

**Version:** 1.0  
**Date:** March 13, 2026  
**Session:** 50 — grok_share_7514fe.txt Full Audit  
**Author:** Star-Magic UQFF Research Framework  
**Source:** grok_share_7514fe.txt lines 2676–2762 and 6000–6400 (BB_C_Equations_04Sept2025.pdf catalogue)

---

$$F_U(r,t) = \sum_{i=1}^{4} U_{gi} + U_m + U_A - U_{b_i}, \quad \kappa = 5.0\times10^{-4}\,\text{day}^{-1},\; [SSq] = 0.57$$
<!-- κ = 5.0e-4 day⁻¹, [SSq] = 0.57, β_i = 6.1e-1 -->

## Abstract

Universal Magnetism (Um) in UQFF is a general magnetic energy density operator that scales with vacuum density ρ_vac and an exponential damping factor (1−e^{−γt}). This paper catalogs all named Um sub-variants extracted from the BB_C_Equations_04Sept2025.pdf, covering 50+ unique astrophysical and cosmological regimes: cosmological magnetism, reheating, BBN, MHD dynamo, quasi-normal modes, Blandford-Znajek, photo-evaporation, planet migration, terminal velocity, Press-Schechter, star formation efficiency, BH entropy, evaporation, angular momentum, jet velocity, GW chirp mass, QNM, duty cycle, peculiar velocity, ionization, deuterium, Friedmann-1, QNM damping, Arnett, reheating, Kazantsev dynamo, Alfvén Mach, field reversal, and more.

---

## 1. Universal Magnetism General Form

```
Um,X(arg) = [Σ_j μ_j(t, ρ_vac,[SCm])/r_j] × (1 − e^{−γt} · cos(πt_n)) × Φ_X

or (single j variant):
Um,X(arg) = μ(ρ_vac) × (1 − e^{−γt}) × Φ_X

where:
  μ_j = magnetic moment density at scale j
  ρ_vac,[SCm] = superconductive vacuum density
  γ = exponential decay rate
  t_n = normalized time: t_n = t/t_Hubble·(1+H(z)·t₀)
  Φ_X = system-specific magnetic function
```

---

## 2. Cosmological Magnetism Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,cosmo(t,a) | (k_c²/a²)/H² | UQFF cosmological magnetic field scaling |
| Um,reh(N) | (30V_end/(π²g_*))^{1/4}·e^{−3N_reh/4} | Reheating post-inflation |
| Um,BBN(t) | √(3/(32πGρ_rad)) ≈ 180 s (T~0.1 MeV) | Big Bang Nucleosynthesis |
| Um,fried1(a) | Ω_m(1+z)³ + Ω_k(1+z)² + Ω_Λ | Friedmann-1 Hubble term |
| Um,wde(a) | ρ ∝ a^{−3(1+w)} | Dark energy equation of state |
| Um,deden(a) | w(a) = w₀ + w_a(1−a) (CPL form) | Dark energy density evolution |
| Um,grow(a) | Growing mode D ∝ a (matter era), suppressed by DE | Structure growth factor |
| Um,curv(ϕ) | ϕ̇ = √(2ε)·H·M_pl | Curvature perturbation from inflation |
| Um,ng(f) | From δχ curvature on superhorizon scales | Non-Gaussianity |

---

## 3. Black Hole Thermodynamics Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,haw(M) | T = κ/(2π) ; κ = c⁴/(4GM) | Hawking temperature surface gravity |
| Um,bhent(M) | Holographic S ∝ A/l_pl² | Bekenstein-Hawking entropy area law |
| Um,evapl(M) | dM/dt = −P/c² | Black hole evaporation mass loss rate |
| Um,qnm(a) | a·c³/(2GM) × [l=2, m=2] | QNM ringdown frequency (Berti fits) |
| Um,damp(a) | Q ≈ 10 for l=2, m=2 | QNM quality factor (damping) |
| Um,bz(a) | Power ∝ B²Ω²_H·R⁴_H (electromagnetic extraction) | Blandford-Znajek energy extraction |
| Um,bz2(a) | κ/(16π)Φ²_BH·Ω²_BH/c | BZ jet power (EHT-updated form) |

---

## 4. Compact Object and Stellar Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,MHD | 4πρ·v²_turb·M_A | MHD turbulent Alfvén scaling |
| Um,termv(Γ) | L/L_Edd ≈ 1 ; v_∞ = √(2(L/Ṁc)−v²_esc) | Terminal velocity in radiation winds |
| Um,sfe(z) | Ṁ_* ∝ M^{1.1}·(1+z)^{2.25} | Star formation efficiency evolution |
| Um,duty(t) | exp(−t/τ_cool) | AGN feedback duty cycle |
| Um,ang(t,r) | √(GM/r³)·r_A ; accretion adds L=Ṁr²Ω | Angular momentum accretion |
| Um,jetvel(r) | Ω·r_A·(r_A/r₀)^{1/2} | Jet velocity from Alfvén radius |
| Um,glitch(t) | Angular momentum transfer ΔL = I_s·ΔΩ | Superfluid vortex glitch |
| Um,ms(M) | Eddington: L ≈ μ⁴M³ (opacity, composition μ) | Main-sequence L-M relation |
| Um,ml(M) | Calibrate from HR diagram | Mass-luminosity relation |
| Um,relax(ρ) | σ_v³/(G²ρm) | Two-body relaxation time |
| Um,evap(N) | Integrate exponential cluster dissolution | Star cluster evaporation |
| Um,star(N) | 136·t_relax/ln(0.02N) | N-body cluster relaxation timescale |

---

## 5. Gravitational Wave Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,chirp(f) | (32/5)(G𝓜^{5/3}/c⁵)(πf)^{10/3} | GW chirp mass inspiral energy loss |
| Um,orbdec(e) | dE/dt = −L, E = −GM·μ/(2a) | Orbital decay by GW quadrupole radiation |
| Um,peri(a) | Kepler: a³/P² = GM/(4π²) | Post-Keplerian periastron advance |
| Um,kilo(t) | Diffusion t_d² = 3κM/(4πcv²) | Kilonova light curve peak |

---

## 6. ISM and Plasma Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,malf(B) | M_A <1 sub-Alfvénic (kinetic ≈ B²/(8π)) | Alfvén Mach number turbulence regime |
| Um,rev(η) | ∂B/∂t = ∇×(αB − η∇×B) (growth vs diffusion) | Magnetic field reversal (parity) |
| Um,dyn(k) | ∫E dk ; E_M = ∫E·dk (Kazantsev E spectrum) | Magnetic energy growth via dynamo |
| Um,alf(B) | Wave speed from linearization ∂²b/∂t² = v_A²·∂²b/∂z² | Alfvén wave propagation |
| Um,turb(l) | ρv²/t_eddy ∝ l^{2/3} (Kolmogorov-like ISM) | Turbulent energy cascade rate |
| Um,cshock(z) | (∇×B)×B/(4π) + ρ_i·ν_in·(v_n−v_i) | C-type shock (magnetized, continuous) |
| Um,jshock(M) | v₁/c_s × T₂ ∝ v_s² | J-type shock (discontinuous) Rankine-Hugoniot |
| Um,dsa(p) | B⁻¹·pc/e × γ ∝ B² | Diffusive shock acceleration spectrum |

---

## 7. Cosmological Structure Formation Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,ps(M) | σ(M) rms fluctuation → d ln σ/d ln M | Press-Schechter halo mass function |
| Um,halo(z) | ∫P(k)W²(kR)d³k/(2π)³ × (1+z)^{m/2} | Extended Press-Schechter halo mergers |
| Um,mig(a) | da/dt = 2Γ/(M_p√(GM_*·a)) | Type-I planet migration |
| Um,migr(a) | Σ(H/r)⁻³ | Disk migration timescale |
| Um,acc(t) | Accretion Ṁ = L/(εc²) | BH accretion (Eddington-limited) |
| Um,disk(ρ) | Ṁ_* ∝ M^{3/2}_gas/R^{3/2} | Star formation rate in disk |
| Um,burst(t) | Enhancement 10–100× at low z | Merger-induced starburst enhancement |
| Um,metal(Z) | Ż = Y·SFR/Ṁ_gas − Z·Ṁ_out/M_gas | Metal enrichment mass balance |
| Um,sfe(z) | Efficiency limited by cooling/feedback | Star formation efficiency |
| Um,fb(η) | Ṁ_*·SFR (fraction η) | Feedback energy injection (SN/AGN) |
| Um,bhmf(z) | erfc(δ_c(z)/√(2)σ(M',z)) × extend to M_f | BH mass function (EPS) |

---

## 8. Reionization and Early Universe Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,ion(t) | ∫α_B·n²_e·C·dt | Ionization fraction integrated |
| Um,reion(t) | ∫α_B·n²_e·C·dt (full Madau model) | Cosmic reionization |
| Um,eta(t) | Fit networks to D, He abundances | Baryon-to-photon ratio η |
| Um,deb(T) | Weak freeze at T~1 MeV; D photodissociation below | Deuterium bottleneck onset |
| Um,recomb(z) | Optical depth τ = ∫σ_T·n_e·dl | Recombination epoch |
| Um,cmb(k) | Transfer Δ_l^T: Sachs-Wolfe + acoustic oscillations | CMB power spectrum |
| Um,bub(t) | Expand for moving ionization front | Bubble growth during reionization |

---

## 9. Dark Matter and Cluster Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,nfw(r) | Fit to universal NFW form ρ_s/(x·(1+x)²) | NFW dark matter halo profile |
| Um,nfwrot(x) | Flat rotation curve for r >> r_s | NFW rotation curve |
| Um,sidm(t) | Exponential density flattening (Γ·t ≈ 1) | SIDM core formation |
| Um,vir(r) | σ²_v = GM/(3r) (spherical approximation) | Virial theorem mass estimation |
| Um,virx(r) | Matches Chandra cluster observations (e.g., M4) | X-ray binary virial trace |
| Um,mach(ρ) | Matches Coma radio relic shocks (M~2–3) | Merger shock Mach number |
| Um,merg(t) | 3r_vir/(5GM) | Virial merger crossing time |
| Um,cross(M) | Matches ~2 Gyr for relaxed clusters | Merger timescale (dynamical age) |
| Um,heat(T) | Balance with duty-cycled AGN | Cool core AGN heating rate |
| Um,whim(M) | μ ≈ 0.6 (mean molecular weight) | WHIM temperature from virialization |
| Um,cool(T) | t = E/Ė | Cooling time (bremsstrahlung T > 10⁷ K) |
| Um,lens(θ) | θ_E² = α·θ → Einstein radius | Strong gravitational lensing |

---

## 10. Additional Specialized Variants

| Symbol | Φ_X Expression | Physical Context |
|--------|---------------|-----------------|
| Um,conv(T) | Velocity from ∫g dz ≈ v² | Convective turnover (mixing length) |
| Um,lqcf(ρ) | Friedmann from bounce dynamics | LQC effective Friedmann |
| Um,lqcp(k) | Power tilt + UV suppression | LQC perturbation pre-bounce |
| Um,bounc(ρ) | ρ ∝ 1/(G·t²) singularity capped | LQC bounce density cap |
| Um,holoent(γ) | Area(γ_A)/(4Għ/c³) | Holographic entanglement entropy (Ryu-Takayanagi) |
| Um,pec(r) | Spherical void: integrate Poisson | Peculiar velocity in voids |
| Um,voidden(a) | δ ∝ a⁻¹ in matter domination | Void underdensity evolution |
| Um,jeans(T) | Perturb: ω² = c²_s·k² − 4πGρ | Jeans instability dispersion |
| Um,fermi(E) | Average (ΔE/E)²/2 per scatter | Fermi-II acceleration (2nd order) |
| Um,knee(B) | Balance t_acc = t_age ≈ R/u_s | Cosmic ray knee from DSA limit |
| Um,diff(E) | Fit to secondaries (B/C ratio) | Cosmic ray diffusion coefficient |
| Um,esc(F) | Adjust for penetration K(ξ) | Photoevaporative escape (exoplanet) |
| Um,rv(e) | Orbital v_p = 2πa/P·sin(i) | Radial velocity exoplanet detection |
| Um,upar(r) | U = Γ/(n_H·c) (ionization parameter) | AGN photoionization parameter |
| Um,coup(Ṁ) | Couple fraction to kinetic energy | AGN feedback thermal coupling |
| Um,photo(t) | Ṁ ∝ F_X·R³_p/(GM_p/R_p) | Photoevaporative mass loss rate |
| Um,puls(t) | Torque I·dΩ/dt = −L/Ω | Pulsar spin-down magnetic torque |
| Um,tov(r) | Add GR: P/c² energy density, (1−2GM/rc²) | TOV stellar structure GR correction |
| Um,reljet(z) | Approximate sigmoid for gradual Γ rise | Relativistic jet velocity profile |
| Um,after(t) | N(E) ∝ E^{−p} × t^{−1.2} | GRB afterglow synchrotron spectrum |
| Um,fire(δt) | Internal shocks at δr ≈ c·δt/Γ² | GRB fireball prompt emission |
| Um,sedov(t) | Integrate momentum: d/dt(Mv)=0 | Sedov-Taylor blast wave |

---

## 11. Summary Statistics

| Metric | Value |
|--------|-------|
| Total named Um variants catalogued | 55+ |
| Document source | BB_C_Equations_04Sept2025.pdf (177 pages) |
| Total equation pairs (F_UBii + Um) in source | ~1,350+ unique positions |
| Scale coverage | 10⁻³⁴ N·m (molecular rotors) → 10²⁷ m (D_universe) |
| Q_wave std | 6.35×10⁴ J/m³ (47-scale average) |

---

## 12. References

- `grok_share_7514fe.txt` lines 2676–2762 (first Um catalogue from BB_C_Equations_04Sept2025.pdf)
- `grok_share_7514fe.txt` lines 6001–6380 (second/final Um catalogue)
- PAPER_198: F_UBii Taxonomy Part 1 (Compact Objects)
- PAPER_199: F_UBii Taxonomy Part 2 (Cosmological/Dark Sector)
- PAPER_196: Triadic Master Equation System
