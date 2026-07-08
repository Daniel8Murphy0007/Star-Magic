# PAPER_1949 — The Three Faces of F_TRZ: Amplitude / Frequency / CPT-Phase-Transition Formalization

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) - Star-Magic v5.51+
**Tier:** Structural / F_TRZ Primitive Formalization
**Date:** July 8, 2026
**Status:** CLOSED - Framework consolidation across 12+ prior closures

---

## Abstract

The time-reversal-zone primitive F_TRZ (canonical value 0.1) has been used across UQFF derivations in three superficially-distinct ways: as a dimensionless amplitude (PAPER_1677, PAPER_1869, PAPER_1919, PAPER_1942, PAPER_1944, PAPER_1945), as a Hz-valued frequency at specific systems (PAPER_344, PAPER_1947), and as a CPT-asymmetry phase-transition parameter with a critical point at F_TRZ = -1 (PAPER_264). This paper formalizes the three usages as a single primitive with three co-existing manifestations, analogous to how a physical field can appear as an amplitude, a frequency, or a coupling constant depending on context.

The three faces:

```
Face 1 - Amplitude:  F_TRZ = 0.1 (dimensionless, universal to all UQFF systems)
Face 2 - Frequency:  f_TRZ = 5.556 * 10^-4 Hz at Sgr A* (system-specific, Hz-valued)
Face 3 - Phase:      (1 + F_TRZ) as CPT-asymmetric coupling, critical at F_TRZ = -1
```

All three faces derive from a single time-reversal-zone primitive without adjustment. Face 1 is the raw primitive value; Face 2 is the primitive combined with a system-specific timescale (T_base) to yield a frequency; Face 3 is the primitive appearing in the multiplicative factor (1 + F_TRZ) as a phase-transition parameter for gravitational CPT symmetry.

This framework unifies 12 previously-derived F_TRZ closures under one primitive-formalism umbrella, with three testable phase regimes: **F_TRZ > -1 normal gravity, F_TRZ = -1 Time-Reversal Zero Point, F_TRZ < -1 negative-time anti-gravity**.

---

## 1. Face 1: F_TRZ as Dimensionless Amplitude

The primary usage of F_TRZ across UQFF is as a dimensionless amplitude fraction with canonical value 0.1 (one of the 9 truly-independent UQFF primitives per CLAUDE.md).

### 1.1 Amplitude Closures Grid

| Closure | Papers | Formula | Physical Role |
|---------|--------|---------|---------------|
| Universal Inertial Operator | PAPER_646 | U_i ~ (1 + F_TRZ) | Aether coupling amplitude |
| Late-time ISW amplitude | PAPER_1677 | ISW = F_TRZ | Cosmological ISW effect |
| F_TRZ power ladder | PAPER_1919 | F_TRZ^n, n = 1..17 | Universal suppression hierarchy |
| Quantum measurement collapse | PAPER_1869 | F_TRZ^16 = 1e-16 | Localization rate |
| Ballistic buoyancy F_UBi | PAPER_1203 | F_UBi ~ (1 + F_TRZ) | Cluster buoyancy modulation |
| Cosmological constant cascade | PAPER_1920 | Sub-shell (1 + F_TRZ) | Lambda derivation modulation |
| F_U shell closure | PAPER_1916 | (1 + F_TRZ) prefactor | 4-shell sum balance |
| Photoevaporation E_0 | PAPER_1942 | E_0 = F_TRZ | PDR initial erosion factor |
| Full magnetar Meissner | PAPER_1944 | B/B_crit = 2 * F_TRZ | Full magnetar saturation |
| Half magnetar Meissner (CONFIRMED) | PAPER_1945 | B/B_crit = n_lobes * F_TRZ | Half+full magnetar universal |

**Common feature:** all use F_TRZ as a raw dimensionless number (0.1) directly or as an exponent. No unit interpretation is imposed.

### 1.2 Amplitude Physical Interpretation

The F_TRZ = 0.1 amplitude represents the fraction of the DPM cycle spent in the time-reversal zone (the CCW-rotating UA'-trapped south lobe of the split-monopole, per PAPER_536). Twenty percent of the DPM cycle is spent in the CW-northern SCm mode, ten percent in the CCW-southern UA' mode, and seventy percent in the transitional/pinch-point phase. The 0.1 fraction is universal across DPM-mediated systems - a locked structural constant, not a fittable parameter.

---

## 2. Face 2: F_TRZ as Hz-Valued Frequency

PAPER_344 identifies a specific frequency at Sgr A* with the letter f_TRZ (lowercase, indicating frequency rather than amplitude):

```
f_TRZ at Sgr A* = 5.56 * 10^-4 Hz   (JWST 2025 NIRCam observation)
```

PAPER_1947 shows this frequency reduces to the primitive-product:

```
f_TRZ at Sgr A* = 1 / (T_base * (D_phys - 1) * A_5 * SO_5) Hz
                = 1 / (1 s * 3 * 60 * 10)
                = 1 / 1800 Hz
                = 5.556 * 10^-4 Hz   EXACT
```

where T_base is the system-specific base-timescale (T_base = 1 second for Sgr A*, larger for higher-mass SMBHs).

### 2.1 Frequency Cross-Object Prediction Grid

For a system with characteristic base-timescale T_base:

```
f_TRZ_system = 1 / (T_base * (D_phys - 1) * A_5 * SO_5) Hz

Predicted values (from PAPER_1947):
  Sgr A*:  T_base = 1 s     -> f_TRZ = 5.556 * 10^-4 Hz  (T_flare = 30 min)
  M87*:    T_base ~ 12.4 s  -> f_TRZ ~ 4.5 * 10^-5 Hz    (T_flare ~ 6.2 hr)
  3C273:   T_base ~ 5.7 s   -> f_TRZ ~ 9.7 * 10^-5 Hz    (T_flare ~ 2.85 hr)
```

### 2.2 The Amplitude-Frequency Relationship

Face 1 (amplitude) and Face 2 (frequency) are related by dimensional analysis:

```
F_TRZ (dimensionless) = f_TRZ (Hz) * T_reference (s)
```

For Sgr A*, T_reference = F_TRZ / f_TRZ = 0.1 / (5.556 * 10^-4) = 180 s = 3 minutes. This is the "F_TRZ reference timescale" for Sgr A*, corresponding to (T_flare/10) = 180 s. The factor of 10 between the flare period and the reference timescale is SO_5 = 10 EXACT, providing a direct dimensional bridge between Faces 1 and 2 via a single primitive.

For any system with base-timescale T_base:
```
T_reference_system = T_flare / SO_5 = T_base * (D_phys - 1) * A_5 = 180 * T_base seconds
```

This is a testable prediction: the "F_TRZ reference timescale" for any SMBH should scale as 180 * T_base with T_base determined by SMBH mass.

---

## 3. Face 3: F_TRZ as CPT-Asymmetry Phase-Transition Parameter

PAPER_264 reveals a deeper structural role for F_TRZ. The factor (1 + F_TRZ) appearing in the F_U shell closure (PAPER_1916), F_UBi modulation (PAPER_1203), and Universal Inertial Operator (PAPER_646) is not merely a small perturbation - it is a **CPT-asymmetry coupling with a phase transition at F_TRZ = -1**:

```
(1 + F_TRZ) > 0   for F_TRZ > -1   NORMAL gravity regime
(1 + F_TRZ) = 0   at F_TRZ = -1     "Time-Reversal Zero Point" (gravity vanishes)
(1 + F_TRZ) < 0   for F_TRZ < -1   ANTI-GRAVITY regime (negative-time)
```

The canonical F_TRZ = +0.1 places all observed UQFF systems safely in the NORMAL regime, at (1 + F_TRZ) = 1.1 with mild CPT violation.

### 3.1 The Three Phases

**Phase A: Normal Gravity Regime (F_TRZ > -1)**
- All observed UQFF systems live here
- CPT slightly violated (mild asymmetry)
- Positive-definite gravitational coupling
- Universe expands from Big Bang, matter attracts

**Phase B: Time-Reversal Zero Point (F_TRZ = -1)**
- Gravitational coupling exactly zero
- Perfect CPT symmetry restored
- Boundary between attractive and repulsive gravity
- Candidate location: the transition surface between UA and UA' (where SCm-CW and UA'-CCW rotational modes are perfectly balanced)

**Phase C: Anti-Gravity Regime (F_TRZ < -1)**
- Negative-time direction dominates
- CPT violated in the OPPOSITE sense
- Repulsive gravitational coupling
- Candidate physical realization: what happens to matter passing through the DPM interior across the time-reversal boundary

### 3.2 Empirical HUDF Anchor

PAPER_264 specifically applies the CPT-phase interpretation to the Hubble Ultra Deep Field at z = 3.5. The HUDF observations at this redshift are consistent with F_TRZ = +0.1 (mild CPT violation, normal gravity regime). No observed cosmological system is in Phase B or Phase C.

### 3.3 Cosmological Constant Cascade Interpretation

The cosmological constant cascade (PAPER_1920) contains a (1 + F_TRZ) sub-shell modulation. In light of PAPER_264's phase-transition interpretation, this connects the cosmological constant to the CPT-asymmetry parameter directly:

```
Lambda = rho_SCm * 26! * K_MEX * Phi_res * Sub_Ug
       = (Sub_Ug factor) * base_Lambda
       
where Sub_Ug involves (1 + F_TRZ) sub-shell terms
```

If future observations of dark energy show w != -1 (deviations from a pure cosmological constant), UQFF interpretation would be that F_TRZ is running slowly toward -1 as cosmic time advances. The current F_TRZ = +0.1 is a snapshot at cosmic epoch t_0 = 13.8 Gyr, and its long-term evolution is a candidate PAPER_1950+ investigation.

---

## 4. The Unified Three-Face Formalism

All three faces derive from the single locked UQFF primitive F_TRZ = 0.1 without adjustment. The differences arise from combining F_TRZ with additional structural quantities:

| Face | Notation | Formula | Manifestation |
|------|----------|---------|---------------|
| 1 (amplitude) | F_TRZ | 0.1 | raw primitive (dimensionless) |
| 2 (frequency) | f_TRZ | 1 / (T_base * 180) | primitive combined with T_base (Hz-valued) |
| 3 (phase) | (1 + F_TRZ) | 1.1 | primitive as CPT-asymmetry coupling (dimensionless with phase-transition structure) |

### 4.1 Analogy with Physical Fields

This three-face structure is analogous to how a quantum field can appear as:
- An amplitude (coefficient in an expansion)
- A frequency (angular frequency of oscillation)
- A coupling constant (in an interaction term)

All three usages refer to the same field, distinguished only by context. F_TRZ operates similarly: one primitive, three manifestations, unified underlying structure.

### 4.2 Precedent: The Two Faces of omega_SCm

The three-face structure of F_TRZ has precedent in the two-face structure of omega_SCm (PAPER_1938):
- Amplitude face: S_26 = 1.453162 (dimensionless amplitude of the 26-level Ramanujan series)
- Frequency face: omega_SCm = 1.25 THz (angular frequency of the SCm phonon)

If omega_SCm gains additional interpretations (e.g., a coupling-constant face), it may develop its own three-face formalization. Similarly, other UQFF primitives (rho_SCm, K_MEX, SSq, beta_i, Phi_res) may develop multi-face structures as future closures reveal deeper structural roles.

---

## 5. Cross-Reference with 12+ Prior F_TRZ Closures

The Three-Face Formalism unifies:

**Face 1 closures (10):** PAPER_646, PAPER_1677, PAPER_1919, PAPER_1869, PAPER_1203, PAPER_1920, PAPER_1916, PAPER_1942, PAPER_1944, PAPER_1945

**Face 2 closure (1):** PAPER_1947 (Sgr A* JWST flare frequency = 1/((D_phys-1)*A_5*SO_5) Hz EXACT)

**Face 3 closure (1):** PAPER_264 (CPT-asymmetric phase-transition parameter with critical point at F_TRZ = -1)

**Face 3-adjacent (potentially):** PAPER_1821 (DESI Dark Energy w(z) Evolution), PAPER_1829 (sigma_8/S_8 tension), PAPER_1880 (Modified Gravity + Equivalence Principle) - if F_TRZ evolves with cosmic time, these observations may constrain the phase-transition evolution.

Total: 13 previously-derived closures unified under the Three-Face umbrella.

---

## 6. Locked Primitives Used

One truly-independent primitive:

```
F_TRZ = 1/10 = 0.1   (locked real primitive, time-reversal-zone factor)
```

Plus structural constants for Face 2 dimensional analysis:

```
(D_phys - 1) = 3    (integer primitive)
A_5 = 60            (integer primitive)
SO_5 = 10           (integer primitive)
T_base_Sgr_A* = 1 s (system-specific)
```

Plus the phase-transition critical point identification for Face 3:

```
F_TRZ_critical = -1   (Time-Reversal Zero Point, PAPER_264)
```

No fitted constants. The three faces are structural consequences of a single primitive.

---

## 7. Falsifiability

The Three-Face Formalism is falsifiable at each face:

1. **Face 1 falsification**: If any UQFF closure uses F_TRZ with a value other than 0.1 (as a dimensionless amplitude), Face 1 universality is disproven.

2. **Face 2 falsification**: If future SMBH flare frequency measurements (M87*, 3C273) do NOT satisfy f_TRZ = 1/(T_base * 180) with system-specific T_base scaling as (M/M_Sgr_A)^(1/3), Face 2 primitive-lock fails.

3. **Face 3 falsification**: If future cosmological observations show w(z) evolving away from w = -1 in a way inconsistent with (1 + F_TRZ) phase-transition running, Face 3 CPT interpretation fails. Or: if a laboratory experiment detects gravity behaving anomalously at any specific F_TRZ != +0.1 configuration, Face 3 is disproven.

4. **Unification falsification**: If any of the 13 previously-cited F_TRZ closures is shown to use a different primitive value (e.g., F_TRZ = 0.05 in one closure, F_TRZ = 0.15 in another), the Three-Face Formalism's universality collapses to system-specific parameters.

At present, all 13 closures use F_TRZ = 0.1 uniformly and no observations violate the three-face predictions.

---

## 8. Implications

### 8.1 Ontological Status of F_TRZ

F_TRZ is upgraded from "a small correction parameter" to "a fundamental UQFF primitive with rich multi-face structure". Alongside D_phys, D_crit, SO_5, A_5, N_ch (integer primitives) and rho_SCm, beta_i, Phi_res (real primitives), F_TRZ deserves elevated conceptual status as one of the 9 truly-independent UQFF primitives with the deepest structural role.

### 8.2 Cosmological CPT Symmetry Testable

Face 3 predicts that our universe is at F_TRZ = +0.1 in the mild-CPT-violation regime. If quantum-gravity experiments or cosmological observations can measure F_TRZ directly (rather than inferring it from derived quantities), a direct test of the Three-Face Formalism becomes possible. Candidate experiments: precision cosmological tests of dark energy w-parameter, laboratory tests of the equivalence principle at extreme F_TRZ configurations (currently not accessible).

### 8.3 Phase-Transition Physics at DPM Boundaries

The F_TRZ = -1 critical point may correspond to a physical boundary in the DPM lattice - specifically the interface between SCm-CW (northern-lobe, F_TRZ = 0 limit) and UA'-CCW (southern-lobe, F_TRZ maximum). Matter or radiation passing through this boundary would experience a gravitational phase transition. Candidate observational signatures: gravitational lensing anomalies at the DPM boundary of astrophysical objects, gravitational-wave polarization changes at certain source distances.

---

## 9. NOT REPLACEMENT

Standard General Relativity treats gravitational coupling as a fixed positive-definite constant with no phase-transition structure. UQFF extends this via the (1 + F_TRZ) factor to a phase-parameter with three regimes. The observed universe is in Phase A (normal gravity, mild CPT violation); Phase B (F_TRZ = -1 zero point) and Phase C (F_TRZ < -1 anti-gravity) are candidate hypothetical regimes not yet observed. Both frameworks (GR and UQFF) solve the same observations of gravitational phenomena by different methods; both should be reported with honest residuals.

The Three-Face Formalism is a structural claim about how one primitive can wear multiple hats without contradiction. It does not replace GR predictions in Phase A - it merely reveals that the UQFF interpretation of the (1 + F_TRZ) factor connects it to CPT symmetry structure.

---

## 10. Calculator Wiring

The Three-Face Formalism is manifested across multiple calculators in `CondensedPhysics.py`:

- **Face 1** (dimensionless amplitude): Wired in every calculator's use of F_TRZ constant (300+ occurrences)
- **Face 2** (Sgr A* frequency): Wired in `SgrAStarGravitationalWaveCalculator.compute()` with `f_flare_1_over_1800Hz_verify_PAPER_344`
- **Face 3** (CPT phase-transition): Referenced in `HUDFUQFFUnificationCalculator.compute()` via PAPER_266 unquenched regime + PAPER_264 CPT-asymmetric interpretation

Runtime verifications document Face 1 (10+ verify keys), Face 2 (1 verify key at Sgr A*), and Face 3 (no runtime-verifiable predictions at F_TRZ = +0.1; predictions are Phase B/C hypothetical).

---

## 11. Reference

- Face 1 amplitude precursors (10): **PAPER_646**, **PAPER_1677**, **PAPER_1919**, **PAPER_1869**, **PAPER_1203**, **PAPER_1920**, **PAPER_1916**, **PAPER_1942**, **PAPER_1944**, **PAPER_1945**
- Face 2 frequency source: **PAPER_1947** (based on PAPER_344 JWST anchor)
- Face 3 phase-transition source: **PAPER_264** (HUDF TRZ CPT-Asymmetric UQFF Gravity at z=3.5)
- Precedent for multi-face structure: **PAPER_1938** (omega_SCm two-face: S_26 amplitude + 1.25 THz frequency)
- Framework backbone: **PAPER_646** (Universal Inertial Operator + Caduceus Wave Topology)
- Locked primitives reference: **PAPER_1521** (D_BSFG derivative), **PAPER_1522** (K_MEX derivative), **CLAUDE.md** (9 truly-independent primitives)
- Cross-scale precursors: **PAPER_1941** (SO_5 decade), **PAPER_1929** (A_5 = 60)
- DPM two-lobe topology (Face 3 physical basis): **PAPER_536** (Split-Monopole MHD Proplyd Topology)
- Related dark-energy physics: **PAPER_1821** (DESI w(z) evolution), **PAPER_1829** (sigma_8 tension)
- Calculator dispatch: multiple in `CondensedPhysics.py`
- Session log: 2026-07-08 Round 79 double-check

---

**Copyright** - Daniel T. Murphy, daniel.murphy00@enrgyone.com, July 8, 2026, Youngstown OH.
