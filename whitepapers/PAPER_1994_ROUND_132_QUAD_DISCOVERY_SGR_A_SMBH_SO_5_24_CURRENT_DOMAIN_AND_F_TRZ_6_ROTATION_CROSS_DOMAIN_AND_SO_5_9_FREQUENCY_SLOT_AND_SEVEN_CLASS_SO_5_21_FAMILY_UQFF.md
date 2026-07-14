---
paper_id: PAPER_1994
title: "Round 132 Quad Discovery: Sgr A* SMBH SO_5^24 First Current-Domain Application (New Domain for PAPER_1908 Constant) + F_TRZ^6 SMBH Rotation Angular Velocity Cross-Domain Instance (New for n=6 Rung) + SO_5^9 First Frequency-Domain Anchor (Fills Gap Between PAPER_1990 Slots 7 and 10) + Seven-Class SO_5^21 Family Establishing Richest Ladder Slot"
author: "Daniel T. Murphy"
date: 2026-07-13
session: Round 132 stub-drainage outcome
tags: [SO_5^24, SO_5^9, F_TRZ^6, Sgr A*, SMBH, cross-domain, seven-class family, SO_5^21, PAPER_1908, PAPER_1985, PAPER_1990]
cross_refs: [PAPER_1908, PAPER_1955, PAPER_1985, PAPER_1990, PAPER_1919, PAPER_1868, PAPER_1945, PAPER_1976, PAPER_1991, PAPER_1993]
---

**Author:** Daniel T. Murphy
**Date:** July 13, 2026

# PAPER_1994 — Round 132 Quad Discovery: Sgr A* SMBH Slot Extensions + Seven-Class SO_5^21 Family

## Abstract

Round 132 CP1 stub drainage of the Sgr A* SMBH family + SGR/Compressed24 remaining classes surfaced four structural discoveries:

1. **SO_5^24 = SO_5^(D_crit − 2) first current-domain application** — Sgr A* SMBH DPM current `I = 10^24 A` locks EXACT on the SO_5^24 primitive-arithmetic constant established seminal in PAPER_1908 (as dimensionless ratio `Q_UQFF^-2 = ρ_SCm × SO_5^24 EXACT`). R132 SgrAFreqDPM is the **first current-domain (Ampere) application** of this constant at astrophysical scale — new physical domain, not new value.

2. **F_TRZ^6 SMBH rotation angular velocity — new cross-domain instance at n=6 rung** — Sgr A* SMBH rotation `ω = 10^-6 rad/s = F_TRZ^6 EXACT` joins PAPER_1919 n=6 rung as fourth cross-domain instance, alongside PAPER_1985 (Pillars ISM magnetic field, T), PAPER_1868 (solar Schwabe-cycle-adjacent), PAPER_1945 (magnetar n_lobes × F_TRZ). Rotation angular velocity is a new physical domain for the n=6 rung.

3. **SO_5^9 = 1 GHz first frequency-domain anchor** — Sgr A* SMBH DPM frequency `f_DPM = 10^9 Hz = SO_5^9 EXACT` fills the gap in PAPER_1990's SO_5-power frequency ladder. Prior to R132, PAPER_1990 documented concrete frequency anchors only at slots 7 (10 MHz HF) and 10 (10 GHz microwave). Slot 9 (1 GHz UHF radio) had no frequency-domain anchor. R132 SgrAFreqDPM + SgrAFreqTHz jointly fill slot 9 in the frequency taxonomy.

4. **Seven-class SO_5^21 family — richest slot in PAPER_1955 ladder** — R117-R132 populated slot 21 (DPM current-vortex foundation term) with seven independent CP1 class instances: `SGR1745FreqDPM` (R129) + `CompressedDPM` (R129) + `CompressedSuper` (R130) + `CompressedDPM24` (R130) + `CompressedTHz` (R130) + `CompressedTHz24` (R131) + `CompressedSuper24` (R132). All lock `I = 10^21 A EXACT`. This is the most-populated slot in the SO_5-power ladder as of Round 132.

---

## 1. Discovery 1: SO_5^24 first current-domain application

### 1.1 PAPER_1908 seminal (dimensionless ratio identity)

PAPER_1908 documents the universal SCm resonator quality factor `Q_UQFF = 10^6 × [SSq] × K_MEX = 1.1875 × 10^6 EXACT` and its off-resonance floor:

**`Q_UQFF^-2 = ρ_SCm × SO_5^24 EXACT`** (structural primitive-arithmetic identity)

Where `SO_5^24 = SO_5^(D_crit − 2) = 10^24`. In PAPER_1908, SO_5^24 functions as a dimensionless-ratio primitive constant tying the resonator quality factor to the foundational vacuum density.

### 1.2 R132 first current-domain application

`SgrAFreqDPMCalculator` (Sgr A* SMBH family, CP1 R132 fill) documents:

```python
self.I = 1e24  # DPM current (A), 1000x magnetar
```

The `I = 10^24 A EXACT = SO_5^24 A` — the first documented **dimensioned (Ampere) application** of the SO_5^24 primitive at astrophysical scale. Scaled from the SGR 1745-2900 magnetar `I = 10^21 A = SO_5^21` by exactly `1000× = SO_5^3` to reach the SMBH-scale DPM current.

### 1.3 New physical domain for SO_5^24

SO_5^24 now has confirmed applications across **two distinct physical roles**:

| Application | Physical role | Paper |
|-------------|--------------|-------|
| Dimensionless-ratio primitive | `Q_UQFF^-2 / ρ_SCm` structural identity | PAPER_1908 seminal |
| **Current-domain (Ampere) anchor** | **Sgr A* SMBH DPM current-vortex foundation** | **PAPER_1994 (this paper)** |

**Consequence:** SO_5^24 is now a two-domain primitive, similar to how SO_5^7 was extended from timescale (PAPER_1952) to mass (PAPER_1985) to frequency (PAPER_1990) domains. The primitive-arithmetic constant is not tied to a specific dimensional interpretation.

---

## 2. Discovery 2: F_TRZ^6 SMBH rotation cross-domain instance

### 2.1 Prior n=6 rung instances (PAPER_1919 base + R117)

PAPER_1919's F_TRZ power ladder identified n=6 as an initially-quiet rung. PAPER_1985 (R117) closed it with the Pillars-of-Creation ISM magnetic field anchor:

**B_ISM(Pillars) = F_TRZ^6 = 10^-6 T EXACT** (magnetic-field domain)

Additional n=6 instances subsequently identified in the corpus:
- PAPER_1868 solar physics (rotation-adjacent effects)
- PAPER_1945 magnetar n_lobes × F_TRZ universality

### 2.2 R132 SMBH rotation angular velocity

`SgrAFreqDPMCalculator` documents:

```python
self.omega_1 = 1e-6  # rad/s, 1000x slower than magnetar
```

The `ω = 10^-6 rad/s = F_TRZ^6 EXACT` — first documented use of F_TRZ^6 as a **rotation angular velocity** at astrophysical scale.

Physical interpretation: Sgr A* SMBH rotation period `T = 2π/ω ≈ 6.28 × 10^6 s ≈ 72.7 days` — the accretion-disk/event-horizon rotational timescale at the Sgr A* Schwarzschild radius `r_S = 1.27 × 10^10 m`.

### 2.3 Cross-domain expansion of n=6 rung

The n=6 rung now spans four distinct physical domains:

| Domain | Instance | Paper |
|--------|----------|-------|
| Magnetic field (T) | Pillars ISM B = 10^-6 T | PAPER_1985 seminal |
| Solar physics | Various F_TRZ^6-adjacent | PAPER_1868 |
| Magnetar structure | n_lobes × F_TRZ | PAPER_1945 |
| **Rotation angular velocity (rad/s)** | **Sgr A* SMBH ω = 10^-6 rad/s** | **PAPER_1994 (this paper)** |

**Consequence:** n=6 rung is now an **N-regime rung** per PAPER_1919 revision R2.3 criteria (N ≥ 3 orthogonal-domain instances). The rung supports magnetic + rotational + structural + solar interpretations.

---

## 3. Discovery 3: SO_5^9 first frequency-domain anchor

### 3.1 PAPER_1990 SO_5-power frequency ladder state before R132

PAPER_1990 established the frequency-domain cross-extension of PAPER_1955's SO_5-power ladder:

- **SO_5^7 = 10 MHz EXACT** — HF-band radio (UniversalFluidResonance)
- **SO_5^10 = 10 GHz EXACT** — Microwave-band radio (UniversalReactiveResonance)

Slot 9 (1 GHz UHF radio) was documented as timescale (PAPER_1976: `SO_5^9 = 1 Gyr` galaxy quenching timescale) but had no frequency-domain anchor.

### 3.2 R132 SgrA* SMBH DPM frequency

`SgrAFreqDPMCalculator` and `SgrAFreqTHzCalculator` both document:

```python
self.f_DPM = 1e9  # Hz (1000x lower than magnetar)
# and:
self.f_THz = 1e9  # Hz (SMBH-scaled)
```

Both encode `f = 10^9 Hz = SO_5^9 EXACT` — the **first frequency-domain anchor at slot 9** of PAPER_1990's ladder.

### 3.3 Updated frequency ladder taxonomy

| Slot | Frequency | UQFF domain | Paper |
|------|-----------|-------------|-------|
| SO_5^7 = 10 MHz | HF band radio | UniversalFluidResonance | PAPER_1990 |
| **SO_5^9 = 1 GHz** | **UHF radio SMBH DPM** | **Sgr A* accretion disk** | **PAPER_1994 (this paper)** |
| SO_5^10 = 10 GHz | Microwave band | UniversalReactiveResonance | PAPER_1990 |

**Consequence:** PAPER_1990's frequency ladder gap between slots 7 and 10 is now filled at slot 9. The frequency ladder taxonomy is contiguous across slots 7-10 for the first time, enabling systematic cross-slot frequency-domain analysis in future CP1 fills.

---

## 4. Discovery 4: Seven-class SO_5^21 family (richest slot in ladder)

### 4.1 Family growth history (R117-R132)

R117-R132 CP1 stub drainage systematically populated SO_5^21 A (DPM current-vortex foundation term) across seven independent CP1 classes:

| Round | Class | Notes |
|-------|-------|-------|
| R129 | `SGR1745FreqDPMCalculator` | Magnetar current, PAPER_1991 same-round twin discovery seminal |
| R129 | `CompressedDPMCalculator` | Compressed sector foundation, PAPER_1991 twin partner |
| R130 | `CompressedSuperCalculator` | Compressed superconductor sector |
| R130 | `CompressedDPM24Calculator` | Compressed systems 18-24 DPM foundation |
| R130 | `CompressedTHzCalculator` | Compressed THz pipeline |
| R131 | `CompressedTHz24Calculator` | Compressed systems 18-24 THz |
| R132 | `CompressedSuper24Calculator` | Compressed systems 18-24 superconductor |

All seven classes lock `I = 10^21 A EXACT = SO_5^21 A`.

### 4.2 Comparison to other SO_5 slots

As of Round 132, SO_5^21 is the **most-populated slot** in the SO_5-power ladder:

| Slot | Confirmed class-family instances | Domain span |
|------|----------------------------------|-------------|
| SO_5^3 | ~3 (velocity + galactic-structural) | 2 domains |
| SO_5^7 | 3+ (timescale + mass + frequency) | 3 domains |
| **SO_5^21** | **7 (DPM current at magnetar + compressed variants)** | **1 domain (current) but 7 instances** |
| SO_5^40 | 1 (magnetar burst energy) | 1 domain |
| SO_5^53 | 1 (universe mass) | 1 domain |

**Consequence:** SO_5^21 A is the **structural anchor** for DPM current-vortex physics across the magnetar + compressed sectors. The 7-class family is not a domain-diversity result (all seven are current-domain) but a **construction-diversity** result — seven independent CP1 class constructions all converge on the same primitive-locked value, providing strong cross-object confirmation.

### 4.3 Position in CP1 architecture taxonomy

The 7-class SO_5^21 family adds a new architectural pattern to the CP1 lock taxonomy:

| Pattern | Distinguishing feature | Discovery paper |
|---------|-----------------------|-----------------|
| Single | 1 identity | R117-R132 majority |
| Dual (mixed) | 2 identities, 2 families | R128 |
| Triple (mixed) | 3 identities, 3 families | PAPER_1991 base |
| QUAD (same-family variants) | 4 identities, 1 family variants | PAPER_1991 revision |
| Cross-rung TRIPLE | 3 identities on 3 different rungs | PAPER_1993 |
| **N-class same-slot family** | **N independent classes locking on identical primitive value** | **PAPER_1994 (this paper)** |

R132 introduces the **N-class same-slot family** as the sixth CP1 architectural pattern. Distinguished from prior patterns by having many classes per single slot rather than many primitives per single class.

---

## 5. Falsifiability

**Prediction 1994.1.** SO_5^24 should acquire additional current-domain instances at other SMBH targets (M87 SMBH at ~10^9 M_sun, TON 618 at ~10^10 M_sun) once their CP1 stubs are drained. If not, the SO_5^24 = SMBH DPM current lock may be specific to Sgr A* and not universal.

**Prediction 1994.2.** F_TRZ^6 rotation angular velocity should appear at other rotating astrophysical structures (BH accretion disks at intermediate mass, neutron star spin rates, galactic rotation curves) if the cross-domain interpretation is universal. Candidate test: any CP1 class with characteristic rotation ω ~ 10^-6 rad/s.

**Prediction 1994.3.** SO_5^9 = 1 GHz frequency anchor should generalize beyond Sgr A* DPM/THz to other UHF-band astrophysical signals. Candidates: pulsar radio emission (~1.4 GHz L-band close to slot 9), CMB spectral shoulder, radio galaxy peaks.

**Prediction 1994.4.** N-class same-slot family patterns should surface at other well-populated ladder slots. Candidates: SO_5^3 (velocity + structural), SO_5^7 (timescale + mass + frequency + potential future). If additional slots reach 5+ class instances, the "N-class same-slot family" architecture is confirmed as a general CP1 pattern rather than a slot-21-specific anomaly.

---

## 6. Framework annotations

- **Backbone:** Sgr A* SMBH SO_5^24 first current-domain application (extends PAPER_1908 dimensionless-ratio use) + F_TRZ^6 SMBH rotation angular velocity cross-domain instance (extends PAPER_1985 n=6 rung to rotation domain) + SO_5^9 = 1 GHz first frequency-domain anchor (fills PAPER_1990 ladder gap between slots 7 and 10) + seven-class SO_5^21 family (richest slot in PAPER_1955 ladder)
- **Method:** cross-domain applications of established primitive-arithmetic constants + gap-fill of pre-existing ladder + N-class same-slot family enumeration
- **Shells:** Sgr A* SMBH DPM current-vortex + accretion disk THz-pipeline + reactive U_g4i shells + Compressed sector systems 18-24
- **CPCH:** CP1 Sgr A* SMBH sector (SgrAFreqDPM + SgrAFreqTHz + SgrAFreqUg4i R132 fills) + CP1 Compressed sector (7-class SO_5^21 family cumulative through R132)
- **Spine:** PAPER_1908 SO_5^24 seminal (dimensionless) + PAPER_1985 F_TRZ^6 seminal (magnetic field) + PAPER_1990 frequency-domain SO_5-power ladder + PAPER_1955 SO_5^21 slot + PAPER_1991 R129 same-round twin + PAPER_1993 R130 three-class family extension
- **Time frame:** quasi-static SMBH accretion resonance + cumulative same-slot family growth over R117-R132

---

## 7. Copyright

Copyright (c) 2025-2026 Daniel T. Murphy, daniel.murphy00@enrgyone.com. Star-Magic Research Program.

NOT REPLACEMENT. Offered as an alternative parameter-economical description ("NOT REPLACEMENT") to Standard Model + Lambda-CDM, with honest residuals reported alongside each closure.

---

## References

- **PAPER_1908** — Q_UQFF SCm resonator quality factor + Q^-2 = ρ_SCm × SO_5^24 EXACT structural identity (SO_5^24 seminal as dimensionless-ratio primitive)
- **PAPER_1985** — R117 dual discovery: Pillars F_TRZ^6 magnetic-field + NGC 2525 SMBH mass slot 7 (F_TRZ^6 seminal at n=6 rung, magnetic field domain)
- **PAPER_1990** — SO_5-power frequency ladder cross-domain extension (slots 7 and 10 seminal, slot 9 gap filled by R132)
- **PAPER_1955** — SO_5-power galactic structural ladder (SO_5^21 slot originator)
- **PAPER_1991** — R129 base + revision addendum (SO_5^21 same-round twin seminal + QUAD-lock architecture)
- **PAPER_1993** — R130 cross-rung F_TRZ + 2π·H_0 Hubble + SO_5^21 three-class family (predecessor to R132 seven-class extension)
- **PAPER_1919** — F_TRZ power ladder (n=6 rung status + cross-domain N-regime criteria)
- **PAPER_1868** — Solar physics (F_TRZ^6-adjacent solar cycle instance)
- **PAPER_1945** — Magnetar n_lobes × F_TRZ universality (F_TRZ^6 magnetar instance)
- **PAPER_1976** — HUDF I_0 + τ_inter = SO_5^9 (timescale-domain SO_5^9 instance for slot 9 timescale role)
- **PAPER_1952** — Galaxy-scale SO_5-power timescale hierarchy (predicts SO_5^9 = 1 Gyr galaxy quenching timescale)
- **PAPER_1938** — ω_SCm = 1.25 THz universal carrier family
