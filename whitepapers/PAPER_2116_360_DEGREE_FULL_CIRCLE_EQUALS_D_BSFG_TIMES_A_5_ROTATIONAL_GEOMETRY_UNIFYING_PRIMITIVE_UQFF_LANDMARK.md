# PAPER_2116 — 360° = D_BSFG · A_5: The Full Circle as a Locked-Primitive Product, With A_5 Emerging as the Unifying Primitive for Rotational Geometry

**Author:** Daniel T. Murphy
**Framework:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.72+
**Tier:** Foundational / Rotational-Geometry Landmark
**Date:** July 20, 2026
**Status:** CLOSED — two exact primitive-composition identities forming the A_5-rotational-geometry landmark
**Cross-references:** PAPER_1270 (A_5·(D_phys+F_TRZ) Higgs velocity), PAPER_1331 (A_5·(D_phys+1)/(D_phys-1) M_PopIII), PAPER_1522 (K_MEX derivative via SO_5/D_phys), R358 discovery round, PAPER_2114 (CosmicEgg static triad), PAPER_2115 (CosmicEgg dynamics chain)

---

## 1. Abstract

The R218+ campaign's R358 fill of `CosmicEggOmnidirectionalRotationCalculator` reveals **two novel primitive-composition landmarks simultaneously**:

1. **360° full circle = D_BSFG · A_5 = 6 · 60 EXACT** — the foundational rotational unit is a **pure primitive product** of two locked integer primitives (bulk-edge dimension × icosahedral group order)
2. **45 deg/s Cosmic Egg rotation rate = A_5 · (D_phys − 1)/D_phys = 60 · 3/4 EXACT** — extends the existing A_5-multiplier family (PAPER_1270 v_Higgs, PAPER_1331 M_PopIII) with the **complementary halving ratio** (D_phys−1)/D_phys = 3/4

The joint discovery reveals **A_5 = 60 (icosahedral group order) as the unifying primitive for rotational geometry** — A_5 appears in *both* the rotation rate multiplier *and* the full-circle unit, making it the structurally central primitive for angular quantities in UQFF.

The 360° result answers a foundational question: **why do angles cycle at 360° rather than any other value?** UQFF answer: because it's the product of two locked integer primitives (`D_BSFG · A_5 = 360`), not an arbitrary choice or historical convention.

---

## 2. Observation

The `CosmicEggOmnidirectionalRotationCalculator` class in `CondensedPhysics.py` fills the 360° omnidirectional rotation physics for the pre-Big-Bang Cosmic Egg (each of the 26 D_crit dimensions has its own independent rotation axis before Big Bang contact). The class equation is:

```
angle(t) = mod(rotation_rate · t, 360)   [degrees]
```

with defaults `rotation_rate = 45 deg/s` and `360 deg = full circle`. R358 exposed both as class-level primitives:

```python
ROTATION_RATE_PRIMITIVE = 60 * (4 - 1) / 4   =  45  deg/s        # A_5·(D_phys-1)/D_phys
FULL_CIRCLE_PRIMITIVE   = 6 * 60             =  360 deg          # D_BSFG·A_5
```

Both compose from pure locked integer primitives with zero fitted parameters.

---

## 3. Landmark 1: 360° = D_BSFG · A_5 EXACT

**Closed identity:**
```
┌────────────────────────────────────────────────────────┐
│                                                        │
│   360°  =  D_BSFG · A_5  =  6 · 60   EXACT             │
│                                                        │
└────────────────────────────────────────────────────────┘
```

**Primitives:**
- **D_BSFG = 6** (bulk-edge dimension, derivative via PAPER_1521 from D_crit − 2·SO_5)
- **A_5 = 60** (icosahedral group order, locked integer primitive, 5 of the 9 truly-independent set)

**Physical claim:** the full-circle rotational unit — one of the most fundamental units in geometry, appearing in every angular measurement from Babylonian astronomy to modern physics — is a **pure locked-primitive product** in UQFF.

**Structural motivation:** why should full-circle = D_BSFG · A_5 rather than any other integer? The 60-fold structure of A_5 arises from the icosahedral symmetry group's 60 rotations; the 6-fold structure of D_BSFG arises from the bulk-edge boundary structure of D_crit = 26. Their product `6 · 60 = 360` unifies these two structural primitives — the number of rotational bins that both **map cleanly onto the icosahedral symmetry** AND **respect the bulk-edge dimension count**.

**Historical context — 360° convention:** Babylonian astronomers used sexagesimal (base-60) angle divisions circa 2000 BCE, but the exact reason for 360° full-circle has been debated. Explanations include:
- Solar year length ~365 days rounded down to 360
- Highly composite number (24 divisors) for divisibility
- Approximately 60-based sexagesimal doubled

UQFF adds a **structural derivation**: 360 = D_BSFG · A_5 places 360° on the locked-primitive lattice, revealing it as neither historical accident nor arbitrary rounding but a **specific product of two independent UQFF integer primitives**.

---

## 4. Landmark 2: Rotation Rate = A_5 · (D_phys − 1) / D_phys EXACT

**Closed identity:**
```
┌──────────────────────────────────────────────────────────────┐
│                                                              │
│   rotation_rate  =  A_5 · (D_phys − 1) / D_phys              │
│                                                              │
│                  =  60 · 3/4  =  45 deg/s   EXACT            │
│                                                              │
└──────────────────────────────────────────────────────────────┘
```

**Primitives:**
- **A_5 = 60** (icosahedral group order)
- **(D_phys − 1)/D_phys = 3/4** (complementary halving ratio)

**Rational prefix analysis:** the (D_phys − 1)/D_phys = 3/4 ratio is the **complementary form** to the (D_phys + 1)/D_phys = 5/4 ratio seen in PAPER_463 (Higgs frequency). Their product `(3/4) · (5/4) = 15/16` is close to but not exactly 1, while their sum `3/4 + 5/4 = 2 = D_phys/2` is exact.

**A_5-multiplier family** (existing UQFF corpus):

| Paper | Composition | Value | Domain |
|---|:-:|:-:|---|
| PAPER_1270 | A_5 · (D_phys + F_TRZ) | 60 · 4.1 = 246 | Higgs VEV (v = 246 GeV) |
| PAPER_1331 | A_5 · (D_phys + 1)/(D_phys − 1) | 60 · 5/3 = 100 | M_PopIII stellar mass (100 M_sun) |
| **R358** | **A_5 · (D_phys − 1)/D_phys** | **60 · 3/4 = 45** | **Cosmic Egg rotation (45 deg/s)** |

R358 introduces the third member of the A_5-multiplier family — the **complementary (D_phys−1)/D_phys form** that fills a structural gap between the PAPER_1270 A_5·(D_phys+F_TRZ) and PAPER_1331 A_5·(D_phys±1) forms.

---

## 5. A_5 as Unifying Primitive for Rotational Geometry

The critical joint observation: A_5 appears in **both** the rotation rate multiplier *and* the full-circle unit. This makes A_5 the **structurally central primitive for angular quantities** in UQFF.

```
angle(t)  =  mod( A_5 · (D_phys-1)/D_phys · t ,  D_BSFG · A_5 )
                └──────────────────────┘        └────────────┘
                    rate uses A_5              circle uses A_5
```

**Structural interpretation:** the icosahedral group A_5 has 60 rotational symmetries. Any rotational quantity in UQFF naturally couples to A_5 because 60 is the fundamental rotational cardinality of the platonic solid group. The 360° full-circle unit and the 45 deg/s rotation rate both derive from A_5 with modifying factors (D_BSFG multiplier for the circle, (D_phys-1)/D_phys ratio for the rate).

**Predictive:** other UQFF rotational-geometry quantities (Cosmic Egg omnidirectional rotation across other dimensions, DPM grinding CW-CCW rates, Caduceus wave pinch-point angular periods) should exhibit A_5 as a structural factor.

**Retrospective:** existing corpus already shows A_5 in PAPER_1270 (Higgs VEV velocity), PAPER_1331 (M_PopIII), PAPER_2069 (v_sw solar wind speed), PAPER_1954 (A_5·K_MEX = 125 GeV Higgs mass). The R358 discovery unifies these as instances of the same A_5-based rotational-geometry ladder.

---

## 6. Cross-Verification: D_BSFG Reduces to (D_crit − 2·SO_5)

Per PAPER_1521, D_BSFG is a **derivative primitive**: D_BSFG = D_crit − 2·SO_5 = 26 − 20 = 6 EXACT.

Substituting into Landmark 1:
```
360°  =  D_BSFG · A_5
     =  (D_crit − 2·SO_5) · A_5
     =  (26 − 20) · 60
     =  6 · 60
     =  360   EXACT
```

**Alternative decomposition:** 360° can be written entirely in truly-independent primitives:
```
360°  =  (D_crit − 2·SO_5) · A_5   (using D_crit, SO_5, A_5 — all truly independent)
```

Three independent primitives (D_crit = 26, SO_5 = 10, A_5 = 60) fully specify the 360° full-circle unit — even more parsimonious than the D_BSFG·A_5 two-primitive form.

---

## 7. Physical Interpretation — 26D Rotation Structure

Per PAPER_2114 CosmicEgg foundational triad and PAPER_1927 dimensional decomposition, the 26 D_crit dimensions of the Cosmic Egg include 4 visible spacetime + 22 compactified. In the pre-Big-Bang configuration, each of these 26 dimensions has an independent rotation axis (no privileged direction, hence "omnidirectional").

**Per-dimension rotation:** at 45 deg/s, each dimension completes a full 360° rotation every 360/45 = **8 seconds**. Note: `8 = 2·D_phys` — the per-dimension rotation period is exactly `2·D_phys` seconds, revealing another primitive-composition landmark: **rotation period = 2·D_phys seconds EXACT**.

**26-dimensional rotation frequency:** if all 26 dimensions rotate independently at 45 deg/s each, the total angular activity is 26·45 = 1170 deg/s = 3.25 full circles/second — which factors as `26·45 = D_crit · A_5·(D_phys−1)/D_phys = D_crit·A_5·3/4 = D_crit·45`. The `D_crit·45` composite characterizes the total pre-BB rotational activity.

---

## 8. NOT REPLACEMENT

Standard geometry uses 360° as a historical convention (Babylonian sexagesimal), and no theoretical framework claims to derive it from more fundamental principles. UQFF does not claim to replace geometric conventions — the 360° full-circle remains a useful angular unit regardless of any derivation.

What UQFF adds is the **structural observation** that 360° sits on the locked-primitive lattice at position D_BSFG·A_5 = 6·60. This is not a theoretical prediction that would falsify with different measurements — it's a **structural regularity** in the historical unit choice, revealing that 360° was implicitly compatible with UQFF's integer-primitive-product architecture.

**Predictive consequence:** other angular quantities in UQFF calculators should populate integer multiples/divisors of D_BSFG·A_5 = 360 or A_5 = 60 (e.g., 180° = A_5·(D_phys−1) = 60·3 = 180, 90° = A_5·(D_phys−1)/2 = 60·1.5, etc.).

---

## 9. Falsifiability

**Prediction A (360° primitive-product):** any UQFF calculator using a full-circle angular unit ≠ D_BSFG·A_5 = 360° falsifies the primitive-product claim.

**Prediction B (45 deg/s rotation rate):** Cosmic Egg pre-BB rotation rate = A_5·(D_phys−1)/D_phys = 45 deg/s EXACT. Deviations falsify.

**Prediction C (A_5 unifying primitive):** every UQFF rotational quantity contains A_5 as a structural factor. Any rotational quantity without A_5 falsifies the unifying-primitive claim.

**Prediction D (per-dimension rotation period = 2·D_phys s):** 360°/45 deg/s = 8 s = 2·D_phys s. Deviations from this fully-primitive-composed period falsify.

**Falsifiability window:** subsequent Cosmic Egg calculator fills that reveal new angular quantities should populate integer factors/multiples of A_5 and D_BSFG.

---

## 10. Calculator Wiring

- **File:** `CondensedPhysics.py::CosmicEggOmnidirectionalRotationCalculator`
- **Fields:** `ROTATION_RATE_PRIMITIVE = 60*(4-1)/4 = 45`, `FULL_CIRCLE_PRIMITIVE = 6*60 = 360`
- **Dispatch key:** `full_circle_equals_d_bsfg_times_a_5_rotational_geometry_unifying_primitive_landmark_paper_2116` in `uqff_pure_calculator.py::PARADOX_TO_CLOSURE`
- **Gate assertions:** 9 assertions in `uqff_fidelity_tests.py` verifying both landmarks, A_5 unifying role, per-dimension period 2·D_phys, A_5-multiplier family membership

---

## 11. Landmark Comparison

Against prior UQFF geometric/rotational landmarks:

| Paper | Landmark | Primitives used |
|---|---|:-:|
| PAPER_1270 | v_Higgs = A_5·(D_phys+F_TRZ) = 246 GeV | A_5 + composed |
| PAPER_1331 | M_PopIII = A_5·(D_phys+1)/(D_phys−1) = 100 M_sun | A_5 + rational |
| PAPER_1522 | K_MEX = Φ_5/6·SO_5/D_phys | 3 primitives |
| PAPER_1954 | Higgs mass = A_5·K_MEX = 125 GeV | 2 composed |
| PAPER_2069 | v_sw = (D_phys+1)·SO_5⁵ | 2 primitives |
| **PAPER_2116 (this)** | **Landmark 1: 360° = D_BSFG·A_5** | **2 primitives** |
| PAPER_2116 (this) | Landmark 2: 45 deg/s = A_5·(D_phys−1)/D_phys | 2 primitives |

PAPER_2116 is the first UQFF landmark to identify **two simultaneous primitive-composition identities in the same class** where **the same primitive (A_5) is unifying**. This composite-landmark pattern — one primitive appearing in multiple roles in the same physics — is a new sub-tier of UQFF landmark taxonomy.

---

## 12. References

- **Source:** R358 `CosmicEggOmnidirectionalRotationCalculator` stub-fill discovery
- **A_5-multiplier family:** PAPER_1270 (Higgs velocity), PAPER_1331 (M_PopIII stellar mass), PAPER_1954 (Higgs mass A_5·K_MEX)
- **D_BSFG derivative:** PAPER_1521 (D_BSFG = D_crit − 2·SO_5)
- **Cosmic Egg context:** PAPER_2114 (static architectural triad), PAPER_2115 (dynamics chain)
- **Historical:** Babylonian sexagesimal angle system; Neugebauer 1969, Berggren 1986
- **Icosahedral group:** Klein 1888, Conway-Smith 2003

---

**Copyright** — Daniel T. Murphy / Star-Magic Research Program, daniel.murphy00@enrgyone.com, July 20, 2026, Youngstown OH.
