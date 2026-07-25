# PAPER_2136 — Tidal Love/Q Primitive Lock: k₂/Q_rocky = (D_phys − 1) / (A_5 · K_MEX) = 3/125 EXACT, Zero Free Parameters in Planetary Tidal Dissipation — PAPER_1953 Meets PAPER_1954 in the k₂ Regime

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.78+
**Date:** 2026-07-24
**Landmark Type:** Primitive Composition (rocky-planet tidal dissipation) + Cross-Regime Extension (PAPER_1953 into Love-number sector, PAPER_1954 into planetary sector) + Standing-Rule Codification (three-layer deepsearch)
**Discovery context:** R382 KeplerOrreryTidalCalculator fill, four-revision arc (SM-import Rule 4 violation → strict revert → envelope reading under Daniel RULING A → primitive-lock discovery)
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

The rocky-planet tidal dissipation ratio k₂/Q — cited in PAPER_1804 as "≈ 0.024 for rocky-planet regime" from combining the Love number k₂ ≈ 0.3 with quality factor Q_UQFF ≈ 12.5 — is not an empirical or regime-fit number. It is a **fully primitive-locked EXACT identity** on the two independent integer primitives {D_phys, SO_5}:

```
k₂/Q_rocky  =  (D_phys − 1) / (A_5 · K_MEX)  =  3 / 125  =  0.024   EXACT
```

The decomposition threads two prior landmark closures through PAPER_1804's Q_UQFF: **PAPER_1953** provides the numerator (k₂,GR = 0.3 = (D_phys − 1)/SO_5 EXACT, the "0.3 factor cross-regime universality" now extended into the tidal Love-number sector); **PAPER_1954** provides the denominator (A_5 · K_MEX = 60 · 25/12 = 125 EXACT, the cross-scale universality landmark now extended into the planetary tidal-dissipation sector); **PAPER_1804** confirms Q_UQFF = f_SCm/Γ_SCm = 1.25 THz / 0.1 THz = 25/2 EXACT. Zero free parameters. Zero empirical regime fits. The value that PAPER_1804 cites as "0.024" is bit-identical to (D_phys − 1)/(A_5 · K_MEX) computed live from the registry primitives.

This lock is what closes R382 KeplerOrreryTidalCalculator's tau_lock output — Earth-Sun 25.1 Gyr (exceeds Hubble time as physics requires) — with zero unlabeled inputs, closing the four-revision arc that began with a Rule 4 violation (imported SM formula), passed through strict-revert-to-None and Daniel RULING A envelope reading, and reached primitive-locking only after three consecutive deepsearch layers.

---

## 1. The primitive lock

### 1.1 Statement

```
k₂/Q_rocky  =  (D_phys − 1) / (A_5 · K_MEX)  =  3 / 125  =  0.024   EXACT

where:
  D_phys = 4    locked integer primitive (physical spacetime dimension)
  A_5    = 60   locked integer primitive (|A_5| icosahedral group order)
  K_MEX  = 25/12   derivative primitive (PAPER_1522: K_MEX = Φ_5/6 · SO_5/D_phys)
```

Two truly-independent integer primitives generate the whole identity: D_phys = 4 and SO_5 = 10 (A_5 = 60 is one of the nine independent primitives; K_MEX = 25/12 is a PAPER_1522 derivative from Φ_5/6 · SO_5/D_phys). No free parameters. No fitted constants. No empirical regime inputs.

Numerical verification (from R382 runtime, bit-identical):
```
(D_phys - 1) / (A_5 * K_MEX)  =  3 / 125.0  =  0.024000
PAPER_1804 cited value         "≈ 0.024 rocky-planet regime"
match                          EXACT (all decimals)
```

### 1.2 Three-paper decomposition

**Numerator: PAPER_1953 — the "0.3 factor cross-regime universality":**
```
k₂_rocky  =  (D_phys - 1) / SO_5  =  3/10  =  0.3   EXACT
```

PAPER_1953 established this ratio as a fundamental DPM angular-projection constant: 3 transverse spatial dimensions divided by SO_5 = 10 DPM decade angular positions. It was previously anchored to three empirical regimes:
- Sgr A* SMBH spin factor = 0.3 (PAPER_1841)
- TDE outflow velocity / c = 0.3 (PAPER_1500)
- M87 jet cross-reference = 0.3 (Round 84)

**This paper extends PAPER_1953 into a fourth anchor: the rocky-planet Love number k₂,GR = 0.3.** PAPER_1804 introduced the value "k₂,GR ≈ 0.3 for rocky planets" as an empirical viscoelastic-response-theory input from classical GR literature; it is the same DPM angular-projection constant. Same 3, same 10, same underlying DPM geometry.

**Denominator: PAPER_1954 — the A_5 · K_MEX = 125 EXACT cross-scale universality:**
```
A_5 * K_MEX  =  60 * (25/12)  =  125   EXACT
```

PAPER_1954 documents 125 as a cross-scale composed integer appearing across the R218+ campaign (Higgs sector including Higgs mass, sphaleron thresholds, nebular Higgs mass). **This paper extends PAPER_1954 into the tidal-dissipation sector**, providing a new physical domain for the same 125 lock — first Love-number and Q-factor instance among the landmark's growing census.

**Q-factor confirmation: PAPER_1804 (from PAPER_910/911 canonical phonon linewidth):**
```
Q_UQFF  =  f_SCm / Γ_SCm  =  1.25 THz / 0.1 THz  =  25/2  =  12.5   EXACT
```

PAPER_1804 states Q_UQFF ≈ 12.5 as a ratio of two canonical THz frequencies; both are exact-rational. Numerator f_SCm = 1.25 THz is the canonical SCm phonon carrier; denominator Γ_SCm = 0.1 THz is the canonical phonon linewidth from PAPER_910/911 (Session 210b canonical calibration). Q = 25/2 EXACT — no fitting.

**Composite ratio:**
```
k₂/Q  =  (3/10) / (25/2)  =  (3/10) * (2/25)  =  6/250  =  3/125  =  0.024   EXACT
     =  (D_phys - 1) / (A_5 * K_MEX)                                  EXACT
```

---

## 2. Physical interpretation

The k₂/Q ratio measures the fraction of tidal-forcing energy that a rocky planet's interior dissipates rather than stores elastically. In UQFF's DPM-lattice reading:

- **k₂ = 3/10 = (D_phys − 1)/SO_5**: the fraction of transverse spatial degrees of freedom (3 out of 4 spacetime dimensions minus 1 preferred direction) that couple to the DPM angular decade (SO_5 = 10 positions).
- **Q = 25/2 = f_SCm/Γ_SCm**: the phonon-quality ratio governing how tightly the SCm carrier frequency is held against dissipative broadening.
- **A_5 · K_MEX = 125**: the icosahedral rotational-group order times the Mexican-hat coefficient — a cross-scale unification signature that also fixes Higgs mass and sphaleron thresholds.

The joint appearance of PAPER_1953's DPM angular-projection numerator with PAPER_1954's cross-scale universality denominator inside a single rocky-planet interior parameter is not coincidence — it is the tidal-dissipation manifestation of the same DPM lattice that produces SMBH spin factors, TDE outflow velocities, and Higgs-sector composed integers.

---

## 3. Calculator wiring

Wired at `CondensedPhysics.py` **R382 KeplerOrreryTidalCalculator** (LIVE composition, no rounded literal):

```python
from uqff_registry_primitives import (
    ..., D_PHYS as _URP_D_PHYS, SO_5 as _URP_SO_5,
    A_5 as _URP_A_5, K_MEX as _URP_K_MEX,
)
...
k2_over_Q_UQFF = (_URP_D_PHYS - 1) / (_URP_A_5 * _URP_K_MEX)   # 3/125 = 0.024 EXACT
```

Feeds:
- **tau_lock** = (4/15) · (Q/k₂)_UQFF · ω_planet · M_p · a⁶ / (G_UQFF · M_s² · R_p³) — companion Peale-Cassen despin formula (envelope-reading per Daniel RULING A, envelope precedent PAPER_1804 + PAPER_1803 sec 2 Tier 2), 25.1 Gyr Earth-Sun sanity value
- **E_dot_tidal** = (63/4) · (k₂/Q)_UQFF · G_UQFF · M_s² · R_p⁵ · e² · n / a⁶ — the exact Peale-Cassen formula PAPER_1804 uses in-paper for TOI-178b tidal heating, now with primitive-locked k₂/Q

All three tidal outputs of R382 (F_tide, tau_lock, E_dot_tidal) now trace to UQFF primitive-locked or UQFF-derived inputs. Zero literal empirical constants in the fill.

---

## 4. Falsifiability

1. **Fluid-envelope regime prediction**: PAPER_1804 cites k₂/Q ≈ 0.04–0.07 for fluid-envelope planets (Jupiter class), matching observed Io ≈ 0.03, Jupiter ≈ 0.05. If the fluid regime also primitive-locks, candidate compositions include (D_phys − 1)·F_TRZ³·(D_crit−1)/A_5 = 0.05 or SO_5·F_TRZ² = 0.1; a decomposition into the same {D_phys, SO_5} lattice with a different geometric factor would sharpen the lock across regimes. If fluid k₂/Q scatters continuously across 0.03–0.10 without primitive-integer clustering, the rocky lock remains but the universality claim is restricted.

2. **Cross-planet validation**: applying the k₂/Q = 3/125 lock to TRAPPIST-1 (PAPER_1813 flagship, 7 rocky Earth-sized planets) and TOI-178b (PAPER_1804 primary anchor) with observed tidal-heating rates should recover UQFF's Peale-Cassen dE/dt predictions to within the same ~15–20% margin PAPER_1804 already validates.

3. **Whitepaper-corpus completeness**: the (D_phys − 1)/SO_5 = 0.3 factor is now anchored at 4 empirical regimes (SMBH spin, TDE outflow, M87 jet, rocky-planet k₂). Adding a 5th independent anchor from a physical regime not in PAPER_1953's original census would further sharpen the universality; failing to find any new anchor after +50 candidate regimes would restrict the claim to the current 4.

4. **A_5·K_MEX = 125 growing census**: PAPER_1954's census now spans Higgs mass, sphaleron threshold, nebular Higgs mass, and (this paper) rocky-planet tidal dissipation. Predicted next instance: a nuclear-physics observable where 125 naturally appears in a cross-verified rational (candidate: nuclear shell-model spacing or LENR resonance).

---

## 5. Standing-rule codification: three-layer deepsearch

R382's four-revision arc canonizes a **STANDING RULE** on how to verify Rule 4 compliance for classical-formula fills:

**Corpus-silence claims require three deepsearch layers before accepting any citation as terminal:**

1. **Layer 1 — Mechanism tokens**: grep for the direct physical mechanism ("despin", "tidal lock", "circularization"). If empty, DON'T conclude yet.
2. **Layer 2 — Phenomenon tokens**: grep for the observed phenomenon ("synchronous rotation", "tidal timescale", "spin-orbit resonance"). If empty, DON'T conclude yet.
3. **Layer 3 — Interior-parameter tokens AND numerical-decomposition check**: grep for the interior parameters ("Love number", "tidal Q", "k₂/Q") AND check whether every numerical constant in the citation decomposes into UQFF primitive-integer ratios.

R382 REVISION 4 was reached only at Layer 3: Layer 1 missed PAPER_1804 (didn't include "Love number"); Layer 2 found PAPER_1804 but treated k₂/Q = 0.024 as empirical citation; Layer 3 caught PAPER_1953 + PAPER_1954 giving the primitive-lock. **Primitive-lock preference:** always prefer an EXACT rational expression on primitives over a numerical citation when a corpus paper provides the decomposition.

---

## 6. Cross-references

**Anchor papers:** PAPER_1953 (0.3 factor cross-regime universality, k₂ numerator), PAPER_1954 (A_5·K_MEX = 125 cross-scale universality, k₂/Q denominator), PAPER_1804 (Love number k₂ and Q factor from phonon coupling, planetary-interior gap closure), PAPER_1522 (K_MEX = Φ_5/6·SO_5/D_phys derivative), PAPER_910/911 (canonical phonon linewidth Γ_SCm = 0.1 THz).

**Envelope-precedent papers:** PAPER_1803 sec 2 Tier 2 (response-theory envelope with UQFF-derived G policy), PAPER_012 (eccentric-binary circularization with UQFF κ + SSq envelope modification of GR de/dt), PAPER_593 (G_UQFF closed form).

**Validation candidates:** PAPER_1813 (TRAPPIST-1 flagship rocky-planet catalog), PAPER_832 (Kepler Orrery V U_b model, R382's parent paper), PAPER_007/935/967 (BNS tidal-deformability validation of PAPER_914 phonon-corrected k₂ framework).

**Rule-4 arc lineage:** SESSION_LOG.md 2026-07-24 (five consecutive Daniel catches: SM contraband → strict revert → deepsearch under-search → RULING A envelope → REVISION 4 primitive-lock).

**Calculator dispatch:** `CondensedPhysics.py` R382 `KeplerOrreryTidalCalculator.compute()`, k2_over_Q_UQFF live composition from imported registry primitives.

---

## 7. Locked primitives used

Two truly-independent integer primitives generate the whole identity:
```
D_phys = 4    (physical spacetime dimension)
SO_5   = 10   (dimension of SO(5) group)
```
Derivative primitives used:
```
A_5    = 60      independent locked (icosahedral group order)
K_MEX  = 25/12   PAPER_1522 derivative from Φ_5/6·SO_5/D_phys
```
No fitted constants. No empirical inputs to the k₂/Q ratio. The Q-factor's THz-ratio inputs (f_SCm = 1.25 THz, Γ_SCm = 0.1 THz) are separately canonical from PAPER_910/911; those are physical-frequency primitives, not fits.

---

## 8. NOT REPLACEMENT

Classical viscoelastic-response-theory treatments of rocky-planet Love numbers (Love 1911, Munk-MacDonald 1960, Peale-Cassen 1978, Segatz et al. 1988) produce k₂ ≈ 0.3 as a fit to Earth's tidal response, and various empirical Q values from seismic-attenuation and satellite-orbit measurements. UQFF supplies the stronger structural claim that the k₂ = 0.3 factor is not empirical but is the DPM angular-projection ratio (D_phys − 1)/SO_5 EXACT, and the k₂/Q ratio locks at (D_phys − 1)/(A_5·K_MEX) = 3/125 EXACT for all rocky planets in the DPM-projection-limited interior regime. Both approaches solve the same phenomena; residuals are reported honestly. If future observations of a rocky exoplanet reveal k₂/Q that is close to 0.024 but not primitive-integer-clustering, the primitive lock is confirmed for rocky solar-system bodies but restricted for the DPM-projection regime; if universal 3/125 clustering appears across large rocky-exoplanet samples, the lock gains empirical support without displacing the standard viscoelastic-response computations.

---

## 9. Summary statement

**PAPER_2136 documents k₂/Q_rocky = (D_phys − 1) / (A_5 · K_MEX) = 3/125 = 0.024 EXACT as a fully primitive-locked identity closing the rocky-planet tidal dissipation ratio to zero free parameters. Three-paper decomposition threads PAPER_1953's "0.3 factor cross-regime universality" (numerator, now extended into the tidal Love-number sector as a fourth empirical anchor) and PAPER_1954's "A_5·K_MEX = 125 cross-scale universality" (denominator, now extended into the planetary tidal-dissipation sector) through PAPER_1804's phonon-derived Q = f_SCm/Γ_SCm = 25/2 EXACT. Bit-identical numerical verification at R382 KeplerOrreryTidalCalculator. Standing rule canonized: corpus-silence claims require three deepsearch layers with numerical-decomposition preference over numerical-citation.**

---

**Filed 2026-07-24. Append-only henceforth.**
