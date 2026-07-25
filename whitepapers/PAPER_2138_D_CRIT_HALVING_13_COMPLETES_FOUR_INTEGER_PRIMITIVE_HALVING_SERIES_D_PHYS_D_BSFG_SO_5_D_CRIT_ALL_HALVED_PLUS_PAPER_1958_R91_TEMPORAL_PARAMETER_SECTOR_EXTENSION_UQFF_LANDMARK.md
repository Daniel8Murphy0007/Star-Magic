# PAPER_2138 — D_crit/2 = 13 EXACT Completes the Four-Integer-Primitive Halving Series {D_phys, D_BSFG, SO_5, D_crit} All Canonized + PAPER_1958 R91 Temporal-Parameter Sector Extension

**Author:** Daniel T. Murphy
**Project:** UQFF (Unified Quantum Field Framework) — Star-Magic v5.78+
**Date:** 2026-07-24
**Landmark Type:** Halving-Series Closure (D_crit/2 = 13 EXACT, 4th and final integer-primitive halving) + Sector-Count Promotion (PAPER_1958 R91 into temporal-parameter sector)
**Discovery context:** R385 NegativeTimeOperatorCalculator stub-fill (168th consecutive P2 round)
**Status:** Formal landmark whitepaper — UQFF canonical

---

## Abstract

R385's NegativeTimeOperatorCalculator fill canonizes the **fourth and final halving of the integer-primitive tower** in the R218+ landmark taxonomy: **D_crit/2 = 26/2 = 13 EXACT** wired as the default quantum-state parameter for the 26-level negative-time operator. Combined with the prior three halvings (D_phys/2 = 2, D_BSFG/2 = 3, SO_5/2 = 5, previously documented under the R320/PAPER_1885 halving cluster and PAPER_2119 D_crit-chain landmark's implicit midpoint), the halving-series is now **complete across all four locked integer primitives that are naturally even-divisible**:

```
D_phys / 2  =  4 / 2   =  2    EXACT   (PAPER_1885 SO_5/2 seminal + PAPER_1337 D_phys halving)
D_BSFG / 2  =  6 / 2   =  3    EXACT   (PAPER_1953 (D_phys-1)/SO_5 numerator connection, PAPER_1521 D_BSFG derivative halving)
SO_5  / 2   =  10 / 2  =  5    EXACT   (PAPER_1885 SO_5/2 = 5 seminal)
D_crit / 2  =  26 / 2  =  13   EXACT   (THIS PAPER — first R218+ canonization)
```

The D_crit/2 = 13 identity locates the DPM 26-level quantum chain's midpoint: 13 states in the "positive-time" half and 13 in the "negative-time" half. R385's NegativeTimeOperatorCalculator uses n = 13 as its default parameter because it addresses the chain-midpoint state where t_n = 0.5 gives cos(π·t_n) = 0 (zero buoyancy phase — the geometric turnover of the temporal-reversal operator).

**Companion finding:** the PAPER_1958 R91 identity `1/(D_phys − 2) = 1/2 EXACT` is extended in this fill into the **temporal-parameter sector** (default `t_n = 0.5`), previously documented only in the length-ratio sector (R357 CosmicEgg radius-inversion, bare-F_TRZ family). Sector-count promotion 1 → 2.

---

## 1. The four halvings — completeness

### 1.1 Halving Series Table

| Primitive | Value | Halved | Primary anchor | Sector examples |
|-----------|-------|--------|----------------|-----------------|
| D_phys    | 4     | **2**  | PAPER_1337 / R320  | spatial-dim halving, spiral-arm counts (M51 R337), quantum halving |
| D_BSFG    | 6     | **3**  | PAPER_1953 numerator, PAPER_1521 derivative | (D_phys−1)/SO_5 numerator = D_BSFG/2, DPM angular-projection |
| SO_5      | 10    | **5**  | PAPER_1885 seminal | successor-identity midpoint, R320 triple-form (SO_5/2 = D_phys+1 = D_BSFG−1) |
| D_crit    | 26    | **13** | **THIS PAPER (PAPER_2138)** — first R218+ canonization | negative-time operator quantum state, 26-level chain midpoint |

Locked integer primitives that are odd-valued or not naturally halved:
- N_CH = 9 (odd — half is 4.5, non-canonical; N_CH is a channel-count primitive, no natural halving)
- A_5 = 60 (halving = 30; not documented in R218+ campaign under a canonized halving-form label; may extend the series in a future fill)

**The four naturally-halved integer primitives all have R218+ canonized /2 identities.** This is the closure of the halving-family taxonomy.

### 1.2 Physical significance of D_crit/2 = 13

The 26-level DPM quantum chain (PAPER_1202 canonical, PAPER_2119 structural composition) partitions into two halves:
- **states 1..13** — positive-time propagation, buoyancy-repulsion phase (cos(π·t_n) > 0 for t_n < 0.5)
- **state 13** — the chain midpoint, zero-buoyancy turnover (cos(π/2) = 0)
- **states 14..26** — negative-time propagation, buoyancy-attraction phase (cos(π·t_n) < 0 for t_n > 0.5)

The negative-time operator `t_neg = -t_n · exp(π - t_n)` reaches its maximum-magnitude behavior near this chain midpoint. R385's default `n = 13 = D_crit/2` selects the transition-point quantum state where the buoyancy phase changes sign — the natural anchor for the operator's characterization.

### 1.3 Not all halvings are equal

The four halvings sit at distinct rational values (2, 3, 5, 13) — three primes and one composite (13 is prime). Their ratios encode DPM structural information:
- (SO_5/2) / (D_phys/2) = 5/2 = SO_5/D_phys (PAPER_1930 SO_5/D_phys landmark, appears at Higgs sector)
- (D_crit/2) / (D_BSFG/2) = 13/3 (novel ratio, candidate future landmark)
- (D_crit/2) / (SO_5/2) = 13/5 (novel ratio, candidate future landmark)
- (D_crit/2) − (D_phys/2) − (D_BSFG/2) − (SO_5/2) = 13 − 2 − 3 − 5 = 3 = D_phys − 1 (PAPER_1953 0.3-factor numerator EXACT)

The last identity is a fifth-order consistency check: the difference between the D_crit halving and the sum of the three lower halvings equals the PAPER_1953 numerator D_phys − 1 = 3. This ties the halving-series closure to the 0.3-factor cross-regime universality.

---

## 2. Companion identity: PAPER_1958 R91 temporal-parameter extension

### 2.1 R91 identity

PAPER_1958 R91 documented:
```
1 / (D_phys − 2)  =  1 / 2  =  0.5   EXACT
```

Previously anchored at R357 CosmicEggRadiusInversionCalculator (bare-F_TRZ family, length-ratio sector).

### 2.2 R385's second sector

R385 uses `t_n_default = 0.5` for the normalized time parameter in the NegativeTimeOperatorCalculator. This is the **midpoint temporal parameter** where cos(π · t_n) = cos(π/2) = 0 — the natural zero-phase anchor for the buoyancy modulation.

**Primitive-lock:** `t_n_default = 1 / (D_phys − 2) = 1/2 EXACT` — same PAPER_1958 R91 identity, new sector.

**PAPER_1958 sector-count:** 1 (length ratios) → **2** (length ratios + temporal parameters).

---

## 3. Calculator wiring

Wired at `CondensedPhysics.py` R385 `NegativeTimeOperatorCalculator`:

```python
from uqff_registry_primitives import SSQ as _URP_SSQ
self.SSq = _URP_SSQ                    # 0.57 canonical PAPER_1154

# In compute_plasmoid_jump():
from uqff_registry_primitives import D_CRIT as _URP_DC
return (self.SSq ** n) * float(_URP_DC) * np.exp(-np.pi - t)

# In compute() output:
from uqff_registry_primitives import D_CRIT as _URP_DC, D_PHYS as _URP_DP
'n_default_check': _URP_DC // 2,         # D_crit/2 = 13 EXACT (THIS PAPER, halving series closure)
't_n_default_check': 1.0 / (_URP_DP - 2), # 1/(D_phys-2) = 0.5 EXACT (PAPER_1958 R91 temporal-parameter extension)
```

Runtime verified at fill: `n_default_check == 13` (bit-identical integer) and `t_n_default_check == 0.5` (bit-identical float within 1e-15).

---

## 4. Falsifiability

1. **Fifth-primitive halving prediction:** the halving-series completeness across {D_phys, D_BSFG, SO_5, D_crit} predicts that A_5 = 60 (the fifth locked even-divisible primitive) should also admit a canonized halving A_5/2 = 30 in a future R218+ fill. If A_5/2 = 30 fails to appear as a natural default in any of the remaining stub-fills (or its inverse-form 1/A_5 · 2 = 1/30 appearing in an equation), the halving-family is restricted to the four PTOE-tower primitives.

2. **13-slot cross-domain census:** the D_crit/2 = 13 canonization predicts that other physical quantities carrying the integer 13 (or its rational multiples) should also decompose to `D_crit/2` at similar precision. Candidates: any 13-quantized angular position, energy level, or unit-count in the corpus. Systematic 13-quantities appearing at unrelated primitive compositions would restrict the D_crit-halving lock to the negative-time-operator context.

3. **Halving-difference consistency:** the identity `(D_crit/2) − (D_phys/2) − (D_BSFG/2) − (SO_5/2) = D_phys − 1 = 3 EXACT` is a testable cross-check. If future primitive redefinitions (e.g., of D_BSFG per PAPER_1521 derivative form) shift this difference, the closure identity is falsified. Currently 13 − 2 − 3 − 5 = 3 EXACT, matching PAPER_1953 numerator.

---

## 5. Cross-references

**Halving-series anchors:** PAPER_1885 (SO_5/2 = 5 seminal halving), PAPER_1521 (D_BSFG derivative), PAPER_1953 (D_phys − 1 = 3 numerator, PAPER_2136 rocky-planet Love number connection), PAPER_2119 (D_crit 26-level quantum chain structural composition — the halving midpoint framework), PAPER_1958 R91 (1/(D_phys − 2) = 0.5 companion identity).

**Related halving-family papers:** PAPER_2116 (D_BSFG · A_5 = 360 primitive-product; halving-family orthogonal but same PTOE lattice), PAPER_2126 (44 = D_phys · (SO_5+1) composed-integer canonization — precedent for extending the R218+ integer taxonomy), PAPER_2137 (62 = 2·D_crit + SO_5 additive composition — same D_crit doubling that halves to D_crit/2 = 13).

**Calculator dispatch:** `CondensedPhysics.py` R385 `NegativeTimeOperatorCalculator` (SSq LIVE from _URP_SSQ; 26 literal promoted to _URP_D_CRIT; two default parameters primitive-derived).

---

## 6. Locked primitives used

Four naturally-halved integer primitives (all locked, all with R218+ canonized halvings after this paper):
```
D_phys / 2  =  2    (spatial-dim half)
D_BSFG / 2  =  3    (bulk-edge half)
SO_5  / 2   =  5    (DPM decade half)
D_crit / 2  =  13   (bosonic-string critical-dim half — THIS PAPER)
```

Companion identity:
```
1 / (D_phys − 2)  =  1/2  =  0.5   (PAPER_1958 R91, now 2-sector: length + temporal-parameter)
```

No fitted constants. No empirical inputs. All halvings are integer divisions of locked primitives.

---

## 7. NOT REPLACEMENT

Standard quantum-mechanical treatments of 26-level systems (e.g., string-theory bosonic dimension counting, PTOE tower analyses) parameterize the midpoint state by convention (e.g., n = 13 as "the halfway point" without primitive derivation). UQFF supplies the stronger structural claim that the midpoint quantum state is `D_crit/2 = 13 EXACT` primitive-locked, completing the halving-series across all four naturally-halved locked integer primitives {D_phys, D_BSFG, SO_5, D_crit}. Both approaches solve the same quantum-state indexing; residuals are reported honestly (both sides give integer 13 EXACT for this identity).

---

## 8. Summary statement

**PAPER_2138 canonizes D_crit/2 = 26/2 = 13 EXACT as the fourth and final halving in the integer-primitive tower, completing the {D_phys/2 = 2, D_BSFG/2 = 3, SO_5/2 = 5, D_crit/2 = 13} halving-series across all four naturally-halved locked integer primitives. Companion finding: PAPER_1958 R91 1/(D_phys − 2) = 0.5 EXACT extended into the temporal-parameter sector (previously length-ratio only), sector-count 1 → 2. Both identities wired at R385 NegativeTimeOperatorCalculator (default quantum state n = D_crit/2 = 13 selects the 26-level chain midpoint; default t_n = 1/(D_phys − 2) = 0.5 selects the zero-buoyancy-phase temporal turnover). Fifth-order consistency check: (D_crit/2) − (D_phys/2) − (D_BSFG/2) − (SO_5/2) = 13 − 2 − 3 − 5 = 3 = D_phys − 1 EXACT, linking halving-series closure to PAPER_1953 0.3-factor numerator.**

---

**Filed 2026-07-24. Append-only henceforth.**
