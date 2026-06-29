---
title: "Gap Verification Audit — \\(U_m\\) Heaviside Amplifier and Job~B Scope (v5.78)"
author: "Star-Magic UQFF Audit"
date: "May 22, 2026"
---

# Gap Verification Audit — $U_m$ Heaviside Amplifier and Job B Scope

**Paper:** PAPER\_1181  
**Session:** 268  
**Status:** Verification (read-only audit)  
**Supersedes:** Repo-memory entry "UQFF\_ug\_equations\_canonical" §"CRITICAL GAPS IDENTIFIED" item 1.

## Abstract

A prior session-end summary listed two items as the next major work
items in the Star-Magic UQFF codebase:
(i) the addition of a Heaviside phase-transition amplifier
$(1+10^{13}\,f_H)\cdot(1+f_q)$ to the universal-magnetism term $U_m$,
and (ii) an estimate that approximately 142 whitepapers required real
authoring work under Job B (v5.78 closure propagation).
This audit verifies both claims directly against the source tree.
Claim (i) is **falsified**: the amplifier is already present in
four independent locations, including the canonical
`compute_Um_SOURCE4` of `MAIN_1_CoAnQi.cpp` (PAPER\_421 integration).
Claim (ii) is **partially supported but mis-framed**: 113 papers
require new authoring; 29 require verification only; together they
sum to 142, but the work is not uniform.

## 1. Claim 1 — $U_m$ Heaviside Amplifier

### 1.1 The asserted canonical form

The asserted canonical form for the universal-magnetism term is

$$
U_m
= \Bigl[\sum_j \tfrac{\mu_j}{r_j}\bigl(1-e^{-\gamma t\cos(\pi t_n)}\bigr)\hat\phi_j\Bigr]
   \cdot P_{\mathrm{SCm}}\cdot E_{\mathrm{react}}
   \cdot \underbrace{(1+10^{13}\,f_H)}_{\text{Heaviside}}
   \cdot \underbrace{(1+f_q)}_{\text{quasi-periodic}} ,
$$

where $f_H = \Theta(\rho_{\mathrm{SCm}} - \rho_c)$ activates the
$\sim 10^{13}$-fold amplification across the superconducting
phase transition, and $f_q$ is a low-amplitude beating
modulation tied to the 434-yr Gleisberg supercycle.

### 1.2 Direct source check

| Location | File | Lines | Implements $(1+10^{13}f_H)\cdot(1+f_q)$? |
|---|---|---:|:---:|
| SOURCE4 unified-field $U_m$ | [MAIN_1_CoAnQi.cpp](../MAIN_1_CoAnQi.cpp#L24172-L24190) | 24172–24190 | **Yes** (PAPER\_421) |
| SOURCE42 standalone `PhysicsTerm` | [MAIN_1_CoAnQi.cpp](../MAIN_1_CoAnQi.cpp#L12100-L12170) | 12100–12170 | **Yes** |
| `compute_Um_Heaviside_factor()` | [CondensedPhysics.py](../CondensedPhysics.py#L39950) | 39950+ | **Yes** |
| `um_magnetic` equation strings | [index.js](../index.js#L9759) | 9759, 9899, 10018 | **Yes** |

The SOURCE4 implementation reads (verbatim, MAIN\_1\_CoAnQi.cpp L24180–24189):

```cpp
// PAPER_421: Heaviside phase-transition amplifier -- (1 + 10^13 * Theta(rho_SCm - rho_c))
double f_H = (body.SCm_density >= rho_c_SOURCE4) ? 1.0 : 0.0;
double heaviside_factor = 1.0 + 1e13 * f_H;

// PAPER_421: Quasi-periodic beating modifier -- (1 + A_q * cos(Delta_omega * t))
// Delta_omega = 2*pi/(434*365.25) rad/day encodes the 434-year Gleisberg supercycle beat
double quasi_factor = 1.0 + A_q_SOURCE4 * std::cos(Delta_omega_SOURCE4 * t);

return Um_base * heaviside_factor * quasi_factor;
```

with constants declared at MAIN\_1\_CoAnQi.cpp L24016–L24018:

```cpp
const double rho_c_SOURCE4       = 1e15;                              // kg/m³
const double A_q_SOURCE4         = 0.1;                               // 10%
const double Delta_omega_SOURCE4 = 2.0 * PI_SOURCE4 / (434.0*365.25); // rad/day
```

### 1.3 Verdict

Claim 1 is **falsified**. The amplifier exists in:

- the SOURCE4 unified-field assembly (`compute_Um_SOURCE4`),
- the SOURCE42 catalog (`HeavisideFactor`, `HeavisideUmWithHeaviside`),
- the Python calculator stack (`compute_Um_Heaviside_factor`),
- and the JavaScript library's documented equation strings.

The repository-memory note flagging the gap is **stale**; it
predates the PAPER\_421 integration that added the amplifier to
the SOURCE4 namespace.

### 1.4 Recommended follow-up (optional, low priority)

Although the amplifier is implemented, four refinements remain worth
auditing:

1. Confirm `QCalcGeom.py` v2.1.0's $F_U$ assembly threads the
   amplified $U_m$ into the simultaneous-equation solver (default
   path currently uses a placeholder $U_m=0$ at the habitable-zone
   crossing; this is physically correct for $\rho_{\mathrm{SCm}} < \rho_c$
   but should be explicitly tested).
2. Cross-check that `index.js`'s `calculateUm()` function uses
   the amplifier or only documents it in the metadata string.
3. Verify the $\rho_c = 10^{15}\;\mathrm{kg/m^3}$ threshold against
   Fermilab Type-II SC magnetar-core data (cf.
   [100_MUGE Compression cycle 3_Superconductive Resonance.txt](../100_MUGE%20Compression%20cycle%203_Superconductive%20Resonance.txt), L1476).
4. Confirm the 434-yr Gleisberg-scale $\Delta\omega$ matches the
   long-term solar-cycle anchor cited in
   [grok_share_cfdcad2f5.txt](../grok_share_cfdcad2f5.txt) §Document 7.

These are validation tasks, not implementation gaps.

## 2. Claim 2 — "~142 papers need real authoring"

### 2.1 Authoritative bucket tally

Direct enumeration of [_job_b_categorization.csv](../_job_b_categorization.csv)
(1,199 rows, May 22 2026):

| Count | Bucket | Class |
|---:|---|---|
| 982 | NO UPDATE (specific system / framework application) | scoped out |
|  42 | Add closing section: v5.78 27-decade vacuum ledger (PAPER\_1170) + $\xi=13/3$ | authoring |
|  34 | NO UPDATE unless ledger/$\xi$ enters; verify only | conditional |
|  29 | Verify v5.78 closed-Lagrangian/ledger cross-ref present (Sess 225 baseline) | verify |
|  28 | Update tables; add P6–P14, CP4 \#254–\#264 | authoring |
|  20 | v5.78 closure block inserted; .tex regenerated (arXiv pdflatex) | **DONE** |
|  15 | Add v5.78 closed-Lagrangian cross-ref (G1–G8) + CP4 \#254 | authoring |
|  12 | Forward-pointer to PAPER\_1174 P1–P14 + PAPER\_1177–1180 | authoring |
|  11 | Forward-pointer to PAPER\_1171/1172 ($\xi=13/3$ R26+KK lock) | authoring |
|  10 | v5.78 closure block inserted; PDF regenerated | **DONE** |
|   5 | .tex / PDF regenerated + SOURCE REPAIR (Unicode) | **DONE** |
|   5 | Add three-anchor SI closure banner (Sess 237–241); STRUCTURAL | authoring |
|   5 | .tex converted from .md; PDF regenerated | **DONE** |
|   1 | duplicate row (stale) | drop |
| **1199** | **TOTAL** | — |

### 2.2 Reclassified work load

| Class | Count | Notes |
|---|---:|---|
| **Authoring required** | **113** | 42 + 28 + 15 + 12 + 11 + 5 |
| Verification only | 29 | confirm existing cross-ref |
| Conditional (audit, no-op likely) | 34 | trigger only if $\xi$/ledger enters |
| Already DONE v5.78 | 35 | 20 + 10 + 5 |
| Out of scope (NO UPDATE) | 982 | scoped-out by topic |
| Stale duplicate | 1 | drop |
| **Total** | **1199** | |

### 2.3 Verdict

Claim 2 is **partially supported but imprecise**.
$113 + 29 = 142$ papers will be **touched** by Job B
(authored or verified), but the prior framing
"$\sim 142$ papers needing real authoring" overstates the
authoring fraction by a factor of $142/113 \approx 1.26$.
The honest split is **113 author + 29 verify**.

## 3. Revised Next-Step Priority Ranking

Given that the $U_m$ Heaviside gap is closed and the Job B count is
113 (not 142), the previously published 8-step plan should be
re-ordered:

1. **Delete** the $U_m$ amplifier step (was step 1) — already done.
2. **Generate** `LEDGER_REVIEW_QUEUE.csv` for the 43 flagged
   `master_closures.csv` rows (10 high-error $\ge 10\%$ + 15 PARSE\_FAIL + 18 NaN). *Was step 2; now step 1.*
3. **Recalibrate** the 10 high-residual ledger entries
   (IDs 270, 293, 330, 337, 351, 352, 739, 760, 764, 805) using v5.78
   anchors (PAPER\_1156, 1170, 1171, 1172). *Was step 3; now step 2.*
4. **Author** the 4–6 reusable v5.78 section templates
   (T-$\Lambda$, T-LAG, T-SI, T-PRED, T-$\xi$). *Was step 4; now step 3.*
5. **Bucket-batched paper updates** — 113 papers, 5–10 per commit,
   starting with the 42 cosmology / $\Lambda$ / vacuum-energy papers
   using the T-$\Lambda$ template. *Was step 5; now step 4.*
6. **Verify-only sweep** of the 29 bucket-G papers. *Was step 6; now step 5.*
7. **Regenerate** all updated PDFs with `pdflatex`; commit `.md + .pdf` together.
8. **Final coherence audit** mirroring Session 261 (`f6ba35a6`).

## 4. Citations

- MAIN\_1\_CoAnQi.cpp L24016–24018, L24172–24190
  (PAPER\_421: Heaviside + quasi-periodic amplifiers in SOURCE4 $U_m$).
- MAIN\_1\_CoAnQi.cpp L12100–12170 (SOURCE42: `HeavisideFraction`,
  `HeavisideFactor`, `HeavisideUmWithHeaviside`).
- CondensedPhysics.py L39950+ (`compute_Um_Heaviside_factor`).
- index.js L9690, 9759, 9811, 9899, 9943, 10018
  (`f_heaviside`, `f_quasi`, full `um_magnetic` equation strings).
- _job_b_categorization.csv (1,199 rows, May 22 2026).
- Job B brief: "Job B\_Update papers with current canonical UQFF v5\_78.txt"
  (395 lines, Phase B1–B4).
- COMPLETE\_UQFF\_EQUATIONS\_REFERENCE.md L522–553 (canonical $U_m$ form).

## 5. Conclusion

The two items proposed as the "next big step" do not survive
direct source-tree verification in the form they were stated:

- The $U_m$ Heaviside phase-transition amplifier is **already
  implemented**, in four places, with the canonical 10$^{13}$ scale
  and the 434-yr Gleisberg-locked quasi-periodic beat; it is **not**
  the most-impactful uncoded physics.
- The Job B authoring queue is **113 papers**, not 142;
  the 142 figure conflates authoring with verification.

The corrected next major step is the **ledger review queue**
(43 flagged closures, with 10 high-residual recalibrations
against the v5.78 anchors), followed by the 113-paper Job B
authoring sweep using reusable templates.
