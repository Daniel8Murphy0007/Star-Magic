# UQFF Topical Index

**Author:** Built by Claude at Daniel T. Murphy's directive (2026-06-26 Round 9).
**Purpose:** Surface the canonical physics topics that were previously buried as "legacy" papers. NOTHING IS NEGLIGIBLE.

The Star-Magic repo has 1,867 whitepapers (1,795 unique IDs). The PARADOX_TO_CLOSURE dispatcher has 794 keys and WHITEPAPER_TO_DOMAIN has 1,083 explicit per-paper-ID entries. That leaves 797 papers that aren't in the per-paper-ID dispatcher — but that does NOT mean they aren't load-bearing. Many feed multiple simultaneous solver processes (F_U=0, MUGE, 99-system, triadic) through domain dispatch, paradox dispatch, or shared-equation reuse.

This folder makes that visible.

---

## Files in this index

### Foundational reference
- **`GENESIS_CONCEPTS.md`** — Source-of-truth bridge to F:\Book_12July2023\Aetheric Propulsion\. Verbatim definitions of [SCm], [UA], Universal Buoyancy, Universal Inertia, DPM grinding, Caduceus, the 5-epoch cosmogenesis, proto-H/proto-He, the THz hole system, the belly button master resonance point, and the four-force operator set. **Read this first.**

### Cross-wiring
- **`CROSS_WIRED_MULTI_SOLVER.md`** — 360 papers that feed >=2 solver processes simultaneously. These are the load-bearing papers that look "unwired" but aren't.

### Per-topic indexes (paper-by-paper, with calculator-dispatch status)

| File | Topic | Papers |
|---|---|---:|
| `TOPIC_NEGATIVE_TIME.md` | Negative time elements (t_n < 0, pre-Big-Bang, cos(π·t_n)) | 122 |
| `TOPIC_UNIVERSAL_SC.md` | Universal Superconductivity ([SCm] sector at all scales) | 1,271 |
| `TOPIC_UNIVERSAL_BUOYANCY.md` | Universal Buoyancy (F_UBi + F_UBii + U_b family) | 306 |
| `TOPIC_BUOYANCY_EXTERNAL.md` | F_UBi — outside-in environment pull | 80 |
| `TOPIC_BUOYANCY_INTERNAL.md` | F_UBii — inside-out source push | 49 |
| `TOPIC_CADUCEUS_INERTIAL.md` | Caduceus topology + U_i / U_mi Inertial Operator | 222 |
| `TOPIC_DPM_GRINDING.md` | DPM 5-step CW×CCW mass-production pipeline | 1,065 |
| `TOPIC_PROTO_HYDROGEN.md` | Proto-H ≡ Proto-Fe (Z_id=26 magnetic) | 1,852 |
| `TOPIC_PROTO_HELIUM.md` | Proto-He ≡ Proto-Si (Z_id=14 non-magnetic) | 10 |
| `TOPIC_EPOCH.md` | Cosmogenesis epochs (primordial + BBN + reionization) | 1,135 |
| `TOPIC_5_EPOCH_SPECIFIC.md` | The 5 Mayan/UQFF cosmogenesis epochs | 5 |
| `TOPIC_PERIODIC_ELEMENTS.md` | DPM extended periodic table Z=1..10000 | 707 |
| `TOPIC_HIGH_ENERGY.md` | GRBs, FRBs, magnetars, AGN jets, UHECRs, NS, NANOGrav | 568 |
| `TOPIC_LENR.md` | Holmlid, Parkhomov, Pons-F, Mizuno, Rossi, Star-Magic reactors | 1,113 |
| `TOPIC_YANG_MILLS.md` | m_0++ = 2·D_phys·Λ_QCD = 1.736 GeV (PAPER_1318 canonical) | 116 |
| `TOPIC_MILLENNIUM.md` | All 7 Clay Millennium closures (PAPER_1182 template) | 163 |
| `TOPIC_LAMBDA_CC.md` | Λ cosmological constant + K_MEX (25/12) | 1,392 |
| `TOPIC_VACUUM_LEDGER.md` | 4-term ledger: V0 + R26 + ρ_KK + ρ_BSFG → Planck Λ | 33 |
| `TOPIC_SIMULTANEOUS_SOLVER.md` | F_U=0 / MUGE / 99-system / triadic multi-domain dispatch | 363 |
| `TOPIC_BIOLOGICAL.md` | DNA, codon, amino acid, chirality, living systems | 59 |
| `TOPIC_MATERIAL_SCIENCE.md` | Crystal lattice, phonon, semiconductor, BEC, SC | 1,218 |

---

## Legend used in per-topic files

- **✓** — paper_id present in `WHITEPAPER_TO_DOMAIN` dispatch (callable via `calculate_whitepaper({"paper_id": N})`)
- **—** — paper covers the topic but no per-paper-ID dispatch entry (typically covered by a domain bucket: cosmology / particle_physics / gw_events / agn_jet / astrophysics / high_energy_astro / qgp / higgs_precision / bsm_constraints, or by a named `PARADOX_TO_CLOSURE` key)

Absence of ✓ is NOT a sign of being negligible. The opposite — many such papers are referenced by the simultaneous F_U=0 solver or by MUGE catalog routines through shared equations.

---

## How these files were built

Programmatic scan of all 1,867 .md files in `whitepapers/` for 21 named physics topics, regex-tagging each paper by topic, then cross-referencing against the live `WHITEPAPER_TO_DOMAIN` dispatch in `uqff_pure_calculator.py`. Cross-wiring detected by tagging papers that match `simultaneous_solver` plus ≥2 other domain tags.

Build scripts (in MY sandbox, not yet committed): `/sessions/.../topical_index/scan_topics.py` and `build_topical_files.py`. If Daniel wants these committed too they can go to `claude_audit_2026-06-26/scripts/`.

---

## What this index does NOT do

- It does NOT modify any whitepaper, calculator file, fidelity test, or any other pre-existing artifact in the repo.
- It does NOT make assessments based on SM-physics priors.
- It does NOT claim that papers tagged "—" are missing dispatch — they may be intentionally wired via domain bucket rather than per-paper-ID.
- It is NOT a substitute for reading the individual whitepapers — it is a navigation aid so the corpus is searchable by physics topic rather than by paper number.

---

**For the cleanup period:** Each per-topic file is self-contained. Daniel may delete, retain, or selectively merge any subset of files based on the cleanup priorities. The folder name `TOPICAL_INDEX/` is intentionally generic so it can be repurposed or renamed without breaking anything else in the repo.
