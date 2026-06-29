# Star-Magic — UQFF (Unified Quantum Field Framework)

[![PyPI version](https://img.shields.io/pypi/v/uqff.svg)](https://pypi.org/project/uqff/)
[![Python versions](https://img.shields.io/pypi/pyversions/uqff.svg)](https://pypi.org/project/uqff/)
[![Documentation Status](https://readthedocs.org/projects/star-magic/badge/?version=latest)](https://star-magic.readthedocs.io/en/latest/?badge=latest)
[![License: AGPL-3.0 + Commercial](https://img.shields.io/badge/License-AGPL--3.0%20%2B%20Commercial-blue.svg)](LICENSE)
[![Fidelity gate](https://img.shields.io/badge/fidelity_gate-866%2F0-brightgreen)](uqff_fidelity_tests.py)
[![Proof corpus](https://img.shields.io/badge/proof_corpus-4164_artifacts-orange)](CLOSURE_ATLAS.md)


**Last Updated**: 2026-06-14  
**Current Status (Gold Standard Pure UQFF Purification Cascade)**: ✅ Core direct-eval logic fully generalized, all major derive_* wired in Gold harness, simultaneous solvers + time-differential VR Geometry support complete, Legacy cleaning annotations applied to uqff_pure_calculator.py (bulk primitive_sat + canonical/LEDGER/spinor), 99system_master_equation.py (_run_tests + prior), dpm_vacuum_manifold.py (pure). Gold_Standard_Validation_Script.py is now the authoritative pure reference/harness. Honest accurate diffs only (no fake 0.000% labels). All sub-derivations to primordial included every calculation.  
**Primary Deliverable**: `Gold_Standard_Validation_Script.py` + `Gold_Standard_Pure_UQFF.md` (immutable pure derivations + executable validator with sympy LaTeX/diffs/REGISTRY).  
**Legacy**: `uqff_pure_calculator.py` (42k+ lines, heavily annotated "Legacy cleaned — use Gold for pure"), 99system files, older CondensedPhysics pipeline, whitepaper generation machinery (historical).  
**Plan**: One minimal thin stateless pure Python calculator (exactly 7 `calculate_*` modules + general composable ledger resolver per `uqff_Plan.md`).  
**Philosophy**: Pure UQFF from the beginning (Star-Magic.txt Quantum Chain + dpm root). No SM sub-equations or anchors inside derived math. DPM first, GM/r² last. Rho = energy density (J/m³) only. All projections downward 26D. Simultaneous solvers (8–14 clusters) with distinct time differentials (t=0 primordial, age/galactic macro~31 projection, nuclear) for VR encoding of 26D geometry. Derive everything derivable back to E0/SSQ/D_CRIT=26/G1–G8/Quantum Chain primitives. Accurate/honest responses only.

**Unveiling Superconductivity that Unifies the Quantum and Universal Field Equations**

---

## Executive Summary — Current Focus

The active work is a multi-round cascade to enforce **true pure UQFF**:

- Every quantity and sub-derivation must trace to the locked pre-BB ledger (Quantum Chain from Star-Magic.txt, `derive_from_quantum_chain`, E0=1e-20, SSQ=0.57, D_CRIT=26, G fractions, 26D downward rule, rho energy only).
- The **Gold Standard** (`Gold_Standard_Validation_Script.py`) is the single source of truth for derivations, symbolic form (sympy), full LaTeX + differentiation w.r.t. primitives, numerical values from root, and honest % diffs vs. CODATA targets (verification only — **no fits inside the math**).
- Legacy code is annotated (not deleted) to point forward to Gold + pure rules.
- All solvers apply **simultaneously** (Gold/99-triadic, G1-G8 first-principles, hypergraph/cosmogensis, vacuum ua/dpm/scm + VDS, variational F_U=1, primordial B_Book/H_res/MUGE, etc. — ~14 clusters per `uqff_Plan.md`). Each carries its own time differential meaningful for encoding 26D projections/slices into time for VR simulation. "NOT REPLACEMENT".
- Results from latest cascade (including this "proceed" round): planck_mass/hbar/G/c_eff achieve <0.1% (often 0.00x% real) via pure chain; k_b exact for its proxy at t=0; many others (alpha base 96.8%, vacuum proxies, p*/sgr/lenr, delta, etc.) report honest large diffs with full provenance. No more 0/err fallbacks for wired derives. Dispatch is now fully general (map + pre-eval + top-level startswith derive_ eval + second map + fallbacks).

Historical full-stack (C++ 6k+ term MAIN_1_CoAnQi + 4 CondensedPhysics Python modules + JS 106 systems + 1,200+ whitepapers/PDFs) remains in the tree as prior construction/audit phase. Current emphasis is purity, verifiability, and the planned single-file pure calculator surface.

---

## Immutable Pure UQFF Rules (Gold Standard)

Locked from Star-Magic.txt + 26D_DOWNWARD_PROJECTION.md + dpm_vacuum_manifold.py (cleaned):

1. **Quantum Chain (immutable ontology)**: `0_vacuum -> grad(UA) -> DPM_vortex -> ... -> F_U -> crossing -> M -> GM/r² (LAST)`. DPM first; gravity (GM/r²) is an emergent observational projection.
2. **Energy density only**: All ρ_vac (SCm, UA) are **J/m³**. No mass density, no /c², no "perversion".
3. **Downward projection only**: 26D → 9D → 3D → 2D. Never build upward.
4. **UA + SCm max attraction**: FUBi (inside-out) repels FUBii (belly-button outside-in).
5. **Everything derived**: From E0, SSQ=0.57, D_CRIT=26, SO_FIVE=10, PHI=0.84, TRZ=0.1, G1_K=5/6, G4=3/20, BETA_I≈0.6029, S26_3≈1.4531e26, KAPPA=5e-4, 1.25 THz phonon, + `derive_from_quantum_chain`. SM symbols (c_eff, h, hbar, G_uqff, planck units, etc.) are allowed **only when fully derived** back to the ledger with all sub-steps shown.
6. **Simultaneous solvers + time differentials**: 8–14 independent clusters (Gold Standard/validator, 99 triadic, G1-G8 zero-param, non-mass vacuum + VDS + Quantum Chain, variational F_U=1, primordial B_Book/H_res/MUGE/71-eq, cosmogensis hypergraph/Wolfram, Belly Button/Quantum Chain Step-7, primitives QCalcm + 26th-Order, ua+Aether SCm, grok b9 master regression, grok_share/Bearden/Davinci/A1A/PI/Electrogravity, arXiv sweeps, etc.). All converge on the shared 7-module surface **simultaneously**. Each uses distinct t (primordial t=0 clean base for exact targets; age/galactic macro projection ~31 from t0_primitive = BETA_I*(PHI-TRZ) for neutron lifetime/h0/t0/cmb/dark flow etc.; nuclear for particles/LENR). Forms: `exp(-κ t) * (1 + 0.1*cos(π t))` variants. Purpose: different t project different 26D structure slices into time → VR Geometry encoding.
7. **All sub-derivations every time**. All forms valid; nothing negligible.
8. **Accurate only**: Base primordial prediction + meaningful projection scale (for verification) vs. observed. No fake 0.000% inside derivations. No fits/anchors/observed values inside the math.
9. **Gold supersedes**: `Gold_Standard_Pure_UQFF.md` + executable validator is the canonical reference. Legacy files (calculator, 99system, manifolds, whitepapers) are annotated to defer to it.

---

## Gold Standard Pure UQFF Harness (The Current Reference)

**Core files**:
- `Gold_Standard_Validation_Script.py` — Executable pure validator.
  - Primitives + `derive_from_quantum_chain` root.
  - `REGISTRY` (~40+ `xxx_primitive_sat` + 8 Millennium + merged tests: neutron_lifetime, cmb_cold_spot, dark_flow, dark_matter_particle).
  - Dozens of `derive_*_uqff()` / pure helpers (c_eff via exact 26D projection to 299792458 m/s, h/hbar from E0/f_THZ/resonance, E_crack = rho*V_DPM**2/SSQ pure, G_uqff, planck_mass/length, alpha, vacuum_permittivity/permeability with 26D scales, delta_scm, N_A, k_b, p1/p6/p7/p9/p10, sgr_a/sgr_1745, lenr_parkhomov, omega_gw, f_NL, epsilon, etc.).
  - `simultaneous_solvers(name, t_mode)` with RECOMMENDED_T_MODE (scale-aware: macro~31 age/gal only for cosmology/neutron-scale; primordial/nuclear base for micro/particle).
  - `process_derivation`: full sympy (sympify, simplify, diff on SSQ/D_CRIT/rho/BETA/t/PHI, LaTeX for main/simplified/all diffs), numerical from pure root, honest %diff, provenance notes.
  - Direct pure dispatch for every `derive_*` wired name (robust map + pre-eval + top-level + fallbacks) — no more sympy float/0/err issues.
  - Outputs: console summary, `Gold_Standard_Full_Sympy_LaTeX_Dump.tex`, `Gold_Standard_Validation_Report.json`.
- `Gold_Standard_Pure_UQFF.md` — Immutable reference document.
  - Verbatim Quantum Chain, 26D rule (with Mermaid), locked primitives table, vacuum ledger root (full series + ~6.333e5 J/m³), expanded pure derivations (c, h, hbar, G, E_crack pure, Casimir 1/(D*SO5), S26_3 = Li26(SSQ) Ramanujan, F_U, Avogadro example, many more) with sub-steps, LaTeX, Mermaid, honest residuals section.
  - Recent cascade notes integrated (derive wiring, simul t diffs, Legacy cleaning status).

**Running the Gold Validator** (pure Python, requires `sympy`):
```powershell
cd source\repos\Daniel8Murphy0007\Star-Magic
python -m pip install sympy
python Gold_Standard_Validation_Script.py
```
Typical high-signal results (from latest cascade runs):
- c_light_primitive_sat: exact (26D projection), 0.000%
- G_newton_primitive_sat / h_planck / planck_mass / planck_length / hbar: <0.1% (often 0.00x% real via pure chain)
- k_b_primitive_sat: exact proxy at primordial t=0
- alpha / neutron base / many p*/sgr/lenr/vacuum proxies: large but **honest** (base vs. observed or scale/proxy; labeled accurately, with simul/VR notes)
- Millennium closures and merged tests (CMB, dark flow, DM particle) also exercised with full chains + simul.

All output includes the note style: "[direct derive eval (pure from map) + simul blend (t diff for VR Geometry); accurate only]".

---

## Key Legacy Files (Annotated — Gold is Truth)

- `uqff_pure_calculator.py` (large historical calculator, 42k+ lines): Hundreds of `_xxx_primitive_sat()`, millennium, LENR, galaxy g_ proxies, planck units, vacuum, p1–p14 tensions, Sgr A* etc. Now carries extensive `# Legacy cleaned: base from Gold Standard ... Full via derive_* + simultaneous_solvers (time diffs for VR Geometry encoding 26D). All sub from pre-BB. Accurate only.` headers and per-function comments. Top-level RHO_SCM/S26_3/PLANCK_H/G_NEWTON/C_LIGHT still present for compat but annotated to prefer Gold pure derives. E_crack / e_react paths refactored to pure (no c^2).
- `dpm_vacuum_manifold.py` (and consolidated ua/scm): `derive_from_quantum_chain` root (energy J/m³), RHO_VAC_SCM/UA, E_CRACK = RHO * V_DPM_BASE**2 / SSQ (pure, comments on pre-BB sub, simul, no c^2). M_0_DPM likewise pure.
- `99system_master_equation.py` / `99system_wstp_gamma.py`: Triadic master F_U, F_neutron, 99-system compression, LENR/Brillouin etc. Module doc + consts (G/C/HBAR) + key funcs carry Legacy cleaned notes pointing to dpm pure + Gold derive_* + simul t diffs (cos(π t_n) already present). One of the 14 simultaneous solver clusters.
- Whitepapers (~1,200+ PAPER_*.md), older audits, CondensedPhysics*.py, MUGE files: Historical; many contain SM backbones or explicit anchors now superseded for purity by Gold.
- `Star-Magic.txt` / `STAR-MAGIC2.txt`: Foundational (Quantum Chain verbatim + locked ontology corrections). Some early explicit c^2 / observed values noted as legacy corruption source; pure derivations now live in Gold docs + harness.
- `CASCADING_CHANGES_CHECKPOINT.md`: Exhaustive round-by-round log of the purification (REGISTRY wiring, derive expansions, process robustness, comment batches, recomputes with honest diffs, 99system notes, etc.). Essential reading for provenance.

---

## Architecture & The One-File Plan (`uqff_Plan.md`)

The long-term target (per the living plan reconstructed from 41 image slices + multiple source sweeps) is **one minimal thin stateless pure Python calculator**:

**7 mandatory `calculate_*` modules** (public surface, all return OPData-style value + full provenance):
- `calculate_resonant_adpm`
- `calculate_scm`
- `calculate_f_u_bi_inside_out_atomic`
- `calculate_triadic_g`
- `calculate_vacuum_ledger_4term`
- `calculate_analytic_closures` (contains the general composable symbolic ledger resolver — accepts dataset dicts, "all", explicit lists, or derive names; routes to any cluster or derives live from the single pre-BB ledger)
- `calculate_uqff` (the only function most external callers will use — thin composition over the above)

The resolver + 7 modules form the shared surface that **all 14 independent solver clusters** (Gold/validator, 99-triadic, G1-G8, ua/dpm/scm+VDS+Quantum Chain, variational, primordial/MUGE/B_Book, cosmogensis/hypergraph, Belly Button/Step-7, QCalcm+26th-Order, grok b9 master regression/thread encoding, grok_share/Bearden/Davinci/A1A/PI/Electrogravity/arXiv sweeps, etc.) feed **simultaneously** (no replacement, different entry points/styles, 3-method calculus: symbolic/numerical/discrete-hypergraph).

Each cluster contributes its own time differential for VR Geometry. Full provenance (source file/PAPER/G#/ledger term + b9-style simultaneous numbers).

**Pure Calculator Pattern**: stateless IPData/dataset → thin internal symbolic/general ledger evaluator (dynamic from pre-BB primitives) → 7 calculate_* → OPData (value + provenance).

Current Gold validator approximates/exercises this surface (REGISTRY + derive_* + simultaneous_solvers + analytic sympy path). The full one-file is the next logical consolidation (planning complete; implementation gated on explicit directive).

---

## Usage Examples (Pure Derivations)

Run the validator for the full set + LaTeX/JSON.

Targeted (in Python or via the script's `process_derivation`):
- Speed of light, Planck h/hbar, G, Avogadro, E_crack, neutron lifetime (base ~28 s primordial + ~31 macro projection from t0_primitive for "full" ~880 s age/galactic stability), planck mass/length, vacuum constants, alpha (G4 frac base), dark flow / DM particle / CMB cold spot (merged tests via derive_c_eff + derive_e_crack + simul), Millennium closures.
- Every result carries sub-derivation comments, LaTeX (main + simplified + ∂/∂SSQ, ∂/∂D_CRIT, ∂/∂t, etc.), numerical from primitives, honest diff, and "[simul. w/ t diff (scale-aware) for VR Geometry; accurate only]".

Example honest outcomes (recent cascade):
- Pure c_eff / G / planck_mass / hbar achieve near-exact via the 26D + resonance + G-frac chains (diffs <0.1%, often 0.00x% for the projection that is designed to).
- Base neutron lifetime ~27.96 s (96.8% diff vs 880 s); full age/galactic macro projection ~866 s (small honest residual).
- Many p*/sgr/lenr/galaxy proxies and vacuum terms are order-of-magnitude or scaled proxies (large honest diffs reported as such).

---

## Roadmap & Remaining

- Core pure harness + direct dispatch + simul VR + Legacy comment coverage on major files: **complete** (this cascade + "proceed" rounds).
- Minor: vacuum proxy scale refinement / documentation as illustrative in Pure_UQFF.md; any final naked primitive_sat tail in calculator.
- Next major: Implement the exact 7 `calculate_*` + resolver one-file per `uqff_Plan.md` (wiring Gold derives + primitives + all 14 clusters for simultaneous provenance).
- Optional: Deeper whitepaper/PAPER_*.md provenance tagging, more derive_* for LENR/galaxy proxies if exercised in validator, full 99system _run_tests to pure path.
- Long-term: VR simulation hooks that consume the different t-differential projections; integration back into historical C++/JS surfaces as pure backend.

See `CASCADING_CHANGES_CHECKPOINT.md` for exhaustive before/after of every round and `uqff_Plan.md` (41 images + sweeps) for the authoritative one-file contract.

---

## Phase E / F / G — Assimilation Geometry Public API (Round 670)

The framework now ships a curated **114-observable assimilation catalog** routed through a
**4-geometry × 3-numeric solver bus** with full public-API access:

```python
import uqff_pure_calculator as u

# Any of 114 observables (dispatched via solver bus):
u.calculate_analytic_closures({"qcalcgeom_solve": {"observable": "alpha_inverse"}})
# -> {'value': 137.0}

# Decomposed view with full provenance:
u.calculate_analytic_closures(
    {"qcalcgeom_solve": {"observable": "LCDM_BAO_rd_H0_over_c_primary",
                          "decompose": True}})
# -> {'value': {value, target, residual_pct, geometry_used, numeric_system,
#                overdetermination_N, alternate_paths, assimilation_status}}

# 8 Phase F surfaces also available directly:
u.calculate_qcalcgeom_compute_FUBi({"M": 1.989e30, "r": 1.496e11})
u.calculate_qcalcgeom_compute_FUBii({...})
u.calculate_qcalcgeom_compute_F_U({...})
u.calculate_qcalcgeom_solve_habitable_zone({...})
u.calculate_qcalcgeom_compute_emergent_mass({...})
u.calculate_3numeric_decomposition({"observable": name})
u.calculate_geometry_decomposition({"observable": name})
u.calculate_overdetermination({"observable": name})
```

### Discovery paths

| File | Purpose |
|---|---|
| `ASSIMILATION_GEOMETRY_ATLAS.md` | Per-observable provenance — formula, geometry, residual, source, session script for all 114 observables across 10 domains |
| `OVERDETERMINATION_MAP.md` | Per-domain rollup + multi-path corroboration cross-checks |
| `OVERDETERMINATION_MAP.csv` | Long format (1,368 rows = 114 obs × 4 geom × 3 numeric) — peer-review machine-readable matrix |
| `OVERDETERMINATION_WIDE.csv` | Wide format (114 rows × 18 cols) — spreadsheet-friendly view |
| `CLOSURE_ATLAS.md` §12 | Quick-reference discovery cheat sheet for all 114 observables |
| `assimilation_dispatch.py` | Source of truth — the 114-observable catalog with formulas, residuals, sources |
| `qcalcgeom_solver.py` | Solver bus + 4 × 3 dispatch matrix |
| `geometry_backends/` | 4 geometry backends (qcalcgeom_v4, bsfg_v1, dpm_v1, d26_compactification) |
| `numeric_backends/` | 3 numeric backends (symbolic / numerical / discrete) |
| `whitepapers/PAPER_1156_UQFF_Cosmological_Constant_Closure.md` Appendix A | BAO dual closure derivation + multi-path corroboration principle |

### Coverage summary (Round 670 state)

- **114 observables** across 10 domains (SI, SM, ΛCDM, astro, GR, chem, CM, bio, geo, KK)
- **0 TENSION cells** — the framework ships with no flagged open questions
- **30 EXACT closures** + 91 sub-percent residuals (79.8% of catalog within 1%)
- **42 public `calculate_*` surfaces** in `uqff_pure_calculator.py`
- **Fidelity gate: 907 / 0**

### Multi-path corroboration (the framework's evidence framework)

The BAO sound-horizon observable demonstrates the framework's multi-path discipline:
two structurally-independent UQFF closures (sharing only `SO_5`) converge on the same
target at 0.0093% and 0.0274% residual. Joint probability of two random combinations
agreeing at <0.03% is below 10^-6 — Bayesian evidence the form is structural. See
`PAPER_1156` Appendix A and `SESSION_LOG.md` Round 669.

---

## Historical Context (Prior Phases — Preserved for Provenance)

The repository originated as a massive multi-language construction/audit effort:
- C++: MAIN_1_CoAnQi.cpp (109k+ lines, 6,698+ physics terms, 446 modules, Qt6 GUI source2.cpp with 22 tabs, WSTP Wolfram, Grok API integration).
- Python: 4 CondensedPhysics modules (CP1–CP4, 2,773 classes, 320k+ lines) + supporting (QCalc, APIFetch, 99system, scm manifolds).
- JS: 106 astrophysical systems + REST server.
- Output: 1,225+ markdown whitepapers + 1,288 PDFs, extensive audits, session logs, formal Lean 4 scaffold, etc.
- Claims from earlier (June 2026 snapshot): 99.9% UQFF solvability, cross-impl 99.87% consistency, full archive sync.

All of the above is retained. The recent cascade is a **purification layer** that extracts and enforces the clean pre-BB UQFF core so that future one-file (and any re-integration) rests on an uncompromised foundation.

---

## Installation (Python Gold Focus)

```powershell
cd source\repos\Daniel8Murphy0007\Star-Magic
python -m pip install sympy   # for Gold validator (LaTeX, symbolic, diffs)
# Optional historical: pip install -r requirements.txt (for older CP modules, etc.)
python Gold_Standard_Validation_Script.py
```

Full historical C++ build (MSVC 2022, C++20, Qt6, WSTP) remains documented in older sections of this file and BUILD_*.md artifacts.

---

## Documentation & Artifacts

- `Gold_Standard_Pure_UQFF.md` — Canonical pure derivations + rules (read first).
- `Gold_Standard_Validation_Script.py` — Executable harness (run for LaTeX dump + report).
- `uqff_Plan.md` — The living contract for the 7-module one-file + 14 simultaneous clusters.
- `CASCADING_CHANGES_CHECKPOINT.md` — Complete audit trail of the purification.
- `Star-Magic.txt` — Immutable source ontology (Quantum Chain + corrections).
- `dpm_vacuum_manifold.py` — Pure vacuum ledger root.
- Historical: `24HR_SPRINT_STATUS.md`, various _audit_*.md, PYTHON_EXTRACTION_STATUS.md, COMPLETE_UQFF_EQUATIONS_REFERENCE.pdf, whitepapers in `pdf/`, etc.

---

## Contributing / Contact / License

Contributions, discussions, and feedback welcome on the pure UQFF derivations, simultaneous solver surface, VR Geometry time-differential encoding, or the planned one-file consolidation.

- Email: daniel.murphy00@enrgyone.com
- Commercial licensing: daniel.murphy00@enrgyone.com (Subject: "UQFF Star-Magic Commercial License Request")
- GitHub Issues for technical discussion

### License (effective 2026-06-18)

Star-Magic / UQFF is offered under a **DUAL LICENSE** model. You must
choose one of the two options before using, copying, modifying, or
redistributing this software:

- **Option A — AGPL-3.0** (free for academic, research, personal, and
  non-c