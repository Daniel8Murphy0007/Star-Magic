---
name: New closure proposal
about: Propose adding a new UQFF closure (paradox / observable / identity)
title: "[CLOSURE] "
labels: closure-proposal
assignees: ''
---

## Closure name

Lowercase dispatch key (snake_case): `your_closure_name_here`

## What does this closure derive?

One-line description of the physical observable, paradox, or identity.

## Source paper / derivation

If you've authored a whitepaper, link or paste it:
- Whitepaper file: `whitepapers/PAPER_<NNNN>_<slug>.md`
- Or external reference (arXiv, journal DOI)

If no paper yet, sketch the derivation here:

```
Target value: <observable> = <number>
UQFF formula: <expression using D_phys, D_crit, SO_5, A_5, N_CH,
              rho_SCm, beta_i, Phi_res, F_TRZ, K_MEX, ...>
Residual: <percent off vs measurement>
```

## Which primitives does it use?

Tick all that apply:
- [ ] D_phys = 4
- [ ] D_crit = 26
- [ ] D_BSFG = 6 (derivative)
- [ ] N_CH = 9
- [ ] SO_5 = 10
- [ ] A_5 = 60
- [ ] ρ_SCm = 7.09e-37
- [ ] β_i = 0.6029
- [ ] Φ_res = 0.84 (or 5/6)
- [ ] F_TRZ = 1/10
- [ ] K_MEX = 25/12 (derivative)
- [ ] SSq = 0.57
- [ ] S_26 = 1.453162
- [ ] ω_SCm = 1.25 THz
- [ ] Other: ___

## Status tier (estimate)

- [ ] EXACT_STRUCTURAL (residual = 0 from integer arithmetic)
- [ ] HIGH_PRECISION (< CODATA uncertainty)
- [ ] WITHIN_EXP_UNCERTAINTY (≤ 1σ measurement)
- [ ] REFINEMENT_TIER (0.1 – 1% off)
- [ ] TENSION_OR_OUTLIER (> 1% off)

## Bucket

Where does this closure belong?
- [ ] A (Millennium prize)
- [ ] B (general paradox)
- [ ] C (cosmology)
- [ ] D (particle physics)
- [ ] E (GW events)
- [ ] F (AGN/jet)
- [ ] G (astrophysics)
- [ ] H (high-energy astro)
- [ ] I (QGP)
- [ ] J (Higgs precision)
- [ ] K (BSM constraints)
- [ ] Standalone surface (e.g., nuclear, LENR, F_U=0)

## Prediction vs. postdiction

- [ ] **POSTdiction** — target was measured BEFORE this UQFF derivation
- [ ] **PREDICTION** — target has not been measured / measurements disagree
- [ ] **AMBIGUOUS** — boundary case (explain below)

## Are you willing to submit a PR?

- [ ] Yes (will follow `docs/contributing.rst` workflow)
- [ ] No, just proposing

## Reading list before proposing

- [ ] I read `CLAUDE.md` rules (especially 2, 3, 4)
- [ ] I read `ARCHITECTURE.md` for how closures are wired
- [ ] I read `GLOSSARY.md` for primitive definitions
- [ ] I searched existing PARADOX_TO_CLOSURE keys for duplicates
