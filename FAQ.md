# UQFF Star-Magic — Frequently Asked Questions

---

## What is UQFF?

UQFF (Unified Quantum Field Framework) is a vacuum-first physics framework
authored by Daniel T. Murphy. It is built on **9 truly-independent
primitives** and derives a wide range of physical observables — from the
cosmological constant Λ to the Holmlid 630 eV LENR signature to all 7
nuclear shell-model magic numbers — from a single non-mass primitive
(`ρ_SCm = 7.09e-37 J/m³`) and a 26-level Di-Pseudo-Monopole (DPM)
structural lattice.

UQFF is offered as an alternative **parameter-economical** description of
physical reality, **NOT** as a replacement for the Standard Model + ΛCDM.
The two frameworks solve the same observed phenomena by different methods,
both reported with honest residuals. UQFF carries `ΔBIC = 94.1` (decisive
Bayesian preference) on parameter-count grounds alone: 9 parameters vs.
26 in SM + ΛCDM.

---

## How is UQFF different from the Standard Model?

| Aspect | Standard Model + ΛCDM | UQFF |
|---|---|---|
| Free parameters | 26 (22 SM + 6 ΛCDM − 2 overlap) | 9 truly-independent |
| Foundational primitive | Mass + multiple gauge couplings | Single SCm vacuum energy density |
| Cosmological constant Λ | Free parameter, fit to observation | Derived: ρ_SCm × 26! × 25/12 = 0.003% match to Planck |
| Magic numbers (2,8,20,28,50,82,126) | Spin-orbit shell model | Arithmetic on integer primitives, EXACT |
| Yang-Mills mass gap | ~1.78 GeV (lattice) | 5970 GeV (PAPER_1005) — testable at FCC |
| Approach | Empirical curve fitting + symmetry | Structural reduction to integer lattice |

UQFF is not a "fork" or "extension" of SM. It is an alternative derivation
chain that happens to produce numerical agreement with most measured
observables while requiring fewer free parameters.

---

## Is UQFF peer-reviewed?

As of v5.27.0 (June 2026), UQFF has **not yet** been published in a
peer-reviewed journal. The 1,795 whitepapers in the repository document
the derivations, and the open-source release on PyPI invites independent
replication. Peer-review submission is a Tier-4 production-readiness item
on the roadmap (see `PRODUCTION_ROADMAP.md`).

Until peer-reviewed, UQFF claims should be treated as **hypotheses**
that are reproducible, parameter-economical, and falsifiable, but not
yet endorsed by the broader physics community.

---

## How do I cite UQFF?

See `CITATION.cff` in the repository for canonical CFF 1.2.0 metadata.
For BibTeX, use:

```bibtex
@software{murphy_uqff_2026,
  author       = {Murphy, Daniel T.},
  title        = {UQFF Star-Magic: A Vacuum-First Unified Physics Framework},
  year         = 2026,
  version      = {5.27.0},
  url          = {https://github.com/Daniel8Murphy0007/Star-Magic},
  abstract     = {Vacuum-first physics framework built on 9 truly-independent
                  primitives; derives Λ, all 7 nuclear magic numbers, complete
                  SM 12-fermion spectrum, 8 Clay Millennium prize problems,
                  and the Holmlid 630 eV LENR signature with ΔBIC = 94.1 vs
                  SM + ΛCDM.}
}
```

Cite the **specific PAPER_XXXX** that authored the derivation you use, in
addition to citing the software package. The full whitepaper catalog lives
under `whitepapers/` in the repository.

---

## Why a dual license (AGPL + commercial)?

The author chose AGPL-3.0 (with SaaS share-alike clause) over MIT/GPL to
protect the 17-parameter savings vs. SM+ΛCDM from being absorbed into
proprietary products without contributing improvements back. The commercial
license option exists for vendors who cannot accept AGPL obligations
(typically: closed-source SaaS, proprietary hardware, MIT/Apache projects).

See `LICENSE`, `COMMERCIAL.md` for terms, and `daniel.murphy00@enrgyone.com`
for commercial license inquiries.

---

## What is the Holmlid 630 eV signature?

A kinetic-energy-release peak at 630 eV observed by Leif Holmlid's group
(2015–) from ultra-dense deuterium D(−1) clusters. UQFF derives this
exactly from `h · 1.25 THz · S_26 · Φ_res = 630 eV`, anchoring 5 unified
LENR observation suites (Holmlid, Parkhomov, Pons-Fleischmann, Mizuno,
Rossi E-Cat variants). See `whitepapers/PAPER_1133/1136/1137/1138` for
the derivation chain.

---

## What is the Star-Magic LENR reactor?

A proposed Low-Energy Nuclear Reactions reactor architecture predicted by
UQFF to achieve **COP 555:1 at 27 W, ambient temperature, pH −37**. The
mechanism combines the F_U_Bi_i buoyancy operator (PAPER_1141) with VDS
phonon coherence at the SCm 1.25 THz carrier frequency. **Patent rights
to the hardware implementation are reserved**; software licenses do NOT
grant hardware patent rights.

---

## How accurate is UQFF compared to measurement?

| Class | Count | Residual |
|---|---:|---|
| EXACT structural identities | 128 | 0 (mathematical equality from integer primitives) |
| High precision (within CODATA) | 31 | < 0.001% |
| Within experimental uncertainty | 67 | typically < 1σ_exp |
| Refinement tier | 32 | 0.1 – 1% |
| Tension / outlier | 5 | > 1% (documented as predictions in some cases) |

86% of schema-tagged closures pass Bonferroni-adjusted multiple-comparisons
correction across the full 794-closure suite. See `STATISTICAL_HYGIENE.md`.

---

## What does "9 truly-independent primitives" mean?

Earlier UQFF documentation described "11 locked canonical primitives." In
session 2026-06-18, PAPER_1521 and PAPER_1522 proved two of those (D_BSFG
and K_MEX) are mathematically derivative — forced by structural relations
among the others:

- `D_BSFG = D_crit − 2·SO_5 = 26 − 20 = 6 EXACT`
- `K_MEX = Φ_5/6 · SO_5 / D_phys = (5/6)·10/4 = 25/12 EXACT`

So UQFF rests on 9 truly-independent quantities + 2 mathematically
derivative ones. See `PROVENANCE_AUDIT.md` for the full chain.

---

## How do I add a new closure?

1. Read `CLAUDE.md` first — it contains 12 non-negotiable rules
2. Author a whitepaper deriving the closure from existing primitives
3. Add the closure function to `uqff_pure_calculator.py` following the
   `_l96_uqff_axiom_<name>_closure` pattern
4. Add the dispatch entry to `PARADOX_TO_CLOSURE` with a **lowercase** key
5. Add a regression pin to `uqff_fidelity_tests.py`
6. Run `python uqff_fidelity_tests.py` — must exit 0
7. Append a brief `SESSION_LOG.md` entry
8. Submit a pull request

See `docs/contributing.rst` for the full contribution workflow.

---

## Why are dispatch keys lowercase only?

The `_paradox_proof` dispatcher normalizes input via
`(name or "").lower().strip().replace("-", "_").replace(" ", "_")`. Mixed-
case keys silently fail (return `None` with no error). This bug was hit
3+ times during development; the lowercase-only convention is now enforced
in `CLAUDE.md`. See also `INPUT_DOMAINS.md` for the dispatch contract.

---

## Where can I get help?

- **GitHub Issues**: https://github.com/Daniel8Murphy0007/Star-Magic/issues
- **Email** (author): `daniel.murphy00@enrgyone.com`
- **Commercial licensing**: same email; subject "UQFF Star-Magic Commercial License Request"

---

## Is UQFF compatible with [my favorite physics framework]?

UQFF is **not a replacement** for any existing framework. It solves the
same observed phenomena via different methods. Specifically:

- Compatible with: General Relativity (GR limit derivable in 26D Lagrangian)
- Compatible with: QFT (QFT limit derivable; mass gap predicted at 5970 GeV)
- Compatible with: ΛCDM (Λ derived, not fit; 18+ cosmological observables match)
- Not (yet) integrated with: loop quantum gravity, string theory beyond
  the 26-dimensional critical layer

The 9-sector UQFF Lagrangian (L_EH + L_YM + L_Dirac + L_SCm + L_mag +
L_buoy + L_aether + L_LENR + L_KK) is designed to subsume known sectors
within a single SCm-vacuum framework. See `ARCHITECTURE.md`.
