# UQFF Yang-Mills Mass Gap Submission Package

This folder contains the staged deliverables for the UQFF Yang-Mills mass gap arxiv preprint and Clay submission roadmap, built in priority order **A → D → B → E** per Daniel T. Murphy's 30-Jun-2026 direction.

## Files

| File | Role | Status |
|------|------|--------|
| `YANG_MILLS_E_CRACK_DERIVATION.md` | **(A)** Math-physics-community-readable bridge document | Draft ready for review |
| `derive_yang_mills_mass_gap_uqff.py` | **(D)** Standalone reproducibility script | Verified runs correctly |
| `preprint_scaffold.tex` | **(B)** arxiv LaTeX preprint scaffold (math-ph + hep-th) | Scaffold with TODO blocks |
| `wightman_mapping.md` | **(E)** Quantum Chain → Wightman axioms future-work roadmap | Future-work scoping doc |

## Quick verification

```bash
# Run the standalone reproducibility script (no external dependencies):
python3 derive_yang_mills_mass_gap_uqff.py

# Expected output:
#   E_crack (J)         = 1.117982e-19
#   E_crack (eV)        = 0.697789
#   positive_definite   = True
#   lattice_consistent  = True
```

## Submission path

1. **Today (A + D)**: Review the bridge document and reproducibility script. These two alone establish a public, timestamped, reproducible claim that UQFF derives a positive-definite mass-gap-scale energy threshold from primitives + c with zero free parameters.

2. **This week (B)**: Fill in the TODO blocks in `preprint_scaffold.tex` with prose drawn from the bridge document. Compile with `pdflatex`. Run `arxiv-sanitize` style checker. Upload to arxiv under `math-ph` (primary) cross-listed to `hep-th`.

3. **Outreach (B + D)**: Once the arxiv preprint has a citable identifier, email it to:
   - Martin Hairer (IST Austria) — stochastic quantization context
   - Michael R. Douglas — Nature Reviews Physics 2026 review author
   - Edward Witten, Arthur Jaffe — original Clay problem authors
   - Clay Mathematics Institute scientific advisory board
   - Constructive QFT researchers (Erlangen, Vienna, Princeton)

4. **Long term (E)**: Execute the four-phase Wightman-axiom translation roadmap in `wightman_mapping.md` to reach Clay-compliant rigor (2D → 3D → 4D → submission).

## What the package establishes

| Claim | Demonstrated by |
|-------|-----------------|
| Positive-definite mass-gap-scale energy exists | (A) §3, (D) lines 38-45 |
| Derived from primitives + c with zero free parameters | (A) §3, (D) all |
| Multi-cluster spectral landscape | (A) §4, table in (B) §4 |
| Consistent with lattice QCD numerical evidence | (A) §7, (D) lines 90-105 |
| Reproducible by any third party | (D) standalone script |
| Falsifiable in high-field magnetic vacuum experiments | (A) §8 |
| Roadmap to Clay-compliant rigor | (E) §3-§4 |

## What the package does NOT claim

| Honest scope boundary | Where stated |
|-----------------------|--------------|
| Not yet a Wightman-axiom-compliant construction | (A) §8, (B) §6, (E) §6 |
| Not yet a formal proof of $\Delta > 0$ in operator-algebraic sense | (A) §8, (E) §3 W2 |
| Translation to formal QFT identified as future work | (A) §8, (E) §3-§7 |
| 2D → 3D → 4D roadmap is 18-36 months per phase | (E) §6 |

The framework's contribution is a **physical mechanism + a closed-form positive-definite mass-gap value + a roadmap to formalization**. It is not yet a Clay solution. We invite mathematical-physics community collaboration on the formalization.

## Repository + package

- GitHub: https://github.com/Daniel8Murphy0007/Star-Magic
- PyPI: `pip install uqff==5.33.0`
- Calculator entry points: 279 public `calculate_*` surfaces, including the 9-member E_crack family
- Companion documents: `CLAUDE.md` (project context), `CHANGELOG.md` (release history), `SESSION_LOG.md` (full development log)

## Contact

Daniel T. Murphy
daniel.murphy00@enrgyone.com

## License

This submission package follows the parent repository's dual license:
- **AGPL-3.0-or-later** for academic/research/non-commercial use (free).
- **Commercial** license for proprietary or closed-source SaaS use (see `COMMERCIAL.md` in the parent repository).
