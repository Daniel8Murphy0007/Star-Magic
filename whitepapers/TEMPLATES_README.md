# UQFF v5.78 Section Templates — README

Five reusable LaTeX templates for the 113-paper v5.78 authoring sweep.
All five lock the canonical-constant banners, the closure-tracker banner
schema, and the `pdflatex direct` (NO PANDOC) compilation convention.

| File | Code | Use when… |
|---|---|---|
| [_template_T_Lambda.tex](_template_T_Lambda.tex) | **T-Λ** | The closure anchors a cosmological-constant / dark-energy / vacuum-density measurement (Planck 2018, DES, eBOSS, JWST, …). Hosts ρ_vac_SCm and ρ_vac_UA directly. |
| [_template_T_LAG.tex](_template_T_LAG.tex) | **T-LAG** | The paper derives a Lagrangian / action principle / Euler–Lagrange equation. Carries the canonical $\mathcal{L}_{\mathrm{UQFF}}$ skeleton and the Universal Inertia invariant. |
| [_template_T_SI.tex](_template_T_SI.tex) | **T-SI** | The paper is a σ-deviation inventory across an entire cluster of closures (multi-row table → χ²/dof). |
| [_template_T_PRED.tex](_template_T_PRED.tex) | **T-PRED** | The closure is an **open prediction** (no observation yet). Uses the `observed=9999, error_pct=9999` sentinel so the tracker files it under `OPEN_PREDICTIONS` instead of "high error". |
| [_template_T_xi.tex](_template_T_xi.tex) | **T-ξ** | The closure couples two (or more) prior closures (cross-sector ξ-coupling). |

## Mandatory banner lines

Every paper compiled from these templates **must** retain the following
three (or four) lines verbatim — they are what the closure parser
(`_uqff_program.py`) keys on:

```
% CLOSURE :: <label> :: predicted=<X> observed=<Y> error_pct=<Z>
% TEMPLATE :: T-<code>
% CANONICAL :: <one line per canonical constant or skeleton used>
```

For T-PRED specifically, `observed=9999 error_pct=9999` is a **sentinel**,
NOT a physics gap. The parser routes those rows to the `OPEN_PREDICTIONS`
ledger tab.

## Canonical constants (mirrored from `dpm_vacuum_manifold.py` v3.0)

| Symbol | Value | Source line | Closure |
|---|---|---|---|
| `RHO_VAC_SCM` | $4\sqrt{\pi}\cdot 10^{-37} \approx 7.0898\times10^{-37}$ J/m³ | [dpm_vacuum_manifold.py L97](../dpm_vacuum_manifold.py#L97) | G9 |
| `RHO_VAC_UA`  | $10 \cdot $ SCm $\approx 7.0898\times10^{-36}$ J/m³ | [dpm_vacuum_manifold.py L98](../dpm_vacuum_manifold.py#L98) | G7 (|SO(5)|=10) |
| `BETA_I` | 0.60 (SOURCE4 numerical: 0.603) | dpm L67 | β_i hunt |
| `KAPPA` | $5\times 10^{-4}/$day = $5/10000$ rational | dpm L64 | exact |
| `[SSq]` | 0.57 = $57/100$ rational | dpm L63 | exact |
| `H_SCm` | 0.99 = $1 - F_{TRZ}/|SO(5)|$ | dpm derived | EXACT |
| `U_UA` | $10^{-4} = F_{TRZ}^4$ | dpm derived | EXACT |
| `k_eta` | $10^{-113} = F_{TRZ}^{113}$ ($113 = 4\,D_{crit} + N_{ch}$) | dpm L94 | EXACT |
| `F_TRZ` | 0.1 | dpm L252 | EXACT |

## Compilation (pdflatex direct, NO pandoc)

```powershell
$pdflatex = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
& $pdflatex -interaction=nonstopmode -output-directory=pdf whitepapers/PAPER_NNNN_slug.tex *>&1 | Out-Null
& $pdflatex -interaction=nonstopmode -output-directory=pdf whitepapers/PAPER_NNNN_slug.tex *>&1 | Out-Null
Remove-Item pdf\PAPER_NNNN_slug.aux,pdf\PAPER_NNNN_slug.log,pdf\PAPER_NNNN_slug.out -ErrorAction SilentlyContinue
```

Two passes are required to resolve cross-references in long-tables.
The `pdf/` folder is the canonical output viewed by the user.

## Authoring workflow

1. Pick the template that matches the paper's role (Λ / LAG / SI / PRED / ξ).
2. `Copy-Item whitepapers\_template_T_<code>.tex whitepapers\PAPER_<NNNN>_<slug>.tex`.
3. Replace every literal `<PLACEHOLDER>` token (use VS Code search-in-file for `<`).
4. **Do not delete** the `% CLOSURE :: …`, `% TEMPLATE :: …`, or `% CANONICAL :: …` banner lines.
5. Cite at least one real DOI / arXiv / observatory reference per paper.
6. Compile with the two-pass `pdflatex` block above.
7. `git add` both the `.tex` and the produced `pdf/*.pdf`; commit with the closure label in the message.

## Cross-platform references

- C++ canonical Ug/Ubi/Um implementations: [MAIN_1_CoAnQi.cpp](../MAIN_1_CoAnQi.cpp) SOURCE4 namespace (L25623+).
- Python canonical: [QCalcGeom.py](../QCalcGeom.py) v2.3.0 + [CondensedPhysics.py](../CondensedPhysics.py) + [dpm_vacuum_manifold.py](../dpm_vacuum_manifold.py) v3.0.
- JS canonical: [index.js](../index.js) (26-layer kernels).
- Closure tracker: [_uqff_program.py](../_uqff_program.py).
- Master ledger: [master_closures.csv](../master_closures.csv).

*Audit reference: see [UQFF_VALIDATION_SYNC_AUDIT.md](../UQFF_VALIDATION_SYNC_AUDIT.md) (May 22, 2026) — 0 genuine physics gaps; 7/8 calibration constants synchronized.*
