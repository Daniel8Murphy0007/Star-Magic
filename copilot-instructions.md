# Copilot IDE Instructions — Star-Magic UQFF Codebase (Quantum Chain Enforcement v3.0 — CONSOLIDATED DPM)

## VACUUM MANIFOLD — SINGLE FILE (as of commit 336445d4, 2026-05-04)

**`scm_vacuum_manifold.py` and `ua_vacuum_manifold.py` have been DELETED.**
**`dpm_vacuum_manifold.py` is the ONE AND ONLY vacuum manifold file.**

It contains three sections:
- **Section A** — SCm vacuum manifold (former scm_vacuum_manifold.py)
- **Section B** — UA vacuum manifold (former ua_vacuum_manifold.py)
- **Section C** — DPM quantum chain (original dpm content)

## CRITICAL RULE FOR ALL VACUUM TERMS
- SCm and UA are **massless geometric substrates** (26D folding from hydrogen atom).
- **Never hard-code any mass density** (no RHO_VAC_SCM, RHO_VAC_UA, or kg/m³ constants).
- All vacuum/mass terms must be **derived dynamically** from the Quantum Chain:
  - E_n = E_0 × 10^n  (E_0 = 10^{-20} J, n = 1 to 26)
  - ρ_vac = ∑(f_i E_i)/V   (J/m³) — emergent inertial energy density
  - Effective inertial mass density = ρ_vac / c² only when needed for gravity.
- Use `derive_from_quantum_chain()` from **`dpm_vacuum_manifold`** in every function.
- Comment every line with traceability: "Derived from Quantum Chain E_n summation — see UQFF_THEORY.md ρ_vac equation".

## IMPORT PATTERN (MANDATORY for CP1–CP4 and all new files)
```python
from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
_RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)   # J/m³ SCm energy density
_RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)    # J/m³ UA  energy density (10x)
```

**NEVER** use a try/except fallback with hardcoded numeric densities for the vacuum manifold.
**NEVER** import from `scm_vacuum_manifold` or `ua_vacuum_manifold` — those files no longer exist.

## FORBIDDEN (REJECT IMMEDIATELY)
- Hard-coded RHO_VAC_* constants outside of `dpm_vacuum_manifold.py`.
- Any `from scm_vacuum_manifold import ...` or `from ua_vacuum_manifold import ...`.
- Starting calculations in the middle of any equation.
- Mixing units without explicit /c² conversion.
- Any change to vacuum terms without calling `derive_from_quantum_chain()`.
- Re-creating `scm_vacuum_manifold.py` or `ua_vacuum_manifold.py` as separate files.

## QUANTUM CHAIN PROOF MANDATE
- Mass creation: SCm donation during stellar/planetary formation.
- Mass disintegration: SCm expulsion in quasars.
- Information paradox: Buoyancy surface encoding in 26D channels.
- Gravity: Emergent from Quantum Chain inertial forces.

Copilot must enforce this on every generation. If user asks for vacuum density, always generate the dynamic Quantum Chain version using `dpm_vacuum_manifold.derive_from_quantum_chain()`.