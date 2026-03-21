#!/usr/bin/env python3
"""
upgrade_early_whitepapers.py

Batch upgrade PAPER_001–099 in whitepapers/ to include:
  • The complete F_U master equation reference (now including 4th λᵢ dissipation term)
  • k-constants: k₁=1.5, k₂=1.2, k₃=1.8, k₄=2.0
  • 4-mode UQFF operational taxonomy
  • PAPER_050  → add 26D implementation acknowledgment section
  • PAPER_200  → add operational status table to Um taxonomy

Run from the Star-Magic root:
    python upgrade_early_whitepapers.py
"""

import os
import re
import sys

WHITEPAPERS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "whitepapers")

# ──────────────────────────────────────────────────────────────────────────────
# Appendix block injected into PAPER_001-099 that lack it
# ──────────────────────────────────────────────────────────────────────────────
UQFF_FRAMEWORK_APPENDIX = """
---

## Appendix: UQFF Production Framework Reference (v4.75+)

> *Added by upgrade_early_whitepapers.py (v4.75). This appendix cross-references
> the production physics constants and master equations to enable reproducibility
> against the current codebase state.*

### A.1 Calibration Constants

| Symbol | Value | Description |
|--------|-------|-------------|
| κ | 5.0 × 10⁻⁴ day⁻¹ | UQFF exponential decay rate |
| [SSq] | 0.57 | Universal Quantized Factor |
| β_i | 0.60–0.61 | Buoyancy coupling coefficient |
| k₁ | 1.5 | Ug1 DPM-dipole coupling |
| k₂ | 1.2 | Ug2 outer-bubble charge coupling |
| k₃ | 1.8 | Ug3 string-rotation coupling |
| k₄ | 2.0 | Ug4 vacuum-concentration coupling |
| η | 10⁻²² | Inertia tensor scale |
| E_react(0) | 10⁴⁶ J | Reference reactive energy |

### A.2 F_U Master Equation (Complete — 4 terms)

$$F_U = U_{g1} + U_{g2} + U_{g3} + U_{g4} + U_{bi} + U_m - \sum_{i=1}^{4}\bigl[\lambda_i \cdot U_i(r,t) \cdot E_{\mathrm{react}}\bigr]$$

| Term | Description | Implementation |
|------|-------------|----------------|
| Ug1 | DPM magnetic dipole | `compute_Ug1_SOURCE4` / `compute_Ug1()` |
| Ug2 | Outer-field bubble (charge-reactivity) | `compute_Ug2_SOURCE4` / `compute_Ug2()` |
| Ug3 | Magnetic string rotation | `compute_Ug3_SOURCE4` / `compute_Ug3()` |
| Ug4 | Vacuum concentration (star-BH) | `compute_Ug4_SOURCE4` / `compute_Ug4()` |
| Ubi | Buoyancy force | `compute_Ubi_SOURCE4` / `compute_Ubi()` |
| Um | Universal Magnetism (Heaviside-amplified) | `compute_Um_SOURCE4` / `compute_Um()` |
| −Σλᵢ·Uᵢ·E_react | 4th dissipation term (PAPER_420) | `compute_FU_SOURCE4` / full pipeline |

**4th dissipation term parameters (PAPER_420):**  
λ₁=10⁻¹⁰, λ₂=10⁻¹², λ₃=10⁻¹¹, λ₄=10⁻¹³ (free parameters, not yet empirically calibrated)

### A.3 Um Heaviside Phase-Transition Amplifier (PAPER_421)

$$U_m^{\mathrm{full}} = U_m^{\mathrm{base}} \times \bigl(1 + 10^{13}\,\Theta(\rho_{SCm} - \rho_c)\bigr) \times \bigl(1 + A_q\cos(\Delta\omega\,t)\bigr)$$

| Symbol | Value | Description |
|--------|-------|-------------|
| ρ_c | 10¹⁵ kg/m³ | SCm critical superconducting density |
| A_q | 0.1 | Quasi-periodic beating amplitude (10%) |
| Δω | 2π/(434·365.25) rad/day | 434-year Gleisberg supercycle |

### A.4 UQFF Four Operational Modes

| Mode | Dominant Term | Primary Use Case |
|------|--------------|-----------------|
| **Compressed** | Ug_sum + Newtonian base | Isolated stellar/BH systems |
| **Resonant** | 5 resonance frequencies (aDPM, aTHz, …) | Multi-scale field interactions |
| **Buoyant** | β_i × Ubi | Expanding nebulae, stellar winds |
| **Superconductive** | Um × (1+10¹³·f_H) | Magnetars, SCm critical-density regime |

*Implementation status: all 4 modes operational in `MAIN_1_CoAnQi.cpp`, `CondensedPhysics.py`, and `CondensedPhysics2.py`.*
"""

# Markers to detect if appendix already present
APPENDIX_MARKERS = [
    "UQFF Production Framework Reference",
    "Appendix: UQFF Production Framework",
    "4-mode UQFF operational taxonomy",
]

# ──────────────────────────────────────────────────────────────────────────────
# PAPER_050 addition: 26D implementation acknowledgment
# ──────────────────────────────────────────────────────────────────────────────
PAPER_050_INSERT_AFTER = "## Abstract"

PAPER_050_NOTE = """
---

> **Implementation Note (v4.75):** This paper presents the full 26-dimensional UQFF
> manifold compactification framework. Current production code (`MAIN_1_CoAnQi.cpp`
> SOURCE115, `CondensedPhysics.py`) operationalizes **4 of 26 dimensions**: 3 spatial
> + 1 temporal (standard spacetime). The remaining 22 compact dimensions are
> analytically described in the 26D polynomial master equations (SOURCE115 § 19-system
> framework) and provide convergent correction terms at scales ≲ 10⁻³⁵ m. Full 26D
> operationalization is a planned future milestone. Results in this paper that
> reference 26D quantities are analytically correct; the numerical evaluations use the
> 4D projection unless explicitly noted.

"""

# ──────────────────────────────────────────────────────────────────────────────
# PAPER_200 addition: operational status table
# ──────────────────────────────────────────────────────────────────────────────
PAPER_200_OPERATIONAL_STATUS = """
---

## 13. Operational Implementation Status

The following table classifies each Um variant by its production status in the
current codebase (`MAIN_1_CoAnQi.cpp`, `CondensedPhysics.py`, `CondensedPhysics2.py`).

| Status | Meaning |
|--------|---------|
| ✅ **Operational** | Fully implemented, returns numerical result for any input dataset |
| 🔄 **Partial** | Core formula implemented; specialized sub-terms use placeholder values |
| 📋 **Reference** | Equation documented for reproducibility; not yet callable from pipeline |

### Core Um Variants

| Symbol | Description | Status | Notes |
|--------|-------------|--------|-------|
| Um (base) | `Σ_j[μ_j/r_j·(1−e^{−γt·cos(πt_n)})·φ̂_j]·P_SCm·E_react` | ✅ **Operational** | `compute_Um_SOURCE4`, `compute_Um()` |
| Um (Heaviside-amplified) | Um_base × (1+10¹³·Θ(ρ_SCm−ρ_c)) × (1+A_q·cos(Δω·t)) | ✅ **Operational** | PAPER_421 — integrated v4.75 |
| Um,BZ | Blandford-Znajek power extraction | ✅ **Operational** | `CondensedPhysics2.py` BZ class |
| Um,haw | Hawking temperature surface gravity | ✅ **Operational** | `bh_thermodynamics_module.py` |
| Um,qnm | QNM ringdown frequency | ✅ **Operational** | `CondensedPhysics.py` QNM class |

### Reference-Only Variants (50+)

All remaining Um,X variants catalogued in §2–§10 of this paper are **📋 Reference**
status. They contain correct analytical expressions derived from BB_C_Equations_04Sept2025.pdf
but are not yet wired into the main computation pipeline.

**Operational upgrade path:** To promote any variant from Reference → Operational,
implement a calculator class in `CondensedPhysics.py` following the `compute_Um()`
interface pattern, then register it in the `PhysicsTermRegistry`.

> *Operational status assessed: v4.75 (January 28, 2026 integration).  
> Um Heaviside amplifier + quasi-periodic modifier now operational (PAPER_421).*
"""


def has_appendix(text: str) -> bool:
    return any(marker in text for marker in APPENDIX_MARKERS)


def get_paper_number(filename: str) -> int | None:
    """Extract numeric paper ID from filename like PAPER_042_..."""
    m = re.match(r'PAPER_(\d+)', filename, re.IGNORECASE)
    if m:
        return int(m.group(1))
    return None


def upgrade_early_papers():
    """Append UQFF framework reference to PAPER_001-099 that lack it."""
    files = sorted(os.listdir(WHITEPAPERS_DIR))
    upgraded = 0
    skipped = 0

    for fname in files:
        if not fname.endswith('.md'):
            continue
        num = get_paper_number(fname)
        if num is None or num < 1 or num > 99:
            continue

        path = os.path.join(WHITEPAPERS_DIR, fname)
        with open(path, 'r', encoding='utf-8', errors='replace') as f:
            text = f.read()

        if has_appendix(text):
            skipped += 1
            continue

        with open(path, 'a', encoding='utf-8') as f:
            f.write(UQFF_FRAMEWORK_APPENDIX)

        print(f"  ✅ Upgraded: {fname}")
        upgraded += 1

    print(f"\nEarly papers: {upgraded} upgraded, {skipped} already had appendix.")
    return upgraded


def upgrade_paper_050():
    """Add 26D implementation acknowledgment note."""
    matches = [
        f for f in os.listdir(WHITEPAPERS_DIR)
        if re.match(r'PAPER_050[_\.]', f) and f.endswith('.md')
    ]
    if not matches:
        print("  ⚠️  PAPER_050 not found — skipping.")
        return

    path = os.path.join(WHITEPAPERS_DIR, matches[0])
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()

    if "Implementation Note (v4.75)" in text or "26 compact dimensions" in text:
        print(f"  ℹ️  {matches[0]}: 26D note already present — skipping.")
        return

    # Insert after the ## Abstract heading (first occurrence)
    abstract_idx = text.find("## Abstract")
    if abstract_idx == -1:
        # Fallback: append
        text += PAPER_050_NOTE
    else:
        # Find end of Abstract section — insert before next ## heading
        after_abstract = abstract_idx + len("## Abstract")
        next_section = text.find("\n## ", after_abstract)
        if next_section == -1:
            text += PAPER_050_NOTE
        else:
            text = text[:next_section] + "\n" + PAPER_050_NOTE + text[next_section:]

    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)

    print(f"  ✅ PAPER_050: 26D implementation note added to {matches[0]}")


def upgrade_paper_200():
    """Add operational status table to Um taxonomy paper."""
    matches = [
        f for f in os.listdir(WHITEPAPERS_DIR)
        if re.match(r'PAPER_200[_\.]', f) and f.endswith('.md')
    ]
    if not matches:
        print("  ⚠️  PAPER_200 not found — skipping.")
        return

    path = os.path.join(WHITEPAPERS_DIR, matches[0])
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        text = f.read()

    if "Operational Implementation Status" in text or "## 13." in text:
        print(f"  ℹ️  {matches[0]}: Operational status section already present — skipping.")
        return

    text += PAPER_200_OPERATIONAL_STATUS

    with open(path, 'w', encoding='utf-8') as f:
        f.write(text)

    print(f"  ✅ PAPER_200: Operational status table appended to {matches[0]}")


def main():
    print("=" * 70)
    print("UQFF Whitepaper Batch Upgrade — v4.75")
    print("=" * 70)

    print("\n[1/3] Upgrading PAPER_001–099 (adding UQFF framework appendix)…")
    upgrade_early_papers()

    print("\n[2/3] Upgrading PAPER_050 (adding 26D implementation note)…")
    upgrade_paper_050()

    print("\n[3/3] Upgrading PAPER_200 (adding Um operational status table)…")
    upgrade_paper_200()

    print("\n✅ Done.\n")


if __name__ == "__main__":
    main()
