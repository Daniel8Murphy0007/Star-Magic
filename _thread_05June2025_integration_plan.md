# Integration Plan: thread_05June2025.txt
**Session 179 Part 3 | 2026-04-02**

## Overview
2 new whitepapers from thread_05June2025.txt (404-line Grok teaching session, June 5, 2025).
~90% of thread content already integrated in Sessions 120–178.

## New Papers
| PAPER | CP4 | Title |
|-------|-----|-------|
| PAPER_734 | #318 | LENR K_n Three-Scenario Calibration Constants |
| PAPER_735 | #319 | U_g2 Electron Shell Energy Eshell = c·νres·h(fSCm)·Ggeo |

## Counters After Integration
| Counter | Before | After |
|---------|--------|-------|
| PAPER count | 733 | 735 |
| CP4 classes | 309 | 311 |
| PDFs | 751 | 753 |
| Version | v5.35 | v5.36 |

---

## Step-by-Step Integration

### Step 1: Create PAPER_734 Whitepaper
File: `whitepapers/PAPER_734_LENR_Kn_ThreeScenario_Calibration_Constants.md`
- LENR K_n document three-scenario kη calibration table
- kη values: 2.75×10^8, 191, 6.06×10^-6
- ktrans = 5.26×10^44 solar corona transmutation
- CVW v2.0.0 compliant

### Step 2: Create PAPER_735 Whitepaper
File: `whitepapers/PAPER_735_Ug2_Electron_Shell_Energy_Eshell.md`
- U_g2 DPM electron shell energy Eshell = c·νres·h(fSCm)·Ggeo
- First explicit UQFF U_g2 equation
- CVW v2.0.0 compliant

### Step 3: Append CP4 Classes #318–319 to CondensedPhysics4.py
- `LENRKnScenarioCalibrationCalculator` (PAPER_734, CP4 #318)
- `Ug2ElectronShellEnergyCalculator` (PAPER_735, CP4 #319)

### Step 4: Generate PDFs
```powershell
# Generate PAPER_734 PDF
pandoc whitepapers/PAPER_734_LENR_Kn_ThreeScenario_Calibration_Constants.md `
  -o pdf/PAPER_734_LENR_Kn_ThreeScenario_Calibration_Constants.pdf `
  --pdf-engine=xelatex --from=markdown

# Generate PAPER_735 PDF
pandoc whitepapers/PAPER_735_Ug2_Electron_Shell_Energy_Eshell.md `
  -o pdf/PAPER_735_Ug2_Electron_Shell_Energy_Eshell.pdf `
  --pdf-engine=xelatex --from=markdown
```

### Step 5: Update Headers
Files to update:
- `HEADER_INTEGRATION_CHECKLIST.md` — update CP4/PAPER/PDF counters
- `ARCHITECTURE_FLOW_DIAGRAM.md` — add session 179 Part 3 footnote
- `CondensedPhysicsAggregator.py` — add imports for new CP4 classes

### Step 6: Git Commit + Push
```powershell
git add -A
git commit -m "Session 179 Part 3: thread_05June2025.txt audit — PAPER_734 LENR K_n 3-Scenario keta + PAPER_735 Ug2 Eshell; CP4 #318-#319 (311 total); 2 PDFs (753 total); 735/1000 (73.5%); v5.36"
git push origin master
```

---

## Integration Notes

### Why PAPER_734 (not already in PAPER_471):
PAPER_471 uses K_η as TARGET η values (1×10^13, 1×10^8, 7×10^-3) in the equation:
  η = Kη·exp(-[SSq]^n·2^6·e^{-π-t})·Um/ρvac
PAPER_734 documents the K_n document's MULTIPLICATIVE kη constants in the DIFFERENT form:
  η = kη·exp(-[SSq]·n/26)·exp(-(π-t)·Um/ρvac,[UA])
The kη values (2.75×10^8, 191, 6.06×10^-6) are calibration multipliers — different from PAPER_471 K_η values.
The ktrans=5.26×10^44 is completely new.

### Why PAPER_735 (not already in HydrogenResonanceShellCalculator):
HydrogenResonanceShellCalculator (CP2) uses:
  H_res = A_res·sin(2πf_res·t) + U_dp·SCm·k_nuc + S_shell
PAPER_735 documents:
  Eshell = c·νres·h(fSCm)·Ggeo
Different equation for U_g2 electron shell energy tied to SCm resonance frequency and geometry.

### Three UQFF Number Systems:
Already in PAPER_429 (CP4 #83) — NO new entry needed.
Thread predates Session 168 when PAPER_429 was integrated.

### ATLAS/CMS Higgs Data:
mH = 125.0±0.30 GeV combined — already in PAPER_718.
kHiggs = 1.79×10^18 — already in PAPER_718.
NO new entry needed.
