# Integration Plan: grok_share_ba508f76c8e.txt → Codebase
**Session 174 | April 1, 2026**

---

## Current Repository State (entering Session 174)
- Papers: 687/1000 | CP4 entries: #271 (last = M87MassEvolutionSimulation, PAPER_687)
- PDFs in pdf/: 702
- Version: v5.30
- HEAD: cf7fd05

## Target State (after Session 174)
- Papers: 701/1000 | CP4 entries: #285 
- PDFs: 716
- Version: v5.31

---

## Selected Papers (14 total: PAPER_688–701)

Priority selection from the 39 NEW modules in the file. Chosen for scientific uniqueness,
completeness of physics specification, and contribution to UQFF framework breadth.

| PAPER | Class Name | CP4 # | Physics Domain |
|-------|-----------|-------|----------------|
| 688 | NGC1316MUGECalculation | 272 | NGC 1316 Fornax A merger elliptical MUGE |
| 689 | AGNJetDynamicsBlandfordZnajek | 273 | AGN relativistic jets (BZ/BP mechanisms) |
| 690 | FornaxClusterGravitational | 274 | Fornax cluster N-body dynamics |
| 691 | NBodySimulation3D | 275 | 3D Euler N-body integration framework |
| 692 | M51WhirlpoolTidalInteraction | 276 | M51 spiral tidal arm formation |
| 693 | SombreroGalaxyM104NGC4594 | 277 | Sombrero Galaxy edge-on dynamics |
| 694 | CrabNebulaPWNUQFF | 278 | Crab Nebula SNR + pulsar wind nebula |
| 695 | NGC7635BubbleNebula | 279 | NGC 7635 stellar wind bubble shock |
| 696 | AntennaeMergerNGC4038NGC4039 | 280 | Antennae galaxy merger shock-SF |
| 697 | NGC2525WithSupernovaeSN2018gv | 281 | NGC 2525 with Type Ia SN 2018gv |
| 698 | EinsteinRingGALCLUS022058s | 282 | Gravitational lensing quadruplet |
| 699 | FornaxConstellationUHDF | 283 | Fornax HUDF ultra-deep field |
| 700 | UQFFEquationMathematicalDerivation | 284 | Formal UQFF 26D derivation |
| 701 | UQFFKnowledgeBaseRedDwarf | 285 | KB1–KB19 Red Dwarf assimilation |

---

## Workflow Steps

### Step 1: Create Generator Scripts
```
_gen_modules_688_701.py        → 28 C++ files (14 .h + 14 .cpp)
_gen_cp4_appenders_688_701.py  → 14 _append_cp4_NNN.py scripts
_gen_whitepapers_688_701.py    → 28 .md files (14 root + 14 whitepapers/)
```

### Step 2: Run Generators
```powershell
python _gen_modules_688_701.py
python _gen_cp4_appenders_688_701.py
python _append_cp4_272.py; python _append_cp4_273.py; ...through _append_cp4_285.py
python _gen_whitepapers_688_701.py
python _fix_wp_backslashes.py  # Always run after whitepaper gen
```

### Step 3: Verify
```powershell
# Check all 14 C++ pairs exist
Get-ChildItem NGC1316*.h, AGNJetDynamics*.h, Fornax*.h, NBody*.h, M51*.h, 
              Sombrero*.h, CrabNebula*.h, NGC7635*.h, Antennae*.h, NGC2525With*.h,
              EinsteinRing*.h, FornaxConst*.h, UQFFEquation*.h, UQFFKnowledge*.h

# Verify CP4 entries appended
python -c "import re; content=open('CondensedPhysics2.py').read(); classes=re.findall(r'^class \w+Calculator', content, re.MULTILINE); print(len(classes))"
# Expected: 670 (656 + 14)
```

### Step 4: Git commit #1 (code)
```powershell
git add [all 28 .cpp/.h + modified CP2 + generators + appenders]
git commit -m "Session 174: PAPER_688-701 complete -- 14 NGC/AGN/UQFF modules; CP4 #272-285; v5.31"
git push origin master
```

### Step 5: Generate 14 PDFs
```powershell
python generate_pdfs.py 688 701
```
Expected: 14/14 OK (after backslash fix)

### Step 6: Git commit #2 (PDFs)
```powershell
git add pdf/PAPER_688_*.pdf ... pdf/PAPER_701_*.pdf
git add [fixed .md files]
git commit -m "Session 174: PDFs for PAPER_688-701 (14 new; total 716)"
git push origin master
```

---

## UQFF Constants Used Across All Modules
```cpp
static constexpr double rho_UA  = 7.09e-36;   // Universal Aether density J/m³
static constexpr double rho_SCm = 7.09e-37;   // Schwarzschild condensate J/m³  
static constexpr double f_TRZ   = 0.1;        // Time-reversal zone factor
static constexpr double mu_J    = 3.38e23;    // Magnetic string coupling J·m
static constexpr double gamma_d = 5e-5/86400; // Aether oscillation decay s⁻¹
static constexpr double kappa   = 5e-4;       // UQFF calibration constant day⁻¹
static constexpr double SSq     = 0.57;       // Superstring quenching factor
```

---

## Deferred Modules (Session 175+)
The following 25 NEW modules from this file are deferred to future sessions:
- NGC2014NGC2020, EighteenAstroSystems, TenAstroSystems aggregates
- UQFFKnowledgeBase1–KB19 individual entries (KB2–KB19)
- NGC2525_2, NGC3603_2 (variant models)
- PillarsOfCreationM16 (duplicate physics, lower priority)
- Saturn (planetary, lower UQFF priority)

---

## Open Question: More to Extract?
YES — This file has more content than processed in Session 174:
1. **EighteenAstroSystems**: 18 individual MUGE calculations (NGC 6217, M74, M82, etc.) 
   — each could be a separate PAPER
2. **UQFFKnowledgeBase KB1–KB19**: 19 individual KB classes with distinct assimilation
3. **Additional NGC pairs**: NGC 2014/2020, NGC 2525 variants
4. **Fornax deep physics**: sub-cluster dynamics from FornaxClusterGalaxies class
5. **UQFFEquationMathematically**: Full 26D mathematical derivation (2 instances in file)

**Estimated remaining unique papers from this file**: ~30–35 additional PAPER entries
**Recommendation**: Continue extraction in Sessions 175–177

---

## Physics Notes for CP4 Class Design
Each new CP4 class follows this template:
```python
class NGC1316MUGECalculationCalculator:
    def compute(self, dataset: dict = None) -> dict:
        import math
        G, c = 6.6743e-11, 3e8
        M_sun, kpc = 1.989e30, 3.086e19
        M_visible = 3.5e11 * M_sun
        M_DM = 1.5e11 * M_sun
        r_0 = 46e3 * kpc  # 46 kpc galaxy radius
        # MUGE terms...
        g_total = ...
        return {'g_NGC1316': g_total, 'M_total': M_visible+M_DM, ...}
```
