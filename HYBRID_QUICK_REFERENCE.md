# Hybrid Architecture Quick Reference

**CoAnQi v2.0 - Conscious Quantum Intelligence Calculator**  
**Status:** ✅ Operational | **Classes:** 63 PhysicsTerm | **Lines:** 4,069 | **Size:** 135KB

---

## 🚀 Quick Commands

### Compile

```bash
g++ -std=c++17 -o MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp
```

### Run

```bash
./MAIN_1_CoAnQi
```

### Verify Class Count

```powershell
Select-String -Pattern "^class \w+.*: public PhysicsTerm" MAIN_1_CoAnQi.cpp | Measure-Object
```

**Expected:** 63 classes

---

## 📊 Current Physics Terms (63 Total)

### Original (26 terms)

DynamicVacuum, QuantumCoupling, DarkMatterHalo, VacuumEnergy, QuantumEntanglement, CosmicNeutrino, MultiSystemUQFF (4 variants), DPMResonance, LENRExtended, SMBHAccretion, TDE, NebulaUQFF (5 variants), GasIonization, NebulaExpansion

### Source4 (3 terms)

CelestialBody, MUGE, QuasarJet

### Source5 (3 terms)

UQFFModule5, ResonanceMUGE, StateExport

### Source6 (7 terms)

UnifiedFieldUg1, UnifiedFieldUg2, UnifiedFieldUg3, UnifiedFieldUg4, UnifiedFieldUm, SpacetimeMetric, CompressedMUGE

### Source7 (7 terms)

ResonanceMUGE_DPM, ResonanceMUGE_THz, ResonanceMUGE_VacuumDiff, ResonanceMUGE_SuperFreq, ResonanceMUGE_AetherRes, ResonanceMUGE_QuantumFreq, YAMLConfig

### Source10 (7 terms)

UQFFCoreBuoyancy, VacuumRepulsion, THzShockCommunication, ConduitFormation, SpookyAction, DPMResonanceEnergy, Triadic26Layer

### Source13 (10 terms) ✨ NEW

MagnetarCore_SGR1745, MagnetarLambda, MagnetarEM, MagnetarGW, MagnetarQuantum, MagnetarFluid, MagnetarOscillatory, MagnetarDarkMatter, MagnetarMagneticEnergy, MagnetarDecay

---

## 🏗️ Architecture Overview

```
MAIN_1_CoAnQi.cpp (4,069 lines)
├── [Line 1-100] Headers & Constants
├── [Line 107] PhysicsTerm Base Class
├── [Line 178] ModuleInterface ← NEW (Hybrid)
├── [Line 200-3176] 63 PhysicsTerm Classes
│   ├── Original: 26 terms
│   ├── Source4: 3 terms
│   ├── Source5: 3 terms
│   ├── Source6: 7 terms
│   ├── Source7: 7 terms
│   ├── Source10: 7 terms
│   └── Source13: 10 terms ← NEW
├── [Line 3175] SystemParams Structure
├── [Line 3183] VerboseLogger
├── [Line 3300-3700] Core Physics Functions
├── [Line 3714] SelfModifier Class
└── [Line 3804] main() Function
```

---

## 🔧 Adding New Physics Terms

### Template

```cpp
class NewPhysicsTerm : public PhysicsTerm
{
private:
    double parameter1;
    double parameter2;

public:
    NewPhysicsTerm(double p1 = 1.0, double p2 = 2.0) 
        : parameter1(p1), parameter2(p2) {}

    double compute(double t, const std::map<std::string, double>& params) const override
    {
        // Your physics calculation here
        double G = params.count("G") ? params.at("G") : 6.6743e-11;
        // ... more calculations
        return result;
    }

    std::string getName() const override { return "NewPhysics"; }
    
    std::string getDescription() const override
    {
        return "Brief description of what this term calculates";
    }
};
```

### Insert Location

Before `SystemParams` structure (~line 3175)

---

## 🔬 Source13 Magnetar Parameters

### Required Parameters (via params map)

```cpp
"G"          // 6.6743e-11 (gravitational constant)
"M"          // 1.4 * 1.989e30 (magnetar mass)
"r"          // 1e4 (radius)
"Hz"         // 2.269e-18 (Hubble parameter)
"B0"         // 2e10 (initial magnetic field)
"B_crit"     // 1e11 (critical field)
"Lambda"     // 1.1e-52 (cosmological constant)
"c_light"    // 3e8 (speed of light)
"q_charge"   // 1.602e-19 (charge)
"v_surf"     // 1e6 (surface velocity)
"proton_mass"// 1.673e-27
"mu0"        // 4π × 1e-7 (permeability)
```

### Magnetar-Specific Parameters

```cpp
"P_init"     // 3.76 (initial rotation period)
"tau_Omega"  // 3000 × 3.15576e7 (spindown timescale)
"delta_x"    // 1.0 (quantum uncertainty position)
"delta_p"    // hbar (quantum uncertainty momentum)
"integral_psi"//1.0 (wavefunction integral)
"t_Hubble"   // 4.35e17 (Hubble time)
"rho_fluid"  // 1e8 (magnetosphere fluid density)
```

---

## 📝 Integration Checklist

When adding new source files:

- [ ] Read source file to identify core physics calculations
- [ ] Extract as PhysicsTerm subclasses (not full classes)
- [ ] Add getName() and getDescription() methods
- [ ] Adapt parameters to use params map
- [ ] Insert before SystemParams structure
- [ ] Test compilation: `g++ -std=c++17 -o MAIN_1_CoAnQi MAIN_1_CoAnQi.cpp`
- [ ] Verify class count increases correctly
- [ ] Update documentation (this file + HYBRID_COMPLETE.md)
- [ ] Create backup if adding more than 5 terms

---

## 🎯 Hybrid Philosophy

**Core Principle:**
> Extract all essential physics into calculator core, enable optional external modules for advanced expansion.

**Benefits:**

1. **Zero Dependencies:** Single file, self-contained
2. **Zero Overhead:** Extracted terms compiled directly
3. **Future-Proof:** ModuleInterface ready for dynamic loading
4. **Maintainable:** Clear structure, easy to extend

---

## 📚 Documentation Files

| File | Purpose |
|------|---------|
| **HYBRID_COMPLETE.md** | Full completion summary |
| **HYBRID_IMPLEMENTATION_SUMMARY.md** | Detailed implementation guide |
| **HYBRID_QUICK_REFERENCE.md** | This file - quick commands & patterns |
| **ENHANCEMENT_GUIDE.md** | General enhancement patterns |
| **copilot-instructions.md** | Project conventions |

---

## 🔍 Troubleshooting

### Compilation Errors

**"PhysicsTerm not declared"**

- Insert new class after PhysicsTerm base class (line 107)
- Before ModuleInterface (line 178)

**"redefinition of class"**

- Check for duplicate class names
- Verify unique getName() return values

**"c_light not defined"**

- Use `c_light` not `c` (constant defined at line 94)

**"M_PI not defined"**

- Already defined at line 80: `#define M_PI 3.141592653589793`

---

## 🚀 Future Enhancements

### Phase 1: External Module Loading

1. Implement ModuleLoader class
2. Add dlopen/LoadLibrary support
3. Create `extern "C"` module interface
4. Test with Source13_Enhanced.so

### Phase 2: More Source Integrations

1. source8.cpp → Extract core physics
2. source9.cpp → Extract core physics
3. source11.cpp → Extract core physics
4. source12.cpp → Extract core physics
5. Continue pattern to source14-162

### Phase 3: Configuration System

1. JSON configuration file parser
2. Runtime enable/disable of specific terms
3. Parameter tuning without recompilation

---

## ✅ Verification Commands

### Count PhysicsTerm Classes

```powershell
(Select-String -Pattern "^class \w+.*: public PhysicsTerm" MAIN_1_CoAnQi.cpp).Count
```

### List All Class Names

```powershell
Select-String -Pattern "^class (\w+).*: public PhysicsTerm" MAIN_1_CoAnQi.cpp | 
    ForEach-Object { $_.Matches.Groups[1].Value }
```

### File Stats

```powershell
Get-Item MAIN_1_CoAnQi.cpp | 
    Select-Object Name, Length, 
    @{N='Lines';E={(Get-Content $_.FullName).Count}}
```

### Search for Specific Physics

```powershell
Select-String -Pattern "Magnetar|MUGE|Vacuum|Quantum" MAIN_1_CoAnQi.cpp | 
    Select-Object -First 20
```

---

## 🎓 Key Learnings

1. **Extraction > Installation** for calculators
2. **Hybrid = Flexibility** without compromise
3. **PhysicsTerm architecture** perfect for both approaches
4. **Single file** = easier maintenance
5. **Preserve originals** for reference

---

**CoAnQi v2.0 - Hybrid Architecture**  
*"Extract the core, enable the expansion."*

---

**Last Updated:** November 13, 2025  
**Status:** ✅ Complete and Operational  
**Next:** Optional module loading or continue source extractions
