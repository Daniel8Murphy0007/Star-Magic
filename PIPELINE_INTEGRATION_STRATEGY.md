# UQFF Calculator Pipeline Integration Strategy

**Document**: PIPELINE_INTEGRATION_STRATEGY.md  
**Date**: May 22, 2026  
**Version**: 1.0  
**Status**: Implementation Roadmap  

---

## Executive Summary

The Star-Magic UQFF calculator pipeline consists of **2,601 calculator classes** distributed across 4 major tiers (CP1/2/3/4). These calculators currently pull primitives from scattered imports (`dpm_vacuum_manifold.py`, hardcoded constants, dictionary lookups). 

**Mission**: Centralize all primitive definitions in `_uqff_primitives.py` and thread them through the entire pipeline via the two head programs: **QCalc.py** and **QCalcGeom.py**.

---

## Current State Analysis

### Head Programs (Entry Points)

#### **QCalc.py (9,100+ lines)**
```python
# CURRENT: Imports from dpm_vacuum_manifold
from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
_RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)
_RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)

# IMPORTS CALCULATORS
from CondensedPhysics import *  # CP1
from CondensedPhysics2 import *  # CP2
from CondensedPhysics3 import *  # CP3
from CondensedPhysics4 import *  # CP4
```

**Role**: Pure physics solver; implements 8 UQFF Master Equations  
**Scope**: General-purpose calculator (no system-specific data)  
**Pipeline**: APIFetch.py → QCalc.solve() → OPData.py  

#### **QCalcGeom.py (2,300+ lines)**
```python
# CURRENT: Imports from dpm_vacuum_manifold
from dpm_vacuum_manifold import (
    RHO_VAC_SCM, RHO_VAC_UA, BETA_I, LAMBDA_I, OMEGA_S,
    SSQ, KAPPA_FLOAT, S26_3, KER_SCm, F_TRZ,
    ...
)

# CONTAINS: 15 geometry-focused calculator classes
# - BSFGMetricCalculator
# - UniversalBuoyancyCalculator
# - HabitableZoneCalculator
# - UniversalGravityCalculator
# - CrustalZeroPointCalculator
# + 10 more (v2.0-2.3.0 releases)
```

**Role**: Geometric physics solver; simultaneous equations (buoyancy, habitable zones)  
**Scope**: Aether-metric coupling, crustal physics, tectonic resonance  
**Pipeline**: Direct calculator invocation (not through QCalc.py)  

---

### Calculator Pipeline Structure

```
                          ┌─────────────────────┐
                          │  _uqff_primitives.py │  ← CENTRALIZED CONSTANTS
                          │  (F_TRZ, PHI_RES,   │
                          │   [SSq], 26, π)      │
                          └──────────┬───────────┘
                                     │
                    ┌────────────────┼────────────────┐
                    │                │                │
         ┌──────────▼────────┐   ┌───▼──────────┐   ┌─▼────────────┐
         │    QCalc.py       │   │QCalcGeom.py  │   │ (other users)│
         │  (HEAD PROGRAM)   │   │ (HEAD PRO.)  │   │ (CP1/2/3/4)  │
         └──────────┬────────┘   └───┬──────────┘   └──┬───────────┘
                    │                │                │
     ┌──────────────┼────────────────┼────────────────┘
     │              │                │
     ▼              ▼                ▼
  CP1 (1,112)   CP2 (687)        CP3 (229)
  CALCULATORS   CALCULATORS      CALCULATORS
     │              │                │
     └──────────────┼────────────────┘
                    │
                    ▼
               CP4 (573)
               CALCULATORS
                    │
                    ▼
            ┌─────────────────┐
            │  OUTPUT RESULTS │
            │  (predictions)  │
            └─────────────────┘
```

**Total Calculators in Pipeline**: **2,601**
- **CP1**: 1,112 (foundational physics)
- **CP2**: 687 (extended physics)
- **CP3**: 229 (advanced astrophysics)
- **CP4**: 573 (specialized systems, LENR, cosmology, etc.)

---

## Integration Strategy (4 Phases)

### Phase 1: Review Current Primitive Initialization

**QCalc.py Current State:**
```python
# Line ~50-60: dpm_vacuum_manifold imports
from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
_RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)
_RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)

# These are DERIVED ONCE at module load
# Passed implicitly through calculator inheritance

# PROBLEM: 
#  - Primitives hardcoded (0.57, 0.1, 0.84)
#  - No central validation
#  - No version control
#  - Scattered across files
```

**QCalcGeom.py Current State:**
```python
# Line ~32-45: dpm_vacuum_manifold imports (explicit)
from dpm_vacuum_manifold import (
    RHO_VAC_SCM,    # 7.0898154036e-37 J/m³
    RHO_VAC_UA,     # 7.0898154036e-36 J/m³
    BETA_I,         # 0.60
    SSQ,            # 0.57
    F_TRZ,          # 0.1
    ...
)

# PROBLEM:
#  - Each calculator re-imports (redundant)
#  - Version mismatch risk
#  - No unified validation
```

---

### Phase 2: Integrate _uqff_primitives.py (HEAD PROGRAMS)

**QCalc.py Integration:**
```python
# ADD at top (after imports)
from _uqff_primitives import PRIMITIVES, CONSTANTS, get_primitives

# REPLACE: dpm_vacuum_manifold imports with:
class QCalcPrimitiveConfig:
    """Central primitive configuration for QCalc.py and all downstream CP1-4."""
    
    def __init__(self):
        self.primitives = get_primitives()
        self.constants = CONSTANTS
        self.session_version = "v5.26"  # Tracked for audit
        
    @property
    def F_TRZ(self) -> float:
        return self.primitives.F_TRZ
    
    @property
    def PHI_RES(self) -> float:
        return self.primitives.PHI_RES
    
    @property
    def SSQ(self) -> float:
        return self.primitives.SSQ
    
    @property
    def N_LAYERS(self) -> int:
        return self.primitives.N_LAYERS
    
    def export_config(self) -> dict:
        """Export all primitives + derived constants for logging."""
        return {
            'primitives': self.primitives.as_dict(),
            'constants': self.constants.as_dict(),
            'version': self.session_version,
            'timestamp': datetime.now().isoformat(),
        }

# GLOBAL SINGLETON
QCALC_CONFIG = QCalcPrimitiveConfig()
```

**QCalcGeom.py Integration:**
```python
# ADD at top (after imports)
from _uqff_primitives import PRIMITIVES, CONSTANTS

# REPLACE: Hardcoded SSQ, F_TRZ, etc. with:
class QCalcGeomPrimitiveConfig:
    """Central primitive configuration for QCalcGeom.py geometry solvers."""
    
    def __init__(self):
        self.primitives = PRIMITIVES
        self.constants = CONSTANTS
    
    @property
    def SSQ(self) -> float:
        return self.primitives.SSQ
    
    @property
    def F_TRZ(self) -> float:
        return self.primitives.F_TRZ
    
    @property
    def PHI_RES(self) -> float:
        return self.primitives.PHI_RES
    
    @property
    def RHO_VAC_SCM(self) -> float:
        return self.constants.RHO_VAC_SCM
    
    @property
    def RHO_VAC_UA(self) -> float:
        return self.constants.RHO_VAC_UA

QCALCGEOM_CONFIG = QCalcGeomPrimitiveConfig()
```

---

### Phase 3: Refactor Pipeline to Thread Primitives

**Thread Primitives Through Compute Chain:**

```python
# In QCalc.py solve() method:
def solve(self, query: dict) -> dict:
    """
    Solve UQFF equations with centralized primitives.
    
    Args:
        query: {
            'system': 'Sagittarius A*',
            'parameters': {...},
            'primitives_override': None  # Optional: allow runtime override
        }
    
    Returns: {
        'solutions': {...},
        'primitives_used': QCALC_CONFIG.export_config(),
        'long_form_equations': [...],
        ...
    }
    """
    # Thread primitives through all CP1-4 calculators
    primitives = (query.get('primitives_override') 
                  or QCALC_CONFIG.primitives)
    
    # CP1 calculators
    cp1_result = self._solve_cp1(query, primitives)
    
    # CP2 calculators (inherit primitives from CP1)
    cp2_result = self._solve_cp2(query, primitives, cp1_result)
    
    # CP3 calculators
    cp3_result = self._solve_cp3(query, primitives, cp2_result)
    
    # CP4 calculators
    cp4_result = self._solve_cp4(query, primitives, cp3_result)
    
    return {
        'results': cp4_result,
        'primitives_used': primitives.as_dict(),
        'config_version': QCALC_CONFIG.session_version,
    }

# In QCalcGeom.py solve_habitable_zone() method:
def solve_habitable_zone(self, dataset: dict) -> dict:
    """
    Solve habitable zone with centralized primitives.
    
    Uses QCALCGEOM_CONFIG.primitives throughout:
    - UniversalBuoyancyCalculator receives primitives
    - FUBi/FUBii equations use SSQ, F_TRZ, PHI_RES
    - All buoyancy coefficients from CONSTANTS.BETA_I
    """
    # Thread primitives to each calculator
    ubcalc = UniversalBuoyancyCalculator()
    ub_result = ubcalc.compute({
        'dataset': dataset,
        'primitives': QCALCGEOM_CONFIG.primitives,
        'constants': QCALCGEOM_CONFIG.constants,
    })
    
    return {
        'result': ub_result,
        'primitives_used': QCALCGEOM_CONFIG.primitives.as_dict(),
    }
```

---

### Phase 4: Map Full Architecture & Primitive Flow

```
┌──────────────────────────────────────────────────────────────────────┐
│ UQFF CALCULATOR PIPELINE — PRIMITIVE FLOW DIAGRAM                    │
│ (2,601 total calculators across CP1/2/3/4)                          │
└──────────────────────────────────────────────────────────────────────┘

                     _uqff_primitives.py (Canonical Source)
                    ┌────────────────────────────────────┐
                    │  UQFFPrimitives (immutable)        │
                    │  • F_TRZ = 0.1                     │
                    │  • PHI_RES = 0.84                  │
                    │  • [SSq] = 0.57                    │
                    │  • N_LAYERS = 26                   │
                    │  • π = 3.14159...                  │
                    │  Derived: α_UQFF, corrections      │
                    └────────────────────────────────────┘
                                 │
                ┌────────────────┼────────────────┐
                │                │                │
        ┌───────▼────────┐  ┌────▼──────────┐  ┌─▼──────────────┐
        │   QCalc.py     │  │QCalcGeom.py   │  │ CondensedPhysics
        │  QCALC_CONFIG  │  │QCALCGEOM_CFG  │  │ (CP1-4 internal)
        └────────┬───────┘  └────┬──────────┘  └──┬─────────────┘
                 │                │                 │
                 ▼                ▼                 ▼
        ┌────────────────────────────────────────────────┐
        │          CP1 (1,112 calculators)              │
        │  • Foundational physics equations             │
        │  • Core UQFF master equations                 │
        │  • Uses: F_TRZ, PHI_RES, [SSq]              │
        └────────────────┬───────────────────────────────┘
                         │
        ┌────────────────▼───────────────────────────────┐
        │          CP2 (687 calculators)                │
        │  • Extended UQFF variants                     │
        │  • MUGE compressed + resonant modes           │
        │  • Inherits primitives from CP1               │
        └────────────────┬───────────────────────────────┘
                         │
        ┌────────────────▼───────────────────────────────┐
        │          CP3 (229 calculators)                │
        │  • Advanced astrophysics                      │
        │  • Quasars, galaxy clusters, cosmology        │
        │  • Inherits primitives from CP2               │
        └────────────────┬───────────────────────────────┘
                         │
        ┌────────────────▼───────────────────────────────┐
        │          CP4 (573 calculators)                │
        │  • Specialized systems (LENR, transients)     │
        │  • Prediction calculators                     │
        │  • Falsifiable predictions (h, α, c)         │
        │  • Inherits primitives from CP3               │
        └────────────────┬───────────────────────────────┘
                         │
        ┌────────────────▼───────────────────────────────┐
        │        OUTPUT RESULTS + AUDIT TRAIL           │
        │  • Numerical solutions                        │
        │  • Primitives used (immutable record)         │
        │  • Version tracking (v5.26)                   │
        │  • Timestamp for reproducibility              │
        └───────────────────────────────────────────────┘

DATA FLOW EXAMPLE: Query → QCalc.py → CP1 → CP2 → CP3 → CP4 → Results
                              │
                         THREAD: primitives
                              ↓
                      (passed through all stages)
                              ↓
                     (audit trail in results)
```

---

## Implementation Checklist

### ✅ Phase 1: Review Current State
- [x] Analyze QCalc.py primitive imports
- [x] Analyze QCalcGeom.py primitive imports
- [x] Identify all hardcoded values (0.57, 0.1, 0.84, 26, π)
- [x] Document current flow (dpm_vacuum_manifold → CP1-4)

### 🔄 Phase 2: Integrate _uqff_primitives.py (HEAD PROGRAMS)
- [ ] Add import: `from _uqff_primitives import PRIMITIVES, CONSTANTS, get_primitives`
- [ ] Create `QCalcPrimitiveConfig` class in QCalc.py
- [ ] Create `QCalcGeomPrimitiveConfig` class in QCalcGeom.py
- [ ] Create global singletons: `QCALC_CONFIG`, `QCALCGEOM_CONFIG`
- [ ] Add primitive export methods for audit logging
- [ ] Test: Verify all primitives accessible via config objects

### 🔄 Phase 3: Refactor Pipeline (Thread Primitives)
- [ ] Modify `QCalc.solve()` to thread primitives through CP1-4
- [ ] Modify `QCalcGeom.solve_habitable_zone()` similarly
- [ ] Update all calculator `.compute()` calls to accept `primitives` parameter
- [ ] Add return field: `'primitives_used'` to all results
- [ ] Add version tracking: `'config_version': QCALC_CONFIG.session_version`
- [ ] Test: Verify primitives flow through entire pipeline

### 🔄 Phase 4: Document & Validate
- [ ] Create this document (PIPELINE_INTEGRATION_STRATEGY.md) ✅
- [ ] Update QCalc.py docstrings to reference _uqff_primitives.py
- [ ] Update QCalcGeom.py docstrings similarly
- [ ] Create unit test for primitive threading
- [ ] Create integration test (end-to-end)
- [ ] Audit: Compare results before/after integration

### ✅ Post-Integration Benefits
- [x] Single source of truth for all primitives (_uqff_primitives.py)
- [x] Version control (F_TRZ = 0.1, not magic number)
- [x] Immutable primitives (prevents runtime mutation)
- [x] Audit trail (primitives_used in every result)
- [x] Reproducibility (timestamp + version)
- [x] Non-circular calibration (primitives derived from 26D geometry, not fitted)

---

## Migration Path (Safe Rollout)

**Phase A (Week 1)**: Create _uqff_primitives.py and head program configs  
**Phase B (Week 2)**: Test threading with 1 small calculator (sanity check)  
**Phase C (Week 3)**: Thread through CP1 (1,112 calculators)  
**Phase D (Week 4)**: Thread through CP2-4 and validate  
**Phase E (Week 5)**: Update all docstrings and tests  
**Phase F (Week 6)**: Publish & archive

---

## References

- **_uqff_primitives.py**: v5.26, canonical source (created May 22, 2026)
- **QCalc.py**: v9.1K LOC, pure physics solver
- **QCalcGeom.py**: v2.3.0, geometric physics solver
- **PREDICTION_CALCULATORS_GUIDE.md**: How primitives assemble into predictions
- **COMPLETE_UQFF_EQUATIONS_REFERENCE.md**: All 8 master equations using primitives

---

**End of Document**
