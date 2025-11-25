# Registry Analysis - November 22, 2025 @ 7:15 PM

## Executive Summary

**Phase 24: Registry Completion Analysis Complete**

### Current State
- **Total PhysicsTerm Classes:** 774 unique classes defined
- **Currently Registered:** 802 term names (includes aliases/helper methods)
- **Unregistered Classes:** 479 classes without registration
- **Wolfram Target:** 5,034+ physics terms (from previous WSTP search)

### Registration Gap
The registration count (802) exceeds the class count (774) because:
1. Some classes are registered with multiple aliases
2. Helper method terms are registered separately
3. The registration system uses descriptive names that differ from class names

### Two-Phase Plan

#### **Phase 1: Register Existing (479 classes)**
Complete registration of all 479 currently unregistered PhysicsTerm classes.

**Calculation:** 774 total - already integrated = 479 to add

#### **Phase 2: Wolfram Expansion (3,781 new classes)**
Create new PhysicsTerm classes to reach Wolfram 5,034+ target.

**Calculation:** 5,034 target - 774 existing - 479 Phase 1 = **3,781 new classes needed**

## Files Generated

| File | Count | Description |
|------|-------|-------------|
| `all_classes.txt` | 774 | All PhysicsTerm class names from MAIN_1_CoAnQi.cpp |
| `registered_terms.txt` | 802 | All registered term names from registerAllPhysicsTerms() |
| `unregistered_classes.txt` | 479 | Classes missing registration calls |

## Sample Unregistered Classes (First 20)

1. Abell2256_BuoyancyTerm
2. Abell2256_CompressedGTerm
3. Abell2256_CompressedTerm
4. Abell2256_DPM_resonanceTerm
5. Abell2256_FTerm
6. Abell2256_IntegrandTerm
7. Abell2256_LENR_termTerm
8. Abell2256_Q_waveTerm
9. Abell2256_ResonantTerm
10. Abell2256_SuperconductiveTerm
11. Abell2256_Ub1Term
12. Abell2256_UiTerm
13. Abell2256_X2Term
14. AetherVacDensity_A_mu_nuTerm
15. AetherVacDensity_Epsilon_factorTerm
16. AetherVacDensity_Rho_vac_ATerm
17. AetherVacDensity_T_s_contributionTerm
18. Andromeda_G_timeTerm
19. Antlia_G_timeTerm
20. Aquarius_G_timeTerm

## Categories of Unregistered Classes

Based on naming patterns:

### Astronomical Systems (computeG time-dependent gravity)
- Pattern: `{System}_G_timeTerm`
- Examples: Andromeda_G_timeTerm, Carina_G_timeTerm, ESO137_G_timeTerm
- Count: ~74 systems

### Astrophysical Methods (DPM, Q_wave, Ub1, Ui, X2, etc.)
- Pattern: `{System}_{Method}Term`
- Examples: Abell2256_DPM_resonanceTerm, VelaPulsar_Q_waveTerm
- Count: ~100 terms (10 systems × 10 methods)

### Helper Methods (Factors, base terms, derived quantities)
- Pattern: `{Module}_{Helper}Term`
- Examples: AetherVacDensity_Rho_vac_ATerm, BlackHole_M_bhInMsunTerm
- Count: ~56 terms

### System-Specific UQFF Terms
- Pattern: `{System}_{Physics}Term`
- Examples: Abell2256_BuoyancyTerm, ESO137_CompressedTerm
- Count: ~249 terms

## Known Constructor Parameter Issues

**3 Classes Commented Out (Require Parameters):**

1. **BuoyancyUQFFTerm** (line 20580)
   - Requires: `system_name` parameter
   - Constructor: `BuoyancyUQFFTerm(const std::string& systemName)`

2. **AstroSystemUQFFTerm** (line 20584)
   - Requires: `system_name` parameter
   - Constructor: `AstroSystemUQFFTerm(const std::string& systemName)`

3. **UQFFMasterTerm** (line 20588)
   - Requires: 7 parameters (system, sfr, wind, mag, f_Ub, M, r)
   - Constructor: `UQFFMasterTerm(string, double, double, double, double, double, double)`

## Next Actions

### Immediate (Phase 1)
1. ✅ Analyze registry gap (COMPLETE)
2. ⏸️ Generate registration calls for 479 unregistered classes
3. ⏸️ Handle constructor parameter requirements
4. ⏸️ Add registrations to registerAllPhysicsTerms() function
5. ⏸️ Update integration_status.json

### Future (Phase 2)
1. Retrieve Wolfram 5,034+ term list from previous WSTP session
2. Categorize terms by physics domain
3. Generate PhysicsTerm class templates
4. Implement compute() methods based on Wolfram physics equations
5. Register all new terms

## User Directive
**"no commits"** - Code changes only, no git operations

## Updated Files
- ✅ `MAIN_1_CoAnQi_integration_status.json` - Phase 24 status, corrected registry metrics
- ✅ `all_classes.txt` - Generated class name list
- ✅ `registered_terms.txt` - Generated registration list
- ✅ `unregistered_classes.txt` - Generated missing registration list
- ✅ `REGISTRY_ANALYSIS_NOV22_7PM.md` - This document

---
**Analysis Complete:** 2025-11-22 @ 19:15:00  
**Phase:** 24 - Registry Completion  
**Status:** Ready for Phase 1 registration additions
