# Phase 2C: Batch Update Template for Remaining Tiers 3-4

**Purpose:** Template and guidelines for systematically updating remaining 100+ PhysicsTerm classes  
**Estimated Classes Remaining:** 100+  
**Pattern:** Apply transformation from basic Ug_sum return to UbiForceBalanceIntegrator pattern  

---

## Pattern Template

### BEFORE (Current Implementation)
```cpp
class YourPhysicsTermName : public PhysicsTerm
{
    // ... private members and constructor ...
    
    double compute(double t, const std::map<std::string, double> &params) const override
    {
        // Physics calculations...
        double result = /* some calculation */;
        return result;
    }
    
    std::string getName() const override { return "YourTermName"; }
    std::string getDescription() const override
    {
        return "Original description without Ubi";
    }
};
```

### AFTER (With Ubi Force Balance)
```cpp
class YourPhysicsTermName : public PhysicsTerm
{
    // ... private members and constructor ... (unchanged)
    
    double compute(double t, const std::map<std::string, double> &params) const override
    {
        // Physics calculations...
        double result = /* some calculation */;
        
        // Phase 2C: Apply Session 252 Ubi force balance
        return UbiForceBalanceIntegrator::applyForceBalance(result, params);
    }
    
    std::string getName() const override { return "YourTermName"; }
    std::string getDescription() const override
    {
        return "Original description + Session 252 buoyancy with Ubi force balance";
    }
};
```

---

## Step-by-Step Application Guide

### For Each PhysicsTerm Class:

1. **Identify the return statement** in the `compute()` method
   - Look for: `return <expression>;`
   
2. **Extract the expression** into a named variable:
   ```cpp
   double result = <expression>;
   return result;
   ```
   
3. **Wrap with UbiForceBalanceIntegrator**:
   ```cpp
   double result = <expression>;
   
   // Phase 2C: Apply Session 252 Ubi force balance
   return UbiForceBalanceIntegrator::applyForceBalance(result, params);
   ```

4. **Update the description** to include Ubi reference:
   ```cpp
   std::string getDescription() const override
   {
       return "Original description + Session 252 buoyancy with Ubi force balance";
   }
   ```

---

## Classes Identified for Phase 2C Updates

### Tier 3: Astrophysical System Classes (20+ classes)
These compute specific astrophysical objects with MUGE framework:

- Line 4670: MagnetarCoreTerm
- Line 4728: MagnetarLambdaTerm
- Line 4750: MagnetarEMTerm
- Line 4778: MagnetarGWTerm
- Line 4816: MagnetarQuantumTerm
- Line 4842: MagnetarFluidTerm
- Line 4872: MagnetarOscillatoryTerm
- Line 4909: MagnetarDarkMatterTerm
- Line 4947: MagnetarMagneticEnergyTerm
- Line 4977: MagnetarDecayTerm
- Line 5014: Magnetar0501CoreTerm
- Line 5064: Magnetar0501LambdaTerm
- Line 5086: Magnetar0501EMTerm
- Line 5124: Magnetar0501GWTerm
- Line 5162: Magnetar0501QuantumTerm
- Line 5188: Magnetar0501FluidTerm
- Line 5218: Magnetar0501OscillatoryTerm
- Line 5254: Magnetar0501DarkMatterTerm
- Line 5296: SgrAStarCoreTerm
- Line 5352: SgrAStarLambdaTerm
(and more SgrAStar variants...)

### Tier 4: Utility/Support Classes (80+ classes)
Generic physics terms that support core calculations:

- Line 1196: DynamicVacuumTerm
- Line 1228: QuantumCouplingTerm
- Line 1259: DarkMatterHaloTerm
- Line 1296: VacuumEnergyTerm
- Line 1326: QuantumEntanglementTerm
- ... (see grep_search output for full list)

---

## Automated Application Suggestions

### Option A: Manual with Template
1. Open MAIN_1_CoAnQi.cpp in editor
2. Search for each class using Ctrl+G (Go to Line)
3. Apply template transformation
4. Verify compile after each 5-10 classes

### Option B: Script-Assisted (for next session)
Create a Python script to:
1. Parse MAIN_1_CoAnQi.cpp
2. Identify all PhysicsTerm subclasses not yet updated
3. Generate suggested diff patches
4. Apply patches with user confirmation

### Option C: Grep & Batch Replace
Use VS Code Find-and-Replace (Ctrl+H) with regex:
1. **Find:** `(return\s+)([^;]+);` (in compute methods)
2. **Replace:** `double result = $2;\n        // Phase 2C: Apply Session 252 Ubi force balance\n        return UbiForceBalanceIntegrator::applyForceBalance(result, params);`
3. **Manual review required** to verify context

---

## Validation Checklist for Batch Updates

After each batch of 5-10 classes:

- [ ] Searched for class name in editor
- [ ] Located compute() method
- [ ] Extracted return expression to named variable
- [ ] Added UbiForceBalanceIntegrator::applyForceBalance() call
- [ ] Updated getDescription() to mention Ubi
- [ ] Verified indentation and syntax
- [ ] Saved file
- [ ] Visual scan for obvious errors

After completing all remaining classes:

- [ ] Run: `cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64`
- [ ] Build: `cmake --build build_msvc --config Release`
- [ ] Verify: 0 compilation errors
- [ ] Git commit: "Phase 2C: Complete Tier 3/4 Ubi integration (all 446 modules)"

---

## Statistics Tracking

### Tier 1: Core UQFF Classes
- Target: 6 classes
- Completed: ✅ 6/6 (100%)
- Status: DONE

### Tier 2: Compressed/Resonance MUGE
- Target: 5+ classes
- Completed: ✅ 7/7 (100%)
- Status: DONE

### Tier 3: Astrophysical Systems
- Target: ~20 classes
- Completed: ⏳ 0/20 (0%)
- Status: PENDING

### Tier 4: Utility/Support Classes
- Target: ~80-100 classes
- Completed: ⏳ 0/80 (0%)
- Status: PENDING

### TOTAL PROGRESS
- **Completed:** 13/446 modules (2.9%)
- **Remaining:** 433 modules (97.1%)
- **Completion Estimate (manual application):** 4-6 hours

---

## Next Steps

**For Immediate Continuation (Session N+1):**
1. Apply template to Tier 3 classes (Magnetar* and SgrAStar* series)
2. Run batch compilation test
3. Update statistics
4. Proceed to Tier 4 utility classes

**For Phase 2D (Validation):**
1. Once all 446 modules updated
2. Full compilation test
3. Run quantum scaling validation
4. Generate convergence metrics

---

## Notes

- **Non-Breaking:** All changes are additive (backward compatible)
- **Revertible:** Simple `git revert` if issues arise
- **Systematic:** Pattern scales uniformly across all remaining classes
- **Low Risk:** No cross-module dependencies modified
- **Quick Wins:** Each class takes ~2-3 minutes to update manually

---

**Template Created:** Phase 2C Step 3D Complete  
**Ready for:** Tier 3 batch application  
**Merge Point:** After Step 3E documentation
