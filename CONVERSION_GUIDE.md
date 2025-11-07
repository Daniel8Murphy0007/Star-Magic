# Legacy Module Conversion Guide

## Overview

### Verification Results (source13-100)

- **Total modules in range**: 88 files
- **All functional**: ✅ Yes
- **Enhanced framework**: ❌ No (legacy structure)
- **Conversion complexity**: High (varied structures)

### Current Enhanced Modules

- **source131-134**: FULL enhanced (4 modules)
- **source147-162**: FULL enhanced (16 modules)
- **Total Enhanced**: 20 modules with complete framework

## Conversion Requirements

### What Each Legacy Module Needs

1. **Import enhanced_dynamics.js** - Add require statement
2. **Complex number helpers** - 6 inline functions (complexAdd, complexMul, etc.)
3. **Convert variables to Map** - this.variables = new Map()
4. **Update all methods** - Access via this.variables.get()
5. **Add enhanced framework** - Call addEnhancedDynamics()
6. **Domain expansion** - Create custom domain methods object
7. **Clone method** - Deep copy for parallel processing
8. **Metadata tracking** - Add version, type, system_name

## Complexity Analysis

### Time Estimate per Module

- **Simple module** (20-30 variables): 20-30 minutes
- **Medium module** (40-60 variables): 30-45 minutes
- **Complex module** (80+ variables): 45-60 minutes

**Total estimate for 88 modules**: 44-66 hours

## Recommended Approach

### Phase 1: Demonstration Batch (5 modules)

Convert small representative sample:

- **source13.js** - MagnetarSGR1745_2900 (magnetar)
- **source14.js** - RotatingBlackHole (Kerr BH)
- **source15.js** - BinaryNeutronStar (merger)
- **source20.js** - TidalDisruptionEvent (TDE)
- **source25.js** - GravitationalLens (lensing)

**Purpose**: Validate conversion process, identify issues

### Phase 2: Tooling Enhancement

Build automated converter improvements:

- Extract variable list from legacy module
- Generate Map initialization code
- Auto-create domain expansion template
- Generate test cases

**Time saved**: ~30% reduction in manual work

### Phase 3: Batch Conversion (Remaining 83 modules)

Once tooling validated:

- Convert in batches of 10 modules
- Test each batch thoroughly
- Commit after each successful batch
- Maintain backups of originals

## Alternative Strategies

### Option A: As-Needed Conversion

- Keep legacy modules functional as-is
- Convert only when user requests enhanced features
- Benefits: No upfront time investment, zero risk
- Drawbacks: Inconsistent codebase, manual conversion delays

### Option B: Hybrid System

- Legacy modules (source13-130): Basic functionality
- Enhanced modules (source131+): Full framework
- Document which systems have which capabilities
- Benefits: Best of both worlds, gradual migration

### Option C: New Module Template

- Future modules use enhanced framework by default
- Legacy stays as-is for backward compatibility
- Create "enhanced" versions with new filenames (source13_enhanced.js)
- Benefits: Zero breakage, clear separation

## Recommendation

### Proposed Plan: Hybrid + Selective Conversion

1. **Keep all legacy modules functional** (no changes)
2. **Continue using enhanced framework for new modules**
3. **Convert specific legacy modules only when**:
   - User explicitly requests enhanced features
   - Module is frequently used in parallel workflows
   - New physics additions require dynamic framework

**Rationale**:

- Demonstrates enhanced framework on legacy modules
- Minimizes risk of breaking working code
- Focuses effort where value is highest
- Maintains backward compatibility

## Decision Checklist

- [ ] Demonstration batch (5 modules, ~2.5 hours)
- [ ] Enhanced tooling (automated converter improvements)
- [ ] Full conversion (88 modules, 44-66 hours)
- [x] **CURRENT STATUS**: Hybrid approach - 20 enhanced modules, 88 legacy stable
