# Commit Plan for 41 Unstaged Files
**Date:** November 20, 2025 @ 4:00 AM  
**Analyzed Files:** 41 (5 modified, 1 deleted, 35 untracked new)

---

## ✅ COMMIT 1: Production Code Updates (RECOMMENDED - SAFE TO COMMIT)

### Files to Include (4 files)
```bash
git add CMakeLists.txt
git add source176_auto_full_uqff.cpp
git add Source48-5.cpp
git add INTEGRATION_TRACKER.csv
```

### Changes Summary
1. **CMakeLists.txt** - Added PhysicsTermExtractor executable with Wolfram WSTP support
2. **source176_auto_full_uqff.cpp** - Bulletproof Wolfram export automation (FINAL VERSION)
3. **Source48-5.cpp** - StarMagicUQFFModule refactored to inline header-only pattern
4. **INTEGRATION_TRACKER.csv** - Added source176 entry (integrated 2025-11-19)

### Commit Command
```bash
git commit -m "Add PhysicsTermExtractor + bulletproof Wolfram export + StarMagicUQFFModule refactoring

- CMakeLists.txt: New PhysicsTermExtractor target with WSTP
- source176_auto_full_uqff.cpp: Final bulletproof version with runtime path detection
- Source48-5.cpp: Convert to inline header-only class definition
- INTEGRATION_TRACKER.csv: Track source176 integration"
```

---

## ✅ COMMIT 2: Documentation Updates (RECOMMENDED - AFTER REVIEW)

### Files to Include (13 files)
```bash
# Updated/corrected documentation
git add SESSION_SUMMARY_CURRENT.md
git add BATCH_6_COMPLETE.md
git add BATCH_EXTRACTION_SUMMARY.md
git add COMPLETE_PHYSICS_INVENTORY.md
git add EXTRACTION_COMPLETE_REPORT.md
git add PHYSICS_TERMS_COMPREHENSIVE_REPORT.md
git add PHYSICS_CLASS_ANALYSIS_SUMMARY.md

# New analysis/tracking files
git add UNSTAGED_CHANGES_ANALYSIS.md
git add COMMIT_PLAN.md
git add QUICK_CLASS_LOOKUP.txt
git add PHYSICS_CLASSES_BY_FILE.txt
```

### Changes Summary
- Updated all session/batch reports with current state notices (807 terms, 26.9%, Batch 15)
- Added historical context markers to inventory documents
- Created comprehensive unstaged changes analysis
- Added commit planning guide (this file)
- Preserved historical extraction data for reference

### Commit Command
```bash
git commit -m "Update documentation: Sync with 807-term state, add historical context

Session/Batch Reports:
- SESSION_SUMMARY_CURRENT.md: Updated to 807 terms (26.9%)
- BATCH_6_COMPLETE.md: Added current status context
- BATCH_EXTRACTION_SUMMARY.md: Updated metrics

Inventory Documents:
- COMPLETE_PHYSICS_INVENTORY.md: Added current integration notice
- EXTRACTION_COMPLETE_REPORT.md: Historical extraction context
- PHYSICS_TERMS_COMPREHENSIVE_REPORT.md: 1505 term inventory with status

Analysis:
- UNSTAGED_CHANGES_ANALYSIS.md: Complete 41-file analysis
- COMMIT_PLAN.md: Commit strategy guide"
```

---

## ✅ COMMIT 3: Data Inventories (OPTIONAL - AFTER VERIFICATION)

### Files to Include (4 CSV files)
```bash
git add COMPLETE_PHYSICS_CLASS_INVENTORY.csv
git add all_physics_classes.csv
git add all_physics_terms_inventory.csv
git add missing_classes_extraction_list.csv
```

### Changes Summary
- **COMPLETE_PHYSICS_CLASS_INVENTORY.csv** - 473 classes catalogued (SourceFile, ClassName, LineNumber, BaseClass, PhysicsType, MethodSignature)
- **all_physics_classes.csv** - Comprehensive class export
- **all_physics_terms_inventory.csv** - Complete physics terms catalog
- **missing_classes_extraction_list.csv** - Gap analysis

### ⚠️ Verification Required
Before committing, verify CSV data matches current MAIN_1_CoAnQi.cpp state:
- Current: 781 PhysicsTerm classes
- CSV reports: 473 classes
- **Action needed:** Re-run inventory extraction on current MAIN_1_CoAnQi.cpp or accept as historical snapshot

### Commit Command
```bash
git commit -m "Add physics class/term inventory data (historical extraction snapshot)

CSV Inventories:
- COMPLETE_PHYSICS_CLASS_INVENTORY.csv: 473 classes catalogued
- all_physics_classes.csv: Comprehensive class export
- all_physics_terms_inventory.csv: Complete terms catalog
- missing_classes_extraction_list.csv: Gap analysis

Note: This data represents extraction phase snapshot. Current active integration
has 781 classes (807 terms) in MAIN_1_CoAnQi.cpp. See UNSTAGED_CHANGES_ANALYSIS.md."
```

---

## 🔴 SKIP: Auto-Generated Code (DO NOT COMMIT)

### Files to Skip (12 C++ generated files)
These are auto-generated artifacts that should be in .gitignore:

```
register_1505_terms.cpp              # AUTO-GENERATED skeleton (TODOs only)
register_all_extracted_classes.cpp   # Empty placeholder (0 classes)
physics_registration.cpp             # Auto-generated registrations
extract_all_physics_terms.cpp        # Extraction output
physics_term_extractor_main.cpp      # Wolfram-powered extractor
generated_500_physics_terms.cpp      # Generated implementations
physics_term_wrappers_all.cpp        # All term wrappers
physics_term_wrappers_top50.cpp      # Top 50 wrappers
extracted_complete_implementations.cpp  # Complete extracted code
extracted_all_missing_classes.h      # Missing class headers
extracted_complete_classes.h         # Complete class headers
```

### Add to .gitignore
```bash
# Add these patterns to .gitignore
register_*_terms.cpp
physics_registration.cpp
generated_*_physics_terms.cpp
physics_term_wrappers_*.cpp
extracted_*.cpp
extracted_*.h
generated_wrappers.txt
physics_term_inventory.txt
```

---

## 🟡 REVIEW: Extraction Scripts (DECISION NEEDED)

### Files to Review (9 scripts)
```
batch_extract_all_classes.ps1         # Batch class extraction
extract_all_500_physics_terms.ps1     # Target 500+ terms
extract_complete_classes.ps1          # Complete extraction pipeline
generate_missing_physics_terms.ps1    # Missing term wrappers
generate_physics_wrappers_simple.ps1  # Simple wrapper generation
extract_all_physics_classes.py        # Python class extractor
create_lookup.py                      # Quick lookup generator
register_all_physics_terms.py         # Auto-generate registrations
wolfram terms search.py               # Wolfram term injection
```

### Options

#### Option A: Keep Scripts (Recommended for Reproducibility)
**Commit if:**
- Scripts may be reused for future extractions
- Useful for documentation/reproducibility
- Team members need extraction tools

```bash
git add batch_extract*.ps1 extract_*.ps1 generate_*.ps1
git add extract_all_physics_classes.py create_lookup.py register_all_physics_terms.py
git add "wolfram terms search.py"

git commit -m "Add extraction/generation scripts for physics term automation

PowerShell:
- batch_extract_all_classes.ps1: Batch automation
- extract_all_500_physics_terms.ps1: Target-based extraction
- extract_complete_classes.ps1: Complete pipeline
- generate_missing_physics_terms.ps1: Gap filling
- generate_physics_wrappers_simple.ps1: Wrapper generation

Python:
- extract_all_physics_classes.py: Class extraction
- create_lookup.py: Quick reference generation
- register_all_physics_terms.py: Registration automation
- wolfram terms search.py: Term injection tool"
```

#### Option B: Skip Scripts (Cleaner Repo)
**Skip if:**
- One-time usage only
- Don't want to maintain scripts
- Prefer manual workflow

```bash
# Add to .gitignore
*.ps1  # or list specific scripts
extract_*.py
generate_*.py
create_lookup.py
register_all_physics_terms.py
```

---

## 🔴 DELETE: Temporary Files (CLEAN UP)

### Files to Delete
```bash
git rm --cached Untitled-1.py  # Temporary Python script
# Review Untitled-2.cpp before deleting (not in unstaged list?)
```

### Delete from Filesystem
```bash
Remove-Item Untitled-1.py
Remove-Item physics_term_inventory.txt
Remove-Item generated_wrappers.txt
```

---

## 🗑️ HANDLE DELETION: INTEGRATION_PROGRESS_REPORT.md

### File Deleted
- **INTEGRATION_PROGRESS_REPORT.md** - Removed (likely superseded by newer tracking)

### Action
```bash
git rm INTEGRATION_PROGRESS_REPORT.md  # Stage deletion
```

Or restore if needed:
```bash
git restore INTEGRATION_PROGRESS_REPORT.md  # Restore from last commit
```

---

## 📋 Recommended Execution Order

### Phase 1: Production Code (Immediate)
```bash
# Commit 1: Production code changes
git add CMakeLists.txt source176_auto_full_uqff.cpp Source48-5.cpp INTEGRATION_TRACKER.csv
git commit -m "Add PhysicsTermExtractor + bulletproof Wolfram export + StarMagicUQFFModule refactoring"
git push origin master
```

### Phase 2: Documentation (After Review)
```bash
# Verify updates are correct, then commit
git add SESSION_SUMMARY_CURRENT.md BATCH_6_COMPLETE.md BATCH_EXTRACTION_SUMMARY.md
git add COMPLETE_PHYSICS_INVENTORY.md EXTRACTION_COMPLETE_REPORT.md PHYSICS_TERMS_COMPREHENSIVE_REPORT.md
git add PHYSICS_CLASS_ANALYSIS_SUMMARY.md UNSTAGED_CHANGES_ANALYSIS.md COMMIT_PLAN.md
git add QUICK_CLASS_LOOKUP.txt PHYSICS_CLASSES_BY_FILE.txt

git commit -m "Update documentation: Sync with 807-term state, add historical context"
git push origin master
```

### Phase 3: Cleanup (Before Commits 1-2)
```bash
# Delete temporary files
Remove-Item Untitled-1.py -ErrorAction SilentlyContinue
git rm --cached Untitled-1.py -ErrorAction SilentlyContinue

# Stage deleted file
git rm INTEGRATION_PROGRESS_REPORT.md

# Update .gitignore
Add-Content .gitignore "`n# Auto-generated files"
Add-Content .gitignore "register_*_terms.cpp"
Add-Content .gitignore "physics_registration.cpp"
Add-Content .gitignore "generated_*_physics_terms.cpp"
Add-Content .gitignore "physics_term_wrappers_*.cpp"
Add-Content .gitignore "extracted_*.cpp"
Add-Content .gitignore "extracted_*.h"
Add-Content .gitignore "generated_wrappers.txt"
Add-Content .gitignore "physics_term_inventory.txt"
Add-Content .gitignore "Untitled-*.py"
Add-Content .gitignore "Untitled-*.cpp"

git add .gitignore
git commit -m "Update .gitignore: Exclude auto-generated files and temporary scripts"
```

### Phase 4: Optional - CSV Data
```bash
# Only if verified against current state
git add COMPLETE_PHYSICS_CLASS_INVENTORY.csv all_physics_classes.csv
git add all_physics_terms_inventory.csv missing_classes_extraction_list.csv
git commit -m "Add physics class/term inventory data (historical snapshot)"
git push origin master
```

### Phase 5: Optional - Scripts
```bash
# Only if keeping for future use
git add batch_extract*.ps1 extract_*.ps1 generate_*.ps1
git add *.py
git commit -m "Add extraction/generation automation scripts"
git push origin master
```

---

## 🎯 Recommended Minimal Commit Strategy

**If you want to commit only essential changes:**

```bash
# Essential production code
git add CMakeLists.txt source176_auto_full_uqff.cpp Source48-5.cpp INTEGRATION_TRACKER.csv
git rm INTEGRATION_PROGRESS_REPORT.md

# Essential documentation updates
git add SESSION_SUMMARY_CURRENT.md BATCH_6_COMPLETE.md BATCH_EXTRACTION_SUMMARY.md
git add UNSTAGED_CHANGES_ANALYSIS.md

# Update .gitignore
# (Add auto-generated file patterns manually)
git add .gitignore

# Commit all at once
git commit -m "Production updates: PhysicsTermExtractor, Wolfram export, docs sync

Code Changes:
- CMakeLists.txt: Add PhysicsTermExtractor target
- source176_auto_full_uqff.cpp: Bulletproof Wolfram export
- Source48-5.cpp: Header-only refactoring
- INTEGRATION_TRACKER.csv: Track source176

Documentation:
- SESSION_SUMMARY_CURRENT.md: Update to 807 terms (26.9%)
- BATCH_6_COMPLETE.md: Add current context
- BATCH_EXTRACTION_SUMMARY.md: Update metrics
- UNSTAGED_CHANGES_ANALYSIS.md: Complete 41-file analysis

Cleanup:
- Remove INTEGRATION_PROGRESS_REPORT.md (superseded)
- Update .gitignore for auto-generated files"

git push origin master
```

---

## Summary Statistics

### What to Commit (Recommended)
- ✅ **Production code:** 4 files (CMakeLists, source176, Source48-5, INTEGRATION_TRACKER.csv)
- ✅ **Documentation:** 11 files (session reports, inventories, analysis)
- 🟡 **CSV data:** 4 files (optional, verify first)
- 🟡 **Scripts:** 9 files (optional, decision needed)
- 🔴 **Auto-generated:** 0 files (skip, add to .gitignore)
- 🗑️ **Cleanup:** 1 deletion (INTEGRATION_PROGRESS_REPORT.md)

### Total Recommended Commits: 1-3
1. **Production Code** (essential) - 4 files
2. **Documentation** (recommended) - 11 files
3. **Optional Data/Scripts** (as needed) - 13 files

### Files to Ignore: 12-13
- All auto-generated C++ code
- Temporary files (Untitled-*.py)
- Generated text outputs

---

**Status:** Ready for user decision on commit strategy.
**Next Action:** Choose minimal (1 commit) or comprehensive (2-3 commits) strategy.
