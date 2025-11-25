# Inbox Dropzone

**Purpose:** Staging area for new physics modules before integration into MAIN_1_CoAnQi.cpp

## Workflow (Phase 0 - PLAN.md)

1. **Drop new source files here** (source*.cpp, new physics modules)
2. **Automated scanning** will detect new files
3. **Validation** checks for:
   - Syntax correctness
   - Physics term completeness
   - Integration compatibility
4. **Integration** into MAIN_1_CoAnQi.cpp via SOURCE blocks
5. **Cleanup** - Files moved to archive after successful integration

## Current Status

- **492 physics terms** integrated (INTEGRATION_TRACKER.csv)
- **446 modules** (SOURCE1-116) active in MAIN_1_CoAnQi.cpp
- **Ready for new modules**

## Usage

Simply copy new physics module files into this folder. The automated workflow will:
- Detect the files
- Validate physics terms
- Update INTEGRATION_TRACKER.csv
- Integrate into main codebase
- Archive the source files

---
**Created:** November 23, 2025  
**Phase:** 0 (PLAN.md workflow preparation)
