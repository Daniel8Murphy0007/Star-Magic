# Workspace Update Checklist
## Files to Update on EVERY Commit

### 🔴 CRITICAL - Update Every Commit

1. **`Star-Magic.md`** (Workspace Status Snapshot)
   - Current line count of MAIN_1_CoAnQi.cpp
   - Latest commit hash and date
   - C++ standard version
   - MSVC compiler version
   - Executable size (MB)
   - Module count (SOURCE1-116)
   - Physics class count
   - Wolfram companion file count
   - Current project phase
   - Build status (✅/❌)
   
2. **`MAIN_1_CoAnQi_integration_status.json`**
   - `last_updated` timestamp
   - `file_statistics.total_lines`
   - `file_statistics.last_modified`
   - `git_status.last_commit`
   - `git_status.commit_date`
   - `git_status.commit_message`
   - `registry_status.registered` (class count)
   - `build_status` field
   - `project_phase` field

3. **`README.md`**
   - Header timestamp (Last Updated)
   - Latest commit reference
   - Progress metrics (physics classes, Wolfram files)
   - Build status line
   - Current phase description

### 🟡 MODERATE - Update on Major Changes

4. **`PLAN.md`**
   - Update when phases change
   - Update when roadmap milestones complete
   - Update when Phase 2-6 architecture decisions made

5. **`INTEGRATION_TRACKER.csv`**
   - Update when new source files integrated
   - Update when module counts change
   - Update when SOURCE blocks added/modified

6. **`BUILD_INSTRUCTIONS_PERMANENT.md`**
   - Update when CMake configuration changes
   - Update when compiler version changes
   - Update when build dependencies change

### 🟢 OCCASIONAL - Update on Specific Events

7. **`COMPLETE_PHYSICS_CLASS_INVENTORY.csv`**
   - Update when new physics classes added (source77-81 batches)
   - Update after major integration sessions

8. **`RESTART_STATUS.md`**
   - Update at milestone checkpoints
   - Update after major architectural changes
   - Update when starting new development sessions

9. **`ENHANCEMENT_GUIDE.md`**
   - Update when 2.0-Enhanced framework patterns change
   - Update when new self-expanding capabilities added

10. **`CMakeLists.txt`**
    - Update when new source files added
    - Update when dependencies change
    - Update when build flags modified

---

## Automated Update Template for Star-Magic.md

```markdown
# Star-Magic.md

## INTEGRATION STATUS ({{ DATE }} @ {{ TIME }})

**Platform:** MAIN_1_CoAnQi.cpp ({{ LINE_COUNT }} lines) - Conscious Quantum Intelligence UQFF Calculator  
**Modules:** {{ MODULE_COUNT }} integrated physics terms across SOURCE1-116 blocks  
**Build:** CMake {{ CMAKE_VERSION }} + MSVC {{ MSVC_VERSION }}, C++{{ CPP_STANDARD }} standard  
**Executable:** MAIN_1_CoAnQi.exe ({{ SIZE_MB }} MB)  
**Commit:** {{ COMMIT_HASH }} (master branch)  
**Framework:** 2.0-Enhanced Self-Expanding with dynamic term registration  
**Repository:** Daniel8Murphy0007/Star-Magic

---

## Current Workspace Snapshot

- **Total Lines:** {{ TOTAL_LINES }}
- **Active Code:** {{ ACTIVE_CODE_LINES }}
- **Physics Classes:** {{ PHYSICS_CLASS_COUNT }} ({{ UQFF_COUNT }} UQFF + {{ WOLFRAM_COUNT }} Wolfram)
- **Wolfram Companions:** {{ WOLFRAM_FILE_COUNT }} files (source{{ WOLFRAM_RANGE }})
- **Build Status:** {{ BUILD_STATUS }}
- **Last Modified:** {{ LAST_MODIFIED }}
- **Current Phase:** {{ PROJECT_PHASE }}

---

## Quick Stats

| Metric | Value |
|--------|-------|
| C++ Standard | C++{{ CPP_STANDARD }} |
| Compiler | MSVC {{ MSVC_VERSION }} |
| Executable Size | {{ EXE_SIZE }} MB |
| Compression | UPX {{ UPX_VERSION }} ({{ COMPRESSION_PCT }}% reduction) |
| SOURCE Modules | {{ SOURCE_MODULE_COUNT }} |
| Astronomical Systems | {{ ASTRO_SYSTEM_COUNT }} |
| Wolfram Integration | ✅ WSTP {{ WSTP_VERSION }} |

---

## Build Command

```powershell
# Configure
cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64

# Build
cmake --build build_msvc --config Release --target MAIN_1_CoAnQi

# Run
.\build_msvc\Release\MAIN_1_CoAnQi.exe
```

---

**For complete theoretical framework, see `UQFF_THEORY.md`**  
**For construction phase roadmap, see `PLAN.md`**  
**For workspace status history, see `UQFF_VALIDATION_CONVERSATION_CAPTURE.md`**
```

---

## PowerShell Script to Extract Values

Save this as `update_workspace_status.ps1`:

```powershell
# Get current values
$lineCount = (Get-Content "MAIN_1_CoAnQi.cpp" -Raw).Split("`n").Count
$commitHash = (git rev-parse --short HEAD)
$commitDate = (git log -1 --format=%cd --date=short)
$fileSize = [math]::Round((Get-Item "MAIN_1_CoAnQi.cpp").Length / 1MB, 2)
$timestamp = Get-Date -Format "MMMM dd, yyyy @ h:mm tt"

# Update Star-Magic.md with current values
Write-Host "Current Workspace Status:"
Write-Host "  Lines: $lineCount"
Write-Host "  Commit: $commitHash"
Write-Host "  Date: $commitDate"
Write-Host "  Size: $fileSize MB"
Write-Host "  Timestamp: $timestamp"
```

---

## Git Pre-Commit Hook (Optional Automation)

Create `.git/hooks/pre-commit`:

```bash
#!/bin/bash
# Auto-update workspace status files before commit

LINE_COUNT=$(wc -l < MAIN_1_CoAnQi.cpp)
COMMIT_HASH=$(git rev-parse --short HEAD)
TIMESTAMP=$(date +"%B %d, %Y @ %I:%M %p")

# Update Star-Magic.md header
sed -i "s/^## INTEGRATION STATUS (.*/## INTEGRATION STATUS ($TIMESTAMP)/" Star-Magic.md
sed -i "s/^\\*\\*Platform:\\*\\* MAIN_1_CoAnQi.cpp (.* lines)/**Platform:** MAIN_1_CoAnQi.cpp ($LINE_COUNT lines)/" Star-Magic.md

# Stage the updated file
git add Star-Magic.md
```

---

## Recommended Workflow

### On Every Commit:
1. Run build: `cmake --build build_msvc --config Release`
2. Get stats: `.\update_workspace_status.ps1`
3. Update `Star-Magic.md` with current values
4. Update `MAIN_1_CoAnQi_integration_status.json` last_updated
5. Update `README.md` header timestamp
6. Commit: `git commit -m "Your message"`
7. Push: `git push origin master`

### Weekly/Monthly:
- Review and update `PLAN.md` roadmap
- Update `INTEGRATION_TRACKER.csv` with new modules
- Snapshot to `WORKSPACE_STATUS_{{DATE}}.md` for historical record

---

**Last Updated:** November 26, 2025
