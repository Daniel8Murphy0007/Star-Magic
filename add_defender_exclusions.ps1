# Star-Magic Workspace — Windows Defender Exclusions
# Run this script as Administrator (UAC prompt will appear)
# Purpose: prevent MsMpEng from scanning every file VS Code opens/writes,
#          which causes CPU spikes and the "nothing responds" freezes.

$WorkspaceRoot = "C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
$VenvPath      = "$WorkspaceRoot\.venv_py314_backup"
$BuildPath     = "$WorkspaceRoot\build_msvc"

Write-Host "`n=== Adding Windows Defender Exclusions for Star-Magic ===" -ForegroundColor Cyan

# --- Path exclusions (entire directories, not scanned at all) ---
$PathExclusions = @(
    $WorkspaceRoot,     # entire repo — largest single win
    $VenvPath,          # Python venv (thousands of .pyd/.dll files)
    $BuildPath,         # MSVC build output (obj, intermediate files)
    "$WorkspaceRoot\build",
    "C:\Users\tmsjd\AppData\Roaming\Code",   # VS Code user data
    "C:\Users\tmsjd\.vscode"                  # VS Code extension data
)

foreach ($path in $PathExclusions) {
    if (Test-Path $path) {
        Add-MpPreference -ExclusionPath $path
        Write-Host "  [PATH] $path" -ForegroundColor Green
    } else {
        Write-Host "  [SKIP] $path (not found)" -ForegroundColor Yellow
    }
}

# --- Process exclusions (these processes' file I/O is never scanned) ---
$ProcessExclusions = @(
    "Code.exe",              # VS Code renderer + main
    "node.exe",              # Pylance, Copilot, extension host
    "python.exe",            # Python interpreter
    "python3.14.exe",        # Explicit version
    "cl.exe",                # MSVC compiler
    "link.exe",              # MSVC linker
    "msbuild.exe",           # Build system
    "cmake.exe",             # CMake
    "git.exe"                # Git (triggers scan on every file it touches)
)

foreach ($proc in $ProcessExclusions) {
    Add-MpPreference -ExclusionProcess $proc
    Write-Host "  [PROC] $proc" -ForegroundColor Green
}

# --- Extension exclusions (file types never scanned) ---
$ExtExclusions = @(
    ".obj",    # MSVC compiler output
    ".pdb",    # Debug symbols
    ".ilk",    # Incremental linker
    ".tlog",   # MSBuild tracking logs
    ".pyc",    # Python bytecode
    ".pyd"     # Python extension modules
)

foreach ($ext in $ExtExclusions) {
    Add-MpPreference -ExclusionExtension $ext
    Write-Host "  [EXT]  $ext" -ForegroundColor Green
}

Write-Host "`n=== Done! Exclusions active immediately — no reboot required ===" -ForegroundColor Cyan
Write-Host "    MsMpEng will no longer scan VS Code, Pylance, or Git file operations." -ForegroundColor White
Write-Host "    Expected: CPU spikes from MsMpEng drop to near zero during builds/edits.`n" -ForegroundColor White
