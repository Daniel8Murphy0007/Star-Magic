# commit_rounds_19_30.ps1 - per-round commits for forensic rounds 19-30
# Pauses between each commit so you can review. Ctrl+C aborts.
$ErrorActionPreference = 'Stop'
Set-Location 'C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic'

function Step {
    param([string]$Number, [string]$Title, [string[]]$Files)
    Write-Host ""
    Write-Host "===================================================="
    Write-Host "Round $Number - $Title"
    Write-Host "----------------------------------------------------"
    Write-Host "Files to stage:"
    foreach ($f in $Files) { Write-Host "  $f" }
    Write-Host "Press Enter to stage + commit (Ctrl+C aborts)..."
    Read-Host | Out-Null
    foreach ($f in $Files) { git add -- $f }
    git commit -m "Round $Number forensic recovery: $Title"
    Write-Host "  --> Round $Number committed."
}

Step "19" "4684f438 LaTeX-3 mangled identifier reversal (24898 fixes)" @(
    "whitepapers/*.md",
    "whitepapers/*.tex"
)

Step "20" "20bb5e4f restore 13 SCm pedagogical sections + utf-8 reconfigure" @(
    "dpm_vacuum_manifold.py",
    "scm_vacuum_manifold.py"
)

Step "23" "a009fadc surgical mojibake fix in CP2 SCm-class range" @(
    "CondensedPhysics2.py"
)

Step "24" "b6756320 restore SCm Super-K comparison + oscillation plot block" @(
    "dpm_vacuum_manifold.py",
    "scm_vacuum_manifold.py"
)

Step "27" "81bebeb6 restore 6 pdf-only pedagogical sections" @(
    "dpm_vacuum_manifold.py",
    "scm_vacuum_manifold.py"
)

Step "29" "9c1c7083 restore scm_latex_exporter.py utility" @(
    "scm_latex_exporter.py"
)

Step "31" "SESSION_LOG entries for rounds 19-30 forensic walk" @(
    "SESSION_LOG.md"
)

Write-Host ""
Write-Host "===================================================="
Write-Host "All 7 commits done. Press Enter to push origin master (Ctrl+C skips)..."
Read-Host | Out-Null
git push origin master
Write-Host "Pushed."
