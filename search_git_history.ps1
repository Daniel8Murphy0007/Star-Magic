# ==============================================================================
# Git History Search for Registration Counts
# Purpose: Find commits where all 6,597 classes were registered
# Compliance: PowerShell 5.1+ / Git 2.0+
# ==============================================================================

Write-Host "`n=== GIT HISTORY REGISTRATION SEARCH ===" -ForegroundColor Cyan
Write-Host "Searching for commits with complete registration (6,597 classes)`n" -ForegroundColor Yellow

# Function to count registrations in a specific commit
function Count-CommitRegistrations {
    param(
        [string]$CommitHash,
        [string]$File
    )
    
    try {
        $content = git show "${CommitHash}:${File}" 2>$null
        if ($LASTEXITCODE -eq 0) {
            $count = ($content | Select-String "core\.registerPhysicsTerm\(" | Measure-Object).Count
            return $count
        }
    } catch {
        return 0
    }
    return 0
}

# Step 1: Get commits modifying MAIN_1_CoAnQi.cpp in last 60 days
Write-Host "[1/5] Fetching commit history (last 60 days)..." -ForegroundColor Green
$commits = git log --since="60 days ago" --oneline --all -- MAIN_1_CoAnQi.cpp

if (-not $commits) {
    Write-Host "  No commits found for MAIN_1_CoAnQi.cpp in last 60 days" -ForegroundColor Red
    Write-Host "  Expanding search to last 6 months..." -ForegroundColor Yellow
    $commits = git log --since="6 months ago" --oneline --all -- MAIN_1_CoAnQi.cpp
}

$commitList = @()
foreach ($line in $commits) {
    if ($line -match '^([a-f0-9]+)\s+(.+)$') {
        $commitList += @{
            Hash = $matches[1]
            Message = $matches[2]
        }
    }
}

Write-Host "  Found: $($commitList.Count) commits`n" -ForegroundColor White

# Step 2: Find commits before 06:25 AM Nov 23, 2025
Write-Host "[2/5] Filtering commits before 2025-11-23 06:25..." -ForegroundColor Green
$beforeTarget = git log --before="2025-11-23 06:25" --oneline -20

$beforeList = @()
foreach ($line in $beforeTarget) {
    if ($line -match '^([a-f0-9]+)\s+(.+)$') {
        $beforeList += @{
            Hash = $matches[1]
            Message = $matches[2]
        }
    }
}

Write-Host "  Found: $($beforeList.Count) commits before destruction time`n" -ForegroundColor White

# Step 3: Count registrations in each commit
Write-Host "[3/5] Analyzing registration counts (this may take a few minutes)..." -ForegroundColor Green
$results = @()

$i = 0
foreach ($commit in $commitList) {
    $i++
    Write-Progress -Activity "Analyzing commits" -Status "Commit $i of $($commitList.Count)" -PercentComplete (($i / $commitList.Count) * 100)
    
    $mainRegs = Count-CommitRegistrations -CommitHash $commit.Hash -File "MAIN_1_CoAnQi.cpp"
    $wolframRegs = Count-CommitRegistrations -CommitHash $commit.Hash -File "wolfram_physics_classes.cpp"
    
    $results += [PSCustomObject]@{
        CommitHash = $commit.Hash
        Message = $commit.Message
        MAINRegistrations = $mainRegs
        WolframRegistrations = $wolframRegs
        TotalRegistrations = $mainRegs + $wolframRegs
    }
}

Write-Progress -Activity "Analyzing commits" -Completed

# Step 4: Find commits with high registration counts
Write-Host "[4/5] Identifying commits with 5,000+ registrations..." -ForegroundColor Green

$highRegCommits = $results | Where-Object { $_.TotalRegistrations -gt 5000 } | Sort-Object TotalRegistrations -Descending

if ($highRegCommits.Count -eq 0) {
    Write-Host "  ⚠️  No commits found with 5,000+ registrations" -ForegroundColor Red
    Write-Host "  Showing top 10 commits by registration count instead:`n" -ForegroundColor Yellow
    $highRegCommits = $results | Sort-Object TotalRegistrations -Descending | Select-Object -First 10
}

Write-Host "`n=== HIGH REGISTRATION COUNT COMMITS ===" -ForegroundColor Cyan
$highRegCommits | Format-Table CommitHash, MAINRegistrations, WolframRegistrations, TotalRegistrations, @{
    Label = "Message"
    Expression = { $_.Message.Substring(0, [Math]::Min(50, $_.Message.Length)) + "..." }
} -AutoSize

# Step 5: Check current state
Write-Host "`n[5/5] Comparing with current state..." -ForegroundColor Green
$currentMain = (Select-String -Path "MAIN_1_CoAnQi.cpp" -Pattern "core\.registerPhysicsTerm\(" | Measure-Object).Count
$currentWolfram = 5703  # Known from wolfram_physics_classes.cpp
$currentTotal = $currentMain + 0  # Wolfram not executing

Write-Host "`n=== CURRENT STATE ===" -ForegroundColor Cyan
Write-Host "MAIN registrations:     $currentMain" -ForegroundColor Yellow
Write-Host "Wolfram registrations:  0 (function not executing)" -ForegroundColor Red
Write-Host "Total registrations:    $currentMain" -ForegroundColor Yellow
Write-Host "Target registrations:   6,597 (894 MAIN + 5,703 Wolfram)`n" -ForegroundColor Green

# Step 6: Identify best restore candidate
$bestCommit = $highRegCommits | Select-Object -First 1

if ($bestCommit -and $bestCommit.TotalRegistrations -gt $currentTotal) {
    Write-Host "=== RECOMMENDED RESTORE CANDIDATE ===" -ForegroundColor Green
    Write-Host "Commit Hash:            $($bestCommit.CommitHash)" -ForegroundColor Cyan
    Write-Host "Message:                $($bestCommit.Message)" -ForegroundColor White
    Write-Host "MAIN Registrations:     $($bestCommit.MAINRegistrations)" -ForegroundColor Yellow
    Write-Host "Wolfram Registrations:  $($bestCommit.WolframRegistrations)" -ForegroundColor Yellow
    Write-Host "Total Registrations:    $($bestCommit.TotalRegistrations)" -ForegroundColor Green
    Write-Host "Improvement:            +$($bestCommit.TotalRegistrations - $currentTotal) registrations`n" -ForegroundColor Magenta
    
    Write-Host "=== RESTORE COMMANDS ===" -ForegroundColor Yellow
    Write-Host "# Preview the file from this commit:"
    Write-Host "git show $($bestCommit.CommitHash):MAIN_1_CoAnQi.cpp | more`n"
    
    Write-Host "# Extract to recovery file:"
    Write-Host "git show $($bestCommit.CommitHash):MAIN_1_CoAnQi.cpp > MAIN_1_CoAnQi_RECOVERED_$($bestCommit.CommitHash).cpp`n"
    
    Write-Host "# Verify registration count:"
    Write-Host "(Select-String -Path 'MAIN_1_CoAnQi_RECOVERED_$($bestCommit.CommitHash).cpp' -Pattern 'core\.registerPhysicsTerm\(').Count`n"
    
    Write-Host "# If correct, replace current file:"
    Write-Host "Copy-Item MAIN_1_CoAnQi.cpp MAIN_1_CoAnQi_BACKUP_$(Get-Date -Format 'yyyyMMdd_HHmmss').cpp"
    Write-Host "Copy-Item MAIN_1_CoAnQi_RECOVERED_$($bestCommit.CommitHash).cpp MAIN_1_CoAnQi.cpp`n"
} else {
    Write-Host "⚠️  No better commit found in history" -ForegroundColor Red
    Write-Host "Current state may already be the best available" -ForegroundColor Yellow
    Write-Host "Recommendation: Generate missing registrations instead of restoring`n" -ForegroundColor Cyan
}

# Save full report
$reportFile = "git_history_report_$(Get-Date -Format 'yyyyMMdd_HHmmss').csv"
$results | Sort-Object TotalRegistrations -Descending | 
    Export-Csv $reportFile -NoTypeInformation

Write-Host "📊 Full report saved to: $reportFile`n" -ForegroundColor Cyan

# Display statistics
Write-Host "=== STATISTICS ===" -ForegroundColor Cyan
$avgRegs = ($results | Measure-Object -Property TotalRegistrations -Average).Average
$maxRegs = ($results | Measure-Object -Property TotalRegistrations -Maximum).Maximum
$minRegs = ($results | Measure-Object -Property TotalRegistrations -Minimum).Minimum

Write-Host "Commits analyzed:       $($results.Count)" -ForegroundColor White
Write-Host "Average registrations:  $([int]$avgRegs)" -ForegroundColor White
Write-Host "Maximum registrations:  $maxRegs" -ForegroundColor Green
Write-Host "Minimum registrations:  $minRegs" -ForegroundColor Red
Write-Host "Current registrations:  $currentMain`n" -ForegroundColor Yellow

Write-Host "✅ COMPLETE - Git history search finished" -ForegroundColor Green
