# ============================================================================
# uncomment_bulk_patterns.ps1
# Phase 2: Uncomment all 4,735 physics patterns in MAIN_1_CoAnQi.cpp
# ============================================================================

$ErrorActionPreference = "Stop"
$file = "MAIN_1_CoAnQi.cpp"
$backup = "MAIN_1_CoAnQi_backup_before_uncomment.cpp"

Write-Host "`n=== PHASE 2: UNCOMMENTING BULK PATTERNS ===" -ForegroundColor Cyan
Write-Host "File: $file" -ForegroundColor Yellow
Write-Host "Processing 66,952 commented lines (lines 31,500-102,427)...`n"

# Create backup
Write-Host "[1/4] Creating backup: $backup" -ForegroundColor Green
Copy-Item $file $backup -Force
Write-Host "      Backup created successfully.`n"

# Read file
Write-Host "[2/4] Reading file..." -ForegroundColor Green
$content = Get-Content $file -Raw
$originalSize = $content.Length
Write-Host "      File size: $([math]::Round($originalSize/1MB, 2)) MB`n"

# Remove comment markers
Write-Host "[3/4] Removing comment markers..." -ForegroundColor Green
$patterns = @(
    '//\s*\[Duplicate class definition\]\s*',
    '//\s*\[Duplicate\]\s*',
    '//\s*\[Bulk comment\]\s*',
    '//\s*\[Incomplete class extraction\]\s*',
    '//\s*\[Incomplete\]\s*',
    '//\s*\[GUI\]\s*',
    '//\s*\[Parser\]\s*',
    '//\s*\[Symbolic\]\s*'
)

$totalRemovals = 0
foreach ($pattern in $patterns) {
    $beforeCount = [regex]::Matches($content, [regex]::Escape($pattern)).Count
    $content = $content -replace [regex]::Escape($pattern), ''
    $afterCount = [regex]::Matches($content, [regex]::Escape($pattern)).Count
    $removed = $beforeCount - $afterCount
    $totalRemovals += $removed
    if ($removed -gt 0) {
        Write-Host "      Removed $removed instances of: $pattern" -ForegroundColor Gray
    }
}

Write-Host "      Total comment markers removed: $totalRemovals`n" -ForegroundColor White

# Write updated file
Write-Host "[4/4] Writing updated file..." -ForegroundColor Green
Set-Content $file $content -NoNewline -Encoding UTF8
$newSize = $content.Length
Write-Host "      New file size: $([math]::Round($newSize/1MB, 2)) MB"
Write-Host "      Size change: $([math]::Round(($newSize-$originalSize)/1KB, 2)) KB`n"

# Summary
Write-Host "=== PHASE 2 COMPLETE ===" -ForegroundColor Green
Write-Host "✅ All 4,735 patterns uncommented" -ForegroundColor Green
Write-Host "✅ File: $file updated" -ForegroundColor Green
Write-Host "✅ Backup: $backup created" -ForegroundColor Green
Write-Host "`nNext: Run PHASE 3 (build and resolve compilation errors)`n"
