# ============================================================================
# uncomment_safe_patterns.ps1
# Safely uncomment net-new physics patterns from MAIN_1_CoAnQi.cpp
# ============================================================================

param(
    [string]$ReportFile = "deduplication_analysis_20251122_152350.json",
    [string]$SourceFile = "MAIN_1_CoAnQi.cpp",
    [int]$BatchNumber = 1,
    [int]$BatchSize = 50
)

$ErrorActionPreference = "Stop"

Write-Host "`n=== SAFE PATTERN UNCOMMENTER ===" -ForegroundColor Cyan
Write-Host "Report: $ReportFile" -ForegroundColor Yellow
Write-Host "Source: $SourceFile" -ForegroundColor Yellow
Write-Host "Batch: $BatchNumber (size=$BatchSize)`n" -ForegroundColor Yellow

# Load analysis report
if (-not (Test-Path $ReportFile)) {
    Write-Host "✗ Report file not found: $ReportFile" -ForegroundColor Red
    exit 1
}

$report = Get-Content $ReportFile | ConvertFrom-Json
$allPatterns = $report.patterns

# Filter net-new patterns
$netNewPatterns = $allPatterns | Where-Object { 
    $_.is_commented -and $_.category -eq 'net_new' 
}

Write-Host "[1/6] Found $($netNewPatterns.Count) net-new patterns" -ForegroundColor Green

# Get batch
$startIndex = ($BatchNumber - 1) * $BatchSize
$endIndex = [Math]::Min($startIndex + $BatchSize, $netNewPatterns.Count)
$batchPatterns = $netNewPatterns[$startIndex..($endIndex-1)]

Write-Host "[2/6] Batch $BatchNumber contains $($batchPatterns.Count) patterns (lines $startIndex-$endIndex)" -ForegroundColor Green

# Create backup
$timestamp = Get-Date -Format "yyyyMMdd_HHmmss"
$backupFile = "MAIN_1_CoAnQi_backup_batch${BatchNumber}_${timestamp}.cpp"
Copy-Item $SourceFile $backupFile -Force
Write-Host "[3/6] Backup created: $backupFile" -ForegroundColor Green

# Load source file
$content = Get-Content $SourceFile

# Process each pattern in batch
$uncommentsApplied = 0
$lineChanges = @{}

foreach ($pattern in $batchPatterns) {
    # Parse line range (format: "start-end")
    $lineRange = $pattern.lines -split '-'
    $startLine = [int]$lineRange[0] - 1  # Convert to 0-indexed
    $endLine = [int]$lineRange[1] - 1
    
    Write-Host "  Uncommenting: $($pattern.name) (lines $($startLine+1)-$($endLine+1))" -ForegroundColor Gray
    
    # Uncomment lines in range
    for ($i = $startLine; $i -le $endLine; $i++) {
        $originalLine = $content[$i]
        $cleanedLine = $originalLine
        
        # Remove // [Duplicate] markers
        $cleanedLine = $cleanedLine -replace '^\s*//\s*\[Duplicate[^\]]*\]\s*', ''
        
        # Remove // comment markers (preserve indentation)
        $cleanedLine = $cleanedLine -replace '^(\s*)//\s*', '$1'
        
        if ($cleanedLine -ne $originalLine) {
            $content[$i] = $cleanedLine
            $uncommentsApplied++
        }
    }
}

Write-Host "[4/6] Applied $uncommentsApplied line uncommits" -ForegroundColor Green

# Save modified file
$content | Set-Content $SourceFile -Encoding UTF8
Write-Host "[5/6] Saved: $SourceFile" -ForegroundColor Green

# Verify compilation
Write-Host "[6/6] Compiling..." -ForegroundColor Yellow
$buildResult = & cmake --build build_msvc --config Release --target MAIN_1_CoAnQi 2>&1

if ($LASTEXITCODE -eq 0) {
    Write-Host "✓ COMPILATION SUCCESSFUL" -ForegroundColor Green
    Write-Host "`n=== BATCH $BatchNumber COMPLETE ===" -ForegroundColor Green
    Write-Host "Uncommented: $($batchPatterns.Count) patterns" -ForegroundColor Green
    Write-Host "Line changes: $uncommentsApplied" -ForegroundColor Green
    Write-Host "Backup: $backupFile`n" -ForegroundColor Green
    
    # Update todo for next batch
    if ($endIndex -lt $netNewPatterns.Count) {
        $remaining = $netNewPatterns.Count - $endIndex
        $nextBatch = $BatchNumber + 1
        Write-Host "Next: Run with -BatchNumber $nextBatch" -ForegroundColor Cyan
        Write-Host "Remaining patterns: $remaining" -ForegroundColor Cyan
    } else {
        Write-Host "All batches complete! Ready for registration phase." -ForegroundColor Cyan
    }
    
    exit 0
} else {
    Write-Host "✗ COMPILATION FAILED" -ForegroundColor Red
    Write-Host "Errors:" -ForegroundColor Red
    Write-Host $buildResult -ForegroundColor Red
    
    # Rollback
    Write-Host "`nRolling back..." -ForegroundColor Yellow
    Copy-Item $backupFile $SourceFile -Force
    Write-Host "✓ Rolled back to: $backupFile" -ForegroundColor Yellow
    Write-Host "`nAnalyze errors and adjust strategy.`n" -ForegroundColor Yellow
    
    exit 1
}
