# PHASE 2: Uncomment bulk patterns (LINE-BY-LINE approach)
# Removes comment prefixes from EVERY line in the commented range

$file = "MAIN_1_CoAnQi.cpp"
$startLine = 31500
$endLine = 102427
$backupFile = "MAIN_1_CoAnQi_backup_before_uncomment_v2.cpp"

Write-Host "`n=== PHASE 2: UNCOMMENTING BULK PATTERNS (PER-LINE) ===" -ForegroundColor Cyan
Write-Host "File: $file"
Write-Host "Processing lines $startLine-$endLine..."

# Create backup
Write-Host "`n[1/4] Creating backup: $backupFile" -ForegroundColor Green
Copy-Item $file $backupFile -Force
Write-Host "      Backup created successfully."

# Read file
Write-Host "`n[2/4] Reading file..." -ForegroundColor Green
$lines = Get-Content $file
$originalSize = (Get-Item $file).Length / 1MB
Write-Host "      File size: $($originalSize.ToString('F1')) MB"
Write-Host "      Total lines: $($lines.Count)"

# Process lines
Write-Host "`n[3/4] Removing comment prefixes from lines $startLine-$endLine..." -ForegroundColor Green
$removedCount = 0
$processedLines = 0

for ($i = 0; $i -lt $lines.Count; $i++) {
    $lineNum = $i + 1  # 1-indexed
    
    # Only process lines in the target range
    if ($lineNum -ge $startLine -and $lineNum -le $endLine) {
        $originalLine = $lines[$i]
        
        # Remove ALL comment prefixes with tags
        $newLine = $originalLine -replace '^//\s*\[Duplicate class definition\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Duplicate\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Bulk comment\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Incomplete class extraction\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Incomplete\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[GUI\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Parser\]\s*', ''
        $newLine = $newLine -replace '^//\s*\[Symbolic\]\s*', ''
        
        if ($newLine -ne $originalLine) {
            $removedCount++
        }
        
        $lines[$i] = $newLine
        $processedLines++
        
        # Progress update every 10,000 lines
        if ($processedLines % 10000 -eq 0) {
            Write-Host "      Progress: $processedLines / $($endLine - $startLine + 1) lines processed, $removedCount prefixes removed" -ForegroundColor DarkGray
        }
    }
}

Write-Host "      Total lines processed: $processedLines"
Write-Host "      Total comment prefixes removed: $removedCount"

# Write file
Write-Host "`n[4/4] Writing updated file..." -ForegroundColor Green
$lines | Set-Content $file -Encoding UTF8
$newSize = (Get-Item $file).Length / 1MB
$sizeDiff = ($newSize - $originalSize) * 1024  # Convert to KB
Write-Host "      New file size: $($newSize.ToString('F1')) MB"
Write-Host "      Size change: $($sizeDiff.ToString('F1')) KB"

Write-Host "`n=== PHASE 2 COMPLETE ===" -ForegroundColor Green
Write-Host "✅ $removedCount comment prefixes removed from $processedLines lines" -ForegroundColor Green
Write-Host "✅ File: $file updated" -ForegroundColor Green
Write-Host "✅ Backup: $backupFile created" -ForegroundColor Green
Write-Host "`nNext: Run PHASE 3 (build and resolve compilation errors)" -ForegroundColor Yellow
