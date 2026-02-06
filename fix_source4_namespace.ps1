# fix_source4_namespace.ps1
# Moves the SOURCE4 namespace from after main() to before main()
# This fixes the build error where SOURCE4 is used inside main() but defined after

$file = "MAIN_1_CoAnQi.cpp"
$backup = "MAIN_1_CoAnQi.cpp.backup_source4fix"

Write-Host "=== SOURCE4 Namespace Fix Script ===" -ForegroundColor Cyan
Write-Host "Reading file: $file"

# Read the file content
$content = Get-Content $file -Raw

# Create backup
Write-Host "Creating backup: $backup"
$content | Set-Content $backup -NoNewline

# Find the SOURCE4 namespace boundaries
# Start: "// ============================================================================"
#        "// SOURCE4: UNIFIED FIELD THEORY"
# End: "} // namespace SOURCE4"

$startMarker = "// ============================================================================`r`n// SOURCE4: UNIFIED FIELD THEORY"
$endMarker = "} // namespace SOURCE4"

# Find start position (look for the comment block before namespace SOURCE4)
$startPattern = '// ============================================================================\r?\n// SOURCE4: UNIFIED FIELD THEORY'
$startMatch = [regex]::Match($content, $startPattern)

if (-not $startMatch.Success) {
    Write-Host "ERROR: Could not find SOURCE4 start marker" -ForegroundColor Red
    exit 1
}

$source4Start = $startMatch.Index
Write-Host "Found SOURCE4 start at position: $source4Start"

# Find end position
$endIndex = $content.IndexOf("} // namespace SOURCE4")
if ($endIndex -lt 0) {
    Write-Host "ERROR: Could not find SOURCE4 end marker" -ForegroundColor Red
    exit 1
}

$source4End = $endIndex + "} // namespace SOURCE4".Length
Write-Host "Found SOURCE4 end at position: $source4End"

# Extract SOURCE4 namespace block
$source4Block = $content.Substring($source4Start, $source4End - $source4Start)
Write-Host "Extracted SOURCE4 block: $($source4Block.Length) characters"

# Find where to insert (before "int main(")
$mainPattern = "// ===========================================================================================`r?`n// MAIN FUNCTION - CoAnQi Interactive Calculator"
$mainMatch = [regex]::Match($content, $mainPattern)

if (-not $mainMatch.Success) {
    # Try alternative pattern
    $mainPattern = "int main\(int argc, char \*argv\[\]\)"
    $mainMatch = [regex]::Match($content, $mainPattern)
}

if (-not $mainMatch.Success) {
    Write-Host "ERROR: Could not find main() function" -ForegroundColor Red
    exit 1
}

$insertPosition = $mainMatch.Index
Write-Host "Will insert before position: $insertPosition"

# Build new content:
# 1. Content before insert position
# 2. SOURCE4 block + newlines
# 3. Content from insert position to SOURCE4 start (excluding SOURCE4)
# 4. Content after SOURCE4 end

$beforeInsert = $content.Substring(0, $insertPosition)
$afterInsertBeforeSource4 = $content.Substring($insertPosition, $source4Start - $insertPosition)
$afterSource4 = $content.Substring($source4End)

# Combine
$newContent = $beforeInsert + $source4Block + "`r`n`r`n" + $afterInsertBeforeSource4 + $afterSource4

Write-Host "New content length: $($newContent.Length)"

# Write the new content to a new file first
$fixedFile = "MAIN_1_CoAnQi_fixed.cpp"
Write-Host "Writing fixed file to: $fixedFile"
[System.IO.File]::WriteAllText($fixedFile, $newContent)

Write-Host "Now replacing original file..."
# Close any handles and replace
Remove-Item $file -Force -ErrorAction SilentlyContinue
Rename-Item $fixedFile $file -Force

Write-Host "=== Fix Complete ===" -ForegroundColor Green
Write-Host "Backup saved to: $backup"
Write-Host "Run 'cmake --build build_msvc --config Release --target MAIN_1_CoAnQi' to verify"
