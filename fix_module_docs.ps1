# Fix Module Documentation
# Properly comments markdown documentation that interferes with compilation

param(
    [string]$ModulesDir = "Core/Modules"
)

Write-Host "Fixing module documentation in $ModulesDir" -ForegroundColor Cyan

$files = Get-ChildItem -Path $ModulesDir -Filter "*.cpp"
$fixedCount = 0

foreach ($file in $files) {
    $content = Get-Content $file.FullName -Raw
    $modified = $false
    
    # Pattern 1: Fix unescaped markdown documentation sections
    # Look for common problematic patterns that appear outside comments
    
    # Check if file has uncommented documentation headers
    if ($content -match '(?m)^(Strengths|Weaknesses|Summary|ASASSN14liUQFFModule Evaluation)\s*:?\s*$') {
        Write-Host "  Fixing: $($file.Name)" -ForegroundColor Yellow
        
        # Find the start of the documentation block (after the last // comment or code)
        # Typically after example usage comments
        if ($content -match '(?s)(.*// Watermark:.*?\n\n)([\s\S]+$)') {
            $beforeDoc = $matches[1]
            $docSection = $matches[2]
            
            # Comment out each line of the documentation
            $commentedDoc = ($docSection -split "`n" | ForEach-Object {
                if ($_ -match '^\s*$') { 
                    "// $_"  # Keep blank lines with //
                } elseif ($_ -match '^//') {
                    $_  # Already commented
                } else {
                    "// $_"  # Add comment marker
                }
            }) -join "`n"
            
            $newContent = $beforeDoc + $commentedDoc
            Set-Content -Path $file.FullName -Value $newContent -NoNewline
            $fixedCount++
            $modified = $true
        }
    }
    
    # Pattern 2: Fix backtick characters in code
    if ($content -match '`') {
        if (-not $modified) {
            Write-Host "  Fixing backticks: $($file.Name)" -ForegroundColor Yellow
        }
        $content = $content -replace '`', ''
        Set-Content -Path $file.FullName -Value $content -NoNewline
        if (-not $modified) { $fixedCount++ }
        $modified = $true
    }
}

Write-Host "`nFixed $fixedCount module files" -ForegroundColor Green
