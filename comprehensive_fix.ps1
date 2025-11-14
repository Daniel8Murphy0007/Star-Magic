# Comprehensive Module Fixes
# 1. Add M_PI support
# 2. Comment all documentation properly

param([string]$ModulesDir = "Core/Modules")

Write-Host "Applying comprehensive fixes to $ModulesDir" -ForegroundColor Cyan

$files = Get-ChildItem -Path $ModulesDir -Filter "*.cpp"
$fixedCount = 0

foreach ($file in $files) {
    $content = Get-Content $file.FullName -Raw
    $modified = $false
    
    # Fix 1: Add M_PI support if needed
    if ($content -match '\bM_PI\b' -and $content -notmatch '#define _USE_MATH_DEFINES') {
        # Insert before first #include
        $content = $content -replace '(#include)', "#define _USE_MATH_DEFINES`n`$1"
        $modified = $true
    }
    
    # Fix 2: Comment out ALL lines between "Module Evaluation" and end of file that aren't already commented
    if ($content -match '(?ms)(Module Evaluation\s*$)') {
        # Find position of "Module Evaluation"
        if ($content -match '(?ms)^(.*)(\n[A-Za-z0-9_]+Module Evaluation\s*$)(.*$)') {
            $before = $matches[1]
            $evalHeader = $matches[2]
            $after = $matches[3]
            
            # Comment the header
            $evalHeader = $evalHeader -replace '([A-Za-z0-9_]+Module Evaluation)', '// $1'
            
            # Comment all non-empty, non-commented lines in the rest
            $afterLines = $after -split "`n"
            $commentedAfter = ($afterLines | ForEach-Object {
                if ($_ -match '^\s*$') {
                    "//"  # Blank line
                } elseif ($_ -match '^\s*//') {
                    $_   # Already commented
                } else {
                    "// $_"  # Add comment
                }
            }) -join "`n"
            
            $content = $before + $evalHeader + "`n" + $commentedAfter
            $modified = $true
        }
    }
    
    if ($modified) {
        Set-Content -Path $file.FullName -Value $content -NoNewline
        Write-Host "  Fixed: $($file.Name)" -ForegroundColor Yellow
        $fixedCount++
    }
}

Write-Host "`nFixed $fixedCount module files" -ForegroundColor Green
