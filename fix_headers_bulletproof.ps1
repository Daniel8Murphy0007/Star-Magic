# Bulletproof M_PI Header Fix
# Replaces mangled headers with validated pattern from source14.cpp

param([string]$ModulesDir = "Core/Modules")

$validatedPattern = @"
#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
"@

Write-Host "Applying bulletproof M_PI headers to $ModulesDir" -ForegroundColor Cyan

$files = Get-ChildItem -Path $ModulesDir -Filter "*.cpp"
$fixedCount = 0

foreach ($file in $files) {
    $content = Get-Content $file.FullName -Raw
    $modified = $false
    
    # Remove ALL existing _USE_MATH_DEFINES and M_PI defines
    $content = $content -replace '(?m)^#define _USE_MATH_DEFINES\s*\n', ''
    $content = $content -replace '(?m)^#ifndef M_PI\s*\n#define M_PI[^\n]*\n#endif\s*\n', ''
    $content = $content -replace '(?m)^#define M_PI[^\n]*\n', ''
    
    # If file uses M_PI or cmath, ensure proper header
    if ($content -match '\bM_PI\b' -or $content -match '#include <cmath>') {
        # Find first #include
        if ($content -match '(?ms)^(.*?)(#include\s+[<"][^>"]+[>"])(.*)$') {
            $before = $matches[1]
            $firstInclude = $matches[2]
            $after = $matches[3]
            
            # Check if it's <cmath>
            if ($firstInclude -match '<cmath>') {
                # Replace with validated pattern
                $content = $before + $validatedPattern + $after
            } else {
                # Insert validated pattern before first include
                $content = $before + $validatedPattern + "`n" + $firstInclude + $after
            }
            $modified = $true
        }
    }
    
    if ($modified) {
        Set-Content -Path $file.FullName -Value $content -NoNewline
        Write-Host "  Fixed: $($file.Name)" -ForegroundColor Green
        $fixedCount++
    }
}

Write-Host "`nFixed $fixedCount modules with bulletproof M_PI headers" -ForegroundColor Cyan
