# Fix Type Conversion Issues in Modules
# Converts integer literals to cdouble for complex arithmetic

param(
    [string]$ModulesDir = "Core/Modules"
)

Write-Host "Fixing type conversions in $ModulesDir" -ForegroundColor Cyan

$files = Get-ChildItem -Path $ModulesDir -Filter "*.cpp"
$fixedCount = 0

foreach ($file in $files) {
    $content = Get-Content $file.FullName -Raw
    $originalContent = $content
    
    # Fix pattern: 2 * variables["..."] -> cdouble(2.0) * variables["..."]
    $content = $content -replace '(\s+cdouble\s+\w+\s*=\s*)(\d+)\s*\*\s*(variables\[)', '$1cdouble($2.0) * $3'
    
    # Fix pattern: 4*a -> cdouble(4.0)*a (in quadratic formulas)
    $content = $content -replace '(\s+)(\d+)\s*\*\s*a([^\w])', '$1cdouble($2.0)*a$3'
    
    # Fix pattern: 2*a (specifically in division contexts)
    $content = $content -replace '\(\s*2\s*\*\s*a\s*\)', '(cdouble(2.0)*a)'
    
    # Fix pattern: b*b -> b*b is OK, but 4*a*c needs fixing
    $content = $content -replace '(\s+)4\s*\*\s*a\s*\*\s*c([^\w])', '$1cdouble(4.0)*a*c$2'
    
    if ($content -ne $originalContent) {
        Set-Content -Path $file.FullName -Value $content -NoNewline
        Write-Host "  Fixed: $($file.Name)" -ForegroundColor Yellow
        $fixedCount++
    }
}

Write-Host "`nFixed $fixedCount module files" -ForegroundColor Green
