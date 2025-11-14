# Fix all remaining `double M_PI =` variable assignments that conflict with M_PI macro
# Replace with `pi_value` and update all uses in the same function

$modulesPath = "Core\Modules"
$files = Get-ChildItem -Path $modulesPath -Filter "*.cpp"

Write-Host "Fixing M_PI variable assignments in all modules..."

foreach ($file in $files) {
    $filePath = $file.FullName
    $content = Get-Content -Path $filePath -Raw
    $modified = $false
    
    # Replace `double M_PI = variables["pi"].real();` with pi_value
    if ($content -match 'double\s+M_PI\s*=\s*variables\["pi"\]\.real\(\);') {
        $content = $content -replace 'double\s+M_PI\s*=\s*variables\["pi"\]\.real\(\);', 
                                    'double pi_value = variables["pi"].real();  // M_PI is a macro'
        $modified = $true
        Write-Host "  Fixed M_PI assignment in $($file.Name)"
    }
    
    # Save if modified
    if ($modified) {
        Set-Content -Path $filePath -Value $content -NoNewline
        Write-Host "    Saved $($file.Name)"
    }
}

Write-Host "`nDone! Now manually check functions to replace M_PI with pi_value where needed."
