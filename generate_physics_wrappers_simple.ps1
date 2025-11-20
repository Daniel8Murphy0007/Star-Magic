# Simple script to scan source files and generate basic PhysicsTerm placeholders
$ErrorActionPreference = "Stop"

Write-Host "Scanning source files for physics classes..." -ForegroundColor Cyan

$classes = @()

foreach ($file in Get-ChildItem -Filter "source*.cpp" | Sort-Object {[int]($_.Name -replace '\D',"")}) {
    Write-Host "  $($file.Name)" -ForegroundColor Gray
    $content = Get-Content $file.FullName -Raw
    
    # Find class declarations
    $matches = [regex]::Matches($content, '(?ms)^\s*class\s+(\w+(?:Module|Term|Calculator|Field))\s*(?::|{)')
    
    foreach ($m in $matches) {
        $className = $m.Groups[1].Value
        
        # Skip base classes
        if ($className -match '^(PhysicsTerm|UQFFModule|ModuleInterface)$') { continue }
        
        $classes += [PSCustomObject]@{
            SourceFile = $file.Name
            ClassName = $className
        }
    }
}

Write-Host "Found $($classes.Count) classes" -ForegroundColor Green

# Generate simple wrapper code
$output = "/* AUTO-GENERATED WRAPPERS */`n`n"

foreach ($class in $classes) {
    $wrapperName = $class.ClassName -replace 'Module$', 'Term'
    $output += "// $wrapperName from $($class.SourceFile)`n"
}

$output | Set-Content ".\generated_wrappers.txt"
Write-Host "Generated: generated_wrappers.txt" -ForegroundColor Green
