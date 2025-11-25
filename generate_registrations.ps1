# ============================================================================
# generate_registrations.ps1
# Phase 4: Generate PhysicsTerm registrations for all classes
# ============================================================================

$ErrorActionPreference = "Stop"
$file = "MAIN_1_CoAnQi.cpp"
$output = "generated_registrations.txt"

Write-Host "`n=== PHASE 4: GENERATING PHYSICS TERM REGISTRATIONS ===" -ForegroundColor Cyan
Write-Host "Scanning: $file`n" -ForegroundColor Yellow

# Read file
$content = Get-Content $file -Raw

# Pattern to find PhysicsTerm subclasses
$classPattern = 'class\s+(\w+)\s*:\s*(?:public\s+)?PhysicsTerm'
$matches = [regex]::Matches($content, $classPattern)

Write-Host "[1/3] Found $($matches.Count) PhysicsTerm classes" -ForegroundColor Green

# Generate registrations
$registrations = @()
$classNames = @{}

foreach ($match in $matches) {
    $className = $match.Groups[1].Value
    
    # Skip base class
    if ($className -eq "PhysicsTerm") { continue }
    
    # Track duplicates
    if ($classNames.ContainsKey($className)) {
        $classNames[$className]++
        # Add suffix for duplicates
        $uniqueName = "${className}_v$($classNames[$className])"
        Write-Host "   Duplicate: $className → registered as $uniqueName" -ForegroundColor Gray
    } else {
        $classNames[$className] = 1
        $uniqueName = $className
    }
    
    $registrations += "    core.registerPhysicsTerm(`"$uniqueName`", std::make_unique<$className>(), `"auto-generated`");"
}

Write-Host "`n[2/3] Generated $($registrations.Count) registration calls" -ForegroundColor Green

# Write to output file
$registrations | Out-File $output -Encoding UTF8

Write-Host "`n[3/3] Output written to: $output" -ForegroundColor Green

# Show summary
Write-Host "`n=== REGISTRATION SUMMARY ===" -ForegroundColor White
Write-Host "Total PhysicsTerms found: $($matches.Count)" -ForegroundColor White
Write-Host "Unique classes: $($classNames.Count)" -ForegroundColor White
Write-Host "Registrations generated: $($registrations.Count)" -ForegroundColor White

# Show sample
Write-Host "`nFirst 10 registrations:" -ForegroundColor Yellow
$registrations | Select-Object -First 10 | ForEach-Object { Write-Host "  $_" -ForegroundColor Gray }

Write-Host "`nLast 10 registrations:" -ForegroundColor Yellow
$registrations | Select-Object -Last 10 | ForEach-Object { Write-Host "  $_" -ForegroundColor Gray }

Write-Host "`n=== PHASE 4 COMPLETE ===" -ForegroundColor Green
Write-Host "✅ Registrations file: $output" -ForegroundColor Green
Write-Host "`nNext: Manually insert registrations into registerAllPhysicsTerms() function`n"
