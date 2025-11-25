# ==============================================================================
# Generate Missing MAIN Registrations Script
# Purpose: Extract 81 unregistered PhysicsTerm classes and generate registration calls
# Compliance: C++20 / MSVC 14.44+ / PowerShell 5.1+
# ==============================================================================

Write-Host "`n=== MISSING REGISTRATION GENERATOR ===" -ForegroundColor Cyan
Write-Host "Extracting unregistered PhysicsTerm classes from MAIN_1_CoAnQi.cpp`n" -ForegroundColor Yellow

# Step 1: Extract all PhysicsTerm class definitions
Write-Host "[1/4] Extracting all PhysicsTerm class definitions..." -ForegroundColor Green
$allClasses = @()
$classMatches = Select-String -Path "MAIN_1_CoAnQi.cpp" -Pattern "^class\s+(\w+)\s*:\s*public\s+PhysicsTerm" -AllMatches

foreach ($match in $classMatches) {
    $className = $match.Matches.Groups[1].Value
    # Remove "Term" suffix to get registration name
    $regName = $className -replace 'Term$', ''
    $allClasses += @{
        ClassName = $className
        RegName = $regName
        LineNumber = $match.LineNumber
    }
}

Write-Host "  Found: $($allClasses.Count) PhysicsTerm classes" -ForegroundColor White

# Step 2: Extract all registered term names
Write-Host "[2/4] Extracting currently registered term names..." -ForegroundColor Green
$registered = @()
$regMatches = Select-String -Path "MAIN_1_CoAnQi.cpp" -Pattern 'core\.registerPhysicsTerm\("([^"]+)"' -AllMatches

foreach ($match in $regMatches) {
    $registered += $match.Matches.Groups[1].Value
}

Write-Host "  Found: $($registered.Count) registered terms" -ForegroundColor White

# Step 3: Find unregistered classes
Write-Host "[3/4] Identifying unregistered classes..." -ForegroundColor Green
$unregistered = @()

foreach ($class in $allClasses) {
    # Check if registration name exists in registered list
    if ($registered -notcontains $class.RegName -and $registered -notcontains $class.ClassName) {
        $unregistered += $class
    }
}

Write-Host "  Found: $($unregistered.Count) unregistered classes`n" -ForegroundColor Red

# Step 4: Generate registration calls
Write-Host "[4/4] Generating registration calls..." -ForegroundColor Green

if ($unregistered.Count -eq 0) {
    Write-Host "`n✅ SUCCESS: All classes are already registered!" -ForegroundColor Green
    exit 0
}

# Generate registration code
$registrationCode = @()
$registrationCode += "    // ========== BATCH 17: MISSING MAIN REGISTRATIONS ($($unregistered.Count) classes) =========="
$registrationCode += "    // Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
$registrationCode += "    // Auto-registered previously undefined PhysicsTerm classes"
$registrationCode += ""

foreach ($class in $unregistered | Sort-Object { $_.ClassName }) {
    $regLine = "    core.registerPhysicsTerm(`"$($class.RegName)`", std::make_unique<$($class.ClassName)>(), `"auto-registered`");"
    $registrationCode += $regLine
}

$registrationCode += ""
$registrationCode += "    // ========== END BATCH 17: $($unregistered.Count) REGISTRATIONS =========="

# Save to file
$outputFile = "missing_main_registrations_$(Get-Date -Format 'yyyyMMdd_HHmmss').txt"
$registrationCode | Out-File $outputFile -Encoding UTF8

Write-Host "`n✅ Generated $($unregistered.Count) registration calls" -ForegroundColor Green
Write-Host "📁 Saved to: $outputFile`n" -ForegroundColor Cyan

# Display summary
Write-Host "=== SUMMARY ===" -ForegroundColor Cyan
Write-Host "Total PhysicsTerm classes:     $($allClasses.Count)" -ForegroundColor White
Write-Host "Currently registered:          $($registered.Count)" -ForegroundColor Green
Write-Host "Unregistered (generated):      $($unregistered.Count)" -ForegroundColor Red
Write-Host "Target total registrations:    $($registered.Count + $unregistered.Count)`n" -ForegroundColor Yellow

# Display first 10 unregistered classes
Write-Host "=== FIRST 10 UNREGISTERED CLASSES ===" -ForegroundColor Cyan
$unregistered | Select-Object -First 10 | ForEach-Object {
    Write-Host "  - $($_.ClassName) (line $($_.LineNumber))" -ForegroundColor Magenta
}

if ($unregistered.Count -gt 10) {
    Write-Host "  ... and $($unregistered.Count - 10) more`n" -ForegroundColor Gray
}

# Display next steps
Write-Host "=== NEXT STEPS ===" -ForegroundColor Yellow
Write-Host "1. Review generated registration calls in: $outputFile"
Write-Host "2. Open MAIN_1_CoAnQi.cpp and locate line 21637 (before registerAllWolframPhysicsTerms)"
Write-Host "3. Insert the generated registration calls"
Write-Host "4. Rebuild: cmake --build build_msvc --config Release"
Write-Host "5. Verify: .\build_msvc\Release\MAIN_1_CoAnQi.exe`n"

# Save detailed report
$reportFile = "missing_registrations_report_$(Get-Date -Format 'yyyyMMdd_HHmmss').csv"
$unregistered | Select-Object ClassName, RegName, LineNumber | 
    Export-Csv $reportFile -NoTypeInformation

Write-Host "📊 Detailed report saved to: $reportFile`n" -ForegroundColor Cyan

Write-Host "✅ COMPLETE - Missing registration generation finished" -ForegroundColor Green
