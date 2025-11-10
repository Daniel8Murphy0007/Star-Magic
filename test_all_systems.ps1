# Comprehensive MAIN_1.cpp System Test Suite
# Tests all 26 predefined astrophysical systems in the UQFF calculator
# Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')

Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "MAIN_1.cpp Comprehensive System Test" -ForegroundColor Cyan
Write-Host "Testing all 26 UQFF Calculator Systems" -ForegroundColor Cyan
Write-Host "========================================`n" -ForegroundColor Cyan

# Define all 26 systems from MAIN_1.cpp line 1198+
$systems = @(
    "ESO 137-001",
    "Black Hole Pairs",
    "SN 1006",
    "Eta Carinae",
    "Galactic Center",
    "Kepler's Supernova Remnant",
    "NGC 1365",
    "Vela Pulsar",
    "ASASSN-14li",
    "El Gordo",
    "Magnetar SGR 1745-2900",
    "Tapestry of Blazing Starbirth NGC 2264",
    "Westerlund 2",
    "Pillars of Creation M16",
    "Rings of Relativity",
    "Chandra Archive Collection",
    "Cassiopeia",
    "3C273",
    "Cen A AGN",
    "UHZ1 AGN",
    "Geminga",
    "GW170817",
    "NGC 1068",
    "PJ352-15",
    "Quasar Survey (Typical)",
    "GSN 069"
)

# Results storage
$results = @()
$successCount = 0
$failCount = 0

# Test each system
foreach ($system in $systems) {
    Write-Host "Testing: $system..." -NoNewline
    
    # Run MAIN_1.exe with system name, no parameter override
    $input = "$system`nn`n"
    $output = ($input | .\build\Debug\MAIN_1.exe 2>&1) | Out-String
    
    # Parse output for key values
    $F_U_Bi_i = $null
    $g_rt = $null
    $L_X = $null
    $B0 = $null
    $omega0 = $null
    
    # Extract F_U_Bi_i (looking for "Original UQFF F_U_Bi_i: <value> N")
    if ($output -match "Original UQFF F_U_Bi_i:\s*([\d.eE+-]+|inf)") {
        if ($matches[1] -eq "inf") {
            $F_U_Bi_i = [double]::PositiveInfinity
        }
        else {
            $F_U_Bi_i = [double]$matches[1]
        }
    }
    
    # Extract g(r,t) (looking for "Compressed g(r,t): <value> m/s^2")
    if ($output -match "Compressed g\(r,t\):\s*([\d.eE+-]+|inf)") {
        if ($matches[1] -eq "inf") {
            $g_rt = [double]::PositiveInfinity
        }
        else {
            $g_rt = [double]$matches[1]
        }
    }
    
    # Extract L_X (looking for "L_X: <value> W")
    if ($output -match "L_X:\s*([\d.eE+-]+)") {
        $L_X = [double]$matches[1]
    }
    
    # Extract B0 (looking for "B0: <value> T")
    if ($output -match "B0:\s*([\d.eE+-]+)") {
        $B0 = [double]$matches[1]
    }
    
    # Extract omega0 (looking for "omega0: <value> s^-1")
    if ($output -match "omega0:\s*([\d.eE+-]+)") {
        $omega0 = [double]$matches[1]
    }
    
    # Determine system category based on characteristics
    $category = "Unknown"
    if ($system -match "Pulsar|Magnetar|Geminga") {
        $category = "Pulsar/Magnetar"
    }
    elseif ($system -match "Galaxy|NGC|Cen A|ESO") {
        $category = "Galaxy"
    }
    elseif ($system -match "Black Hole|AGN|Quasar|3C273|PJ352|GSN") {
        $category = "Black Hole/AGN"
    }
    elseif ($system -match "Cluster|El Gordo|Westerlund|Tapestry") {
        $category = "Cluster"
    }
    elseif ($system -match "SN |Supernova|Cassiopeia|Kepler") {
        $category = "Supernova Remnant"
    }
    elseif ($system -match "GW|Merger") {
        $category = "Merger Event"
    }
    elseif ($system -match "Eta Carinae") {
        $category = "Massive Star"
    }
    elseif ($system -match "Pillars|Nebula") {
        $category = "Nebula"
    }
    elseif ($system -match "ASASSN|TDE") {
        $category = "Transient"
    }
    elseif ($system -match "Rings|Einstein") {
        $category = "Gravitational Lens"
    }
    elseif ($system -match "Chandra") {
        $category = "Archive/Generic"
    }
    
    # Check if test passed
    if ($F_U_Bi_i -ne $null -and $g_rt -ne $null) {
        Write-Host " [OK]" -ForegroundColor Green
        $successCount++
        $status = "PASS"
    }
    else {
        Write-Host " [FAIL]" -ForegroundColor Red
        $failCount++
        $status = "FAIL"
    }
    
    # Store result
    $results += [PSCustomObject]@{
        System   = $system
        Category = $category
        F_U_Bi_i = $F_U_Bi_i
        g_rt     = $g_rt
        L_X      = $L_X
        B0       = $B0
        omega0   = $omega0
        Status   = $status
    }
}

# Display summary table
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "RESULTS SUMMARY" -ForegroundColor Cyan
Write-Host "========================================`n" -ForegroundColor Cyan

# Display results by category
$categories = $results | Group-Object -Property Category | Sort-Object Name

foreach ($cat in $categories) {
    Write-Host "`n--- $($cat.Name) ($($cat.Count) systems) ---" -ForegroundColor Yellow
    
    foreach ($result in $cat.Group) {
        $statusColor = if ($result.Status -eq "PASS") { "Green" } else { "Red" }
        Write-Host "  $($result.System): " -NoNewline
        Write-Host "[$($result.Status)]" -ForegroundColor $statusColor
        
        if ($result.F_U_Bi_i -ne $null) {
            Write-Host "    F_U_Bi_i = $($result.F_U_Bi_i.ToString('E2')) N" -ForegroundColor Gray
        }
        if ($result.g_rt -ne $null) {
            Write-Host "    g(r,t) = $($result.g_rt.ToString('E2')) m/s^2" -ForegroundColor Gray
        }
        if ($result.L_X -ne $null) {
            Write-Host "    L_X = $($result.L_X.ToString('E2')) W" -ForegroundColor Gray
        }
    }
}

# Compute dynamic ranges
$validResults = $results | Where-Object { $_.F_U_Bi_i -ne $null -and $_.g_rt -ne $null }

if ($validResults.Count -gt 0) {
    # Filter out inf values for range calculation
    $finiteF = $validResults | Where-Object { [Math]::IsFinite($_.F_U_Bi_i) }
    $finiteG = $validResults | Where-Object { [Math]::IsFinite($_.g_rt) }
    
    if ($finiteF.Count -gt 0) {
        $minF = ($finiteF | Measure-Object -Property F_U_Bi_i -Minimum).Minimum
        $maxF = ($finiteF | Measure-Object -Property F_U_Bi_i -Maximum).Maximum
        $F_range = [Math]::Log10($maxF / $minF)
        
        Write-Host "`n========================================" -ForegroundColor Cyan
        Write-Host "DYNAMIC RANGE ANALYSIS" -ForegroundColor Cyan
        Write-Host "========================================" -ForegroundColor Cyan
        Write-Host "F_U_Bi_i Range: $($minF.ToString('E2')) to $($maxF.ToString('E2')) N" -ForegroundColor White
        Write-Host "  Dynamic Range: $([Math]::Round($F_range, 2)) orders of magnitude" -ForegroundColor White
    }
    
    if ($finiteG.Count -gt 0) {
        $minG = ($finiteG | Measure-Object -Property g_rt -Minimum).Minimum
        $maxG = ($finiteG | Measure-Object -Property g_rt -Maximum).Maximum
        $G_range = [Math]::Log10($maxG / $minG)
        
        Write-Host "g(r,t) Range: $($minG.ToString('E2')) to $($maxG.ToString('E2')) m/s^2" -ForegroundColor White
        Write-Host "  Dynamic Range: $([Math]::Round($G_range, 2)) orders of magnitude" -ForegroundColor White
    }
}

# Category breakdown
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "CATEGORY BREAKDOWN" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan

$categoryStats = $results | Group-Object -Property Category | Sort-Object Count -Descending
foreach ($stat in $categoryStats) {
    Write-Host "$($stat.Name): $($stat.Count) systems" -ForegroundColor White
}

# Final status
Write-Host "`n========================================" -ForegroundColor Cyan
Write-Host "FINAL STATUS" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Total Systems Tested: $($systems.Count)" -ForegroundColor White
Write-Host "Passed: $successCount" -ForegroundColor Green
Write-Host "Failed: $failCount" -ForegroundColor Red
Write-Host "Success Rate: $([Math]::Round(($successCount / $systems.Count) * 100, 2))%" -ForegroundColor $(if ($successCount -eq $systems.Count) { "Green" } else { "Yellow" })

# Export to CSV for analysis
$csvPath = "test_results_$(Get-Date -Format 'yyyyMMdd_HHmmss').csv"
$results | Export-Csv -Path $csvPath -NoTypeInformation
Write-Host "`nResults exported to: $csvPath" -ForegroundColor Cyan

Write-Host "`n========================================`n" -ForegroundColor Cyan
