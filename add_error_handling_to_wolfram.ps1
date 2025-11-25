# Add Try-Catch Error Handling to Wolfram Registration Function
# This will wrap each registration in try-catch to identify the 14 failing classes

Write-Host "=== ADDING ERROR HANDLING TO WOLFRAM REGISTRATIONS ===" -ForegroundColor Cyan
Write-Host ""

$wolframFile = ".\wolfram_physics_classes.cpp"
$backupFile = ".\wolfram_physics_classes_backup_before_error_handling.cpp"

# Backup original file
Write-Host "[1/4] Creating backup..." -ForegroundColor Yellow
Copy-Item $wolframFile $backupFile -Force
Write-Host "      Backup created: $backupFile" -ForegroundColor Green
Write-Host ""

# Read the file
Write-Host "[2/4] Reading wolfram_physics_classes.cpp (132,533 lines)..." -ForegroundColor Yellow
$content = Get-Content $wolframFile -Raw
Write-Host "      File loaded successfully" -ForegroundColor Green
Write-Host ""

# Find the registration function and wrap each registration in try-catch
Write-Host "[3/4] Wrapping each registration in try-catch block..." -ForegroundColor Yellow

# Replace the function header to add error tracking
$newFunctionStart = @'
void registerAllWolframPhysicsTerms(CalculatorCore& core) {
    int successCount = 0;
    int failureCount = 0;
    std::vector<std::string> failedRegistrations;
    
    // Lambda to safely register with error handling
    auto safeRegister = [&](const std::string& name, auto termFactory, const std::string& category) {
        try {
            core.registerPhysicsTerm(name, termFactory(), category);
            successCount++;
        } catch (const std::exception& e) {
            failureCount++;
            failedRegistrations.push_back(name + " [" + std::string(e.what()) + "]");
            g_logger.log("FAILED to register: " + name + " - " + std::string(e.what()), 1);
        } catch (...) {
            failureCount++;
            failedRegistrations.push_back(name + " [unknown error]");
            g_logger.log("FAILED to register: " + name + " - unknown error", 1);
        }
    };
    
    g_logger.log("=== WOLFRAM REGISTRATION WITH ERROR HANDLING ===", 1);
    
'@

# Replace all core.registerPhysicsTerm calls with safeRegister calls
$content = $content -replace 'void registerAllWolframPhysicsTerms\(CalculatorCore& core\) \{', $newFunctionStart
$content = $content -replace 'core\.registerPhysicsTerm\(', 'safeRegister('

# Add summary logging at the end of the function (before closing brace)
$summaryLog = @'

    // Summary report
    g_logger.log("=== WOLFRAM REGISTRATION SUMMARY ===", 1);
    g_logger.log("Successful registrations: " + std::to_string(successCount), 1);
    g_logger.log("Failed registrations: " + std::to_string(failureCount), 1);
    if (failureCount > 0) {
        g_logger.log("=== FAILED REGISTRATIONS LIST ===", 1);
        for (const auto& failed : failedRegistrations) {
            g_logger.log("  - " + failed, 1);
        }
    }
    g_logger.log("=== END WOLFRAM REGISTRATION ===", 1);
}
'@

# Find the last line of the function (closing brace) and add summary before it
$lastBracePattern = '(?s)(.*core\.registerPhysicsTerm.*\);)\s*\}'
if ($content -match $lastBracePattern) {
    $content = $content -replace '\}\s*$', $summaryLog
}

Write-Host "      Error handling added successfully" -ForegroundColor Green
Write-Host ""

# Write the modified content
Write-Host "[4/4] Writing modified file..." -ForegroundColor Yellow
Set-Content $wolframFile -Value $content -NoNewline
Write-Host "      File updated successfully" -ForegroundColor Green
Write-Host ""

Write-Host "=== COMPLETE ===" -ForegroundColor Cyan
Write-Host "Modified file: $wolframFile" -ForegroundColor White
Write-Host "Backup file: $backupFile" -ForegroundColor White
Write-Host ""
Write-Host "Next step: Rebuild the project to see which 14 registrations fail" -ForegroundColor Yellow
Write-Host "Command: cmake --build build_msvc --config Release --target MAIN_1_CoAnQi" -ForegroundColor Gray
