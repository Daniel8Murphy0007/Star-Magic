# Comprehensive batch extraction of ALL 140 missing classes
# This will extract complete class definitions and create PhysicsTerm wrappers

$missingList = Import-Csv "missing_classes_extraction_list.csv"
$totalClasses = $missingList.Count
$successCount = 0
$failCount = 0
$totalMethodsExtracted = 0

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "BATCH CLASS EXTRACTION" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Total classes to extract: $totalClasses`n" -ForegroundColor Yellow

# Prepare output files
$headerContent = @"
// extracted_all_missing_classes.h
// AUTO-GENERATED: Complete extraction of ALL 140 missing classes
// Total methods: ~1,214
// Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")

#ifndef EXTRACTED_ALL_MISSING_CLASSES_H
#define EXTRACTED_ALL_MISSING_CLASSES_H

#include <map>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <memory>
#include <iostream>

using cdouble = std::complex<double>;

"@

$wrapperContent = @"
// physics_term_wrappers_all.cpp
// AUTO-GENERATED: PhysicsTerm wrappers for all 140 extracted classes
// Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")

#include "extracted_all_missing_classes.h"

"@

$registrationContent = @"
// register_all_extracted_classes.cpp
// AUTO-GENERATED: Registration function for all extracted classes
// Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")

void registerAllExtractedClasses(CalculatorCore& core) {
    g_logger.log("Registering all extracted physics classes...", 1);
    
"@

$progressInterval = [Math]::Max(1, [Math]::Floor($totalClasses / 20))

foreach($item in $missingList) {
    $className = $item.Class
    $fileName = $item.File
    $methodCount = [int]$item.MethodCount
    
    if(($successCount + $failCount) % $progressInterval -eq 0) {
        $pct = [Math]::Round((($successCount + $failCount) / $totalClasses) * 100, 1)
        Write-Host "Progress: $pct% ($successCount extracted, $failCount skipped)" -ForegroundColor Gray
    }
    
    if(-not (Test-Path $fileName)) {
        Write-Host "  SKIP: $className - File not found: $fileName" -ForegroundColor DarkGray
        $failCount++
        continue
    }
    
    $content = Get-Content $fileName -Raw
    
    # Try to extract class definition
    $classPattern = "class\s+$className\s*[^{]*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\};"
    $classMatch = [regex]::Match($content, $classPattern, [System.Text.RegularExpressions.RegexOptions]::Singleline)
    
    if(-not $classMatch.Success) {
        # Try alternate pattern without trailing semicolon
        $classPattern = "class\s+$className\s*[^{]*\{((?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*)\}"
        $classMatch = [regex]::Match($content, $classPattern, [System.Text.RegularExpressions.RegexOptions]::Singleline)
    }
    
    if($classMatch.Success) {
        $classDefinition = $classMatch.Value
        $actualMethods = ([regex]::Matches($classDefinition, 'compute\w+\s*\(')).Count
        
        # Add class to header
        $headerContent += "// ========== $className from $fileName ==========`n"
        $headerContent += "// Methods: $actualMethods (expected: $methodCount)`n"
        $headerContent += $classDefinition + "`n`n"
        
        # Create PhysicsTerm wrapper
        $wrapperContent += @"
// Wrapper for $className
class ${className}_Term : public PhysicsTerm {
private:
    std::shared_ptr<$className> instance;
public:
    ${className}_Term() : instance(std::make_shared<$className>()) {
        setMetadata("version", "1.0");
        setMetadata("source", "$fileName");
        setMetadata("class", "$className");
        setMetadata("methods", "$actualMethods");
    }
    
    double compute(double t, const SystemParams& params) const override {
        // TODO: Implement primary computation
        return 0.0;
    }
    
    std::string getDescription() const override {
        return "$className from $fileName";
    }
};

"@
        
        # Add to registration
        $registrationContent += "    core.registerPhysicsTerm(`"$className`", std::make_unique<${className}_Term>(), `"extracted`");`n"
        
        $successCount++
        $totalMethodsExtracted += $actualMethods
    }
    else {
        $failCount++
    }
}

# Finalize files
$headerContent += "`n#endif // EXTRACTED_ALL_MISSING_CLASSES_H`n"
$registrationContent += "`n    g_logger.log(`"Registered $successCount extracted classes with $totalMethodsExtracted methods`", 1);`n}`n"

# Write files
$headerContent | Out-File "extracted_all_missing_classes.h" -Encoding UTF8
$wrapperContent | Out-File "physics_term_wrappers_all.cpp" -Encoding UTF8
$registrationContent | Out-File "register_all_extracted_classes.cpp" -Encoding UTF8

# Summary
Write-Host "`n========================================" -ForegroundColor Green
Write-Host "BATCH EXTRACTION COMPLETE" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host "Classes processed: $totalClasses" -ForegroundColor White
Write-Host "Successfully extracted: $successCount" -ForegroundColor Green
Write-Host "Skipped/Failed: $failCount" -ForegroundColor Yellow
Write-Host "Total methods extracted: $totalMethodsExtracted" -ForegroundColor Cyan
Write-Host "`nOutput files:" -ForegroundColor White
Write-Host "  extracted_all_missing_classes.h ($successCount classes)" -ForegroundColor Gray
Write-Host "  physics_term_wrappers_all.cpp ($successCount wrappers)" -ForegroundColor Gray
Write-Host "  register_all_extracted_classes.cpp (registration function)" -ForegroundColor Gray
Write-Host ""

# Generate summary report
$summaryReport = @"
# Batch Class Extraction Summary
Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")

## Statistics
- Total classes targeted: $totalClasses
- Successfully extracted: $successCount
- Failed/Skipped: $failCount
- Success rate: $([Math]::Round(($successCount / $totalClasses) * 100, 1))%
- Total methods extracted: $totalMethodsExtracted

## Integration Status
- Previously integrated in MAIN_1_CoAnQi.cpp: 291 terms
- Newly extracted: $totalMethodsExtracted terms
- **New total: $($291 + $totalMethodsExtracted) physics terms**

## Next Steps
1. Review extracted_all_missing_classes.h for compilation errors
2. Implement compute() methods in wrappers
3. Add #include "extracted_all_missing_classes.h" to MAIN_1_CoAnQi.cpp
4. Call registerAllExtractedClasses(core) in main()
5. Recompile and test

"@

$summaryReport | Out-File "BATCH_EXTRACTION_SUMMARY.md" -Encoding UTF8
Write-Host "Summary report: BATCH_EXTRACTION_SUMMARY.md`n" -ForegroundColor Cyan
