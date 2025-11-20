# Extract complete class definitions from source files and create PhysicsTerm wrappers
# Target: Extract all 140 missing classes with 1,214 methods

$targetClasses = @(
    @{Class="UQFFBuoyancyModule"; File="Source155.cpp"; Methods=102},
    @{Class="UQFFModule5"; File="Source5.cpp"; Methods=63},
    @{Class="NGC346UQFFModule"; File="source81.cpp"; Methods=36},
    @{Class="SurfaceMagneticFieldModule"; File="source121.cpp"; Methods=34},
    @{Class="UQFFBuoyancyCNBModule"; File="Source156.cpp"; Methods=30},
    @{Class="Abell2256UQFFModule"; File="source134.cpp"; Methods=30},
    @{Class="LagoonNebulaUQFFModule"; File="Source143.cpp"; Methods=30},
    @{Class="CentaurusAUQFFModule"; File="source133.cpp"; Methods=29},
    @{Class="MUGEModule"; File="source86.cpp"; Methods=24},
    @{Class="FluidSolver"; File="source4.cpp"; Methods=24}
)

$outputHeader = @"
// extracted_complete_classes.h
// AUTO-GENERATED: Complete class extractions from source files
// These classes will be wrapped as PhysicsTerm subclasses

#ifndef EXTRACTED_COMPLETE_CLASSES_H
#define EXTRACTED_COMPLETE_CLASSES_H

#include <map>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <memory>

// Forward declarations
class PhysicsTerm;

"@

$outputImplementation = @"
// extracted_complete_implementations.cpp
// AUTO-GENERATED: PhysicsTerm wrappers for extracted classes

"@

$totalExtracted = 0
$totalMethods = 0

foreach($target in $targetClasses) {
    $className = $target.Class
    $fileName = $target.File
    $expectedMethods = $target.Methods
    
    Write-Host "Extracting $className from $fileName..." -ForegroundColor Cyan
    
    if(-not (Test-Path $fileName)) {
        Write-Host "  WARNING: File not found: $fileName" -ForegroundColor Yellow
        continue
    }
    
    $content = Get-Content $fileName -Raw
    
    # Extract class definition
    $classPattern = "class\s+$className\s*\{(.*?)(?=\n\s*class\s|\n#endif|\z)"
    $classMatch = [regex]::Match($content, $classPattern, [System.Text.RegularExpressions.RegexOptions]::Singleline)
    
    if($classMatch.Success) {
        $classBody = $classMatch.Value
        $methodCount = ([regex]::Matches($classBody, 'compute\w+\s*\(')).Count
        
        Write-Host "  Found class definition ($methodCount compute methods)" -ForegroundColor Green
        
        # Add to header
        $outputHeader += "// ========== $className from $fileName ==========`n"
        $outputHeader += "// Expected methods: $expectedMethods, Found: $methodCount`n"
        $outputHeader += $classBody + "`n`n"
        
        # Create PhysicsTerm wrapper
        $wrapperClass = @"
// PhysicsTerm wrapper for $className
class ${className}_PhysicsTerm : public PhysicsTerm {
private:
    std::unique_ptr<$className> instance;
    std::string systemName;
public:
    ${className}_PhysicsTerm(const std::string& system = "default") 
        : instance(std::make_unique<$className>()), systemName(system) {
        setMetadata("version", "1.0");
        setMetadata("source", "$fileName");
        setMetadata("wrapped_class", "$className");
        setMetadata("method_count", "$methodCount");
    }
    
    double compute(double t, const SystemParams& params) const override {
        // Call primary compute method from wrapped class
        // TODO: Map params to class-specific parameters
        return 0.0; // Placeholder
    }
    
    std::string getDescription() const override {
        return "$className from $fileName - $methodCount compute methods";
    }
};

"@
        $outputImplementation += $wrapperClass
        
        $totalExtracted++
        $totalMethods += $methodCount
    }
    else {
        Write-Host "  ERROR: Could not extract class definition" -ForegroundColor Red
    }
}

$outputHeader += "`n#endif // EXTRACTED_COMPLETE_CLASSES_H`n"

# Write output files
$outputHeader | Out-File "extracted_complete_classes.h" -Encoding UTF8
$outputImplementation | Out-File "extracted_complete_implementations.cpp" -Encoding UTF8

Write-Host "`n========================================" -ForegroundColor Green
Write-Host "EXTRACTION COMPLETE" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Green
Write-Host "Classes extracted: $totalExtracted" -ForegroundColor Cyan
Write-Host "Total methods: $totalMethods" -ForegroundColor Cyan
Write-Host "`nOutput files:" -ForegroundColor White
Write-Host "  extracted_complete_classes.h" -ForegroundColor Gray
Write-Host "  extracted_complete_implementations.cpp" -ForegroundColor Gray
Write-Host ""
