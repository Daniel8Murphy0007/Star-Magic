# ===================================================================
# Phase 4: Generate C++ PhysicsTerm classes from validated entities
# ===================================================================
# Input: wolfram_extraction/validated_entities.json
# Output: sourceXXX_wolfram.cpp files with proper C++ class definitions
# ===================================================================

param(
    [string]$InputJson = "wolfram_extraction\validated_entities.json",
    [string]$OutputDir = "wolfram_extraction\generated_classes",
    [switch]$DryRun = $false,
    [switch]$SingleFile = $false  # If true, generate single consolidated file
)

# Create output directory
if (-not (Test-Path $OutputDir)) {
    New-Item -ItemType Directory -Path $OutputDir | Out-Null
}

# Load validated entities
if (-not (Test-Path $InputJson)) {
    Write-Error "Validated entities file not found: $InputJson"
    Write-Host "Run Phase 3 (wolfram_extraction_phase3.ps1) first to generate validation data."
    exit 1
}

$validatedData = Get-Content $InputJson -Raw | ConvertFrom-Json

# ===================================================================
# C++ Class Generation Functions
# ===================================================================

function Get-SafeCppIdentifier {
    param([string]$Name)
    
    # Convert to safe C++ identifier
    $safe = $Name -replace '[^a-zA-Z0-9_]', '_'
    $safe = $safe -replace '^(\d)', '_$1'  # Prefix digit with underscore
    $safe = $safe -replace '__+', '_'      # Collapse multiple underscores
    return $safe.Trim('_')
}

function Get-CategoryEnum {
    param([string]$Type)
    
    switch ($Type.ToLower()) {
        "constant" { return "PhysicsTermCategory::FUNDAMENTAL_CONSTANTS" }
        "particle" { return "PhysicsTermCategory::PARTICLE_PHYSICS" }
        "system"   { return "PhysicsTermCategory::ASTROPHYSICAL_SYSTEMS" }
        "quantity" { return "PhysicsTermCategory::QUANTUM_MECHANICS" }
        default    { return "PhysicsTermCategory::OTHER" }
    }
}

function Generate-ConstantClass {
    param(
        [string]$Name,
        [object]$Entity
    )
    
    $className = "Wolfram_$(Get-SafeCppIdentifier $Name)"
    $category = Get-CategoryEnum "constant"
    $value = if ($Entity.value) { $Entity.value } else { "0.0" }
    $unit = if ($Entity.unit) { $Entity.unit } else { "dimensionless" }
    $description = if ($Entity.description) { $Entity.description } else { $Name }
    
    # Escape special characters in description
    $description = $description -replace '\\', '\\\\'
    $description = $description -replace '"', '\"'
    
    return @"
// Wolfram Knowledgebase: $Name
class $className : public PhysicsTerm {
public:
    $className() : PhysicsTerm("$Name", $category) {}
    
    double calculate(const SystemParams& params) const override {
        // Canonical value from Wolfram Knowledgebase
        return $value;
    }
    
    std::string getDescription() const override {
        return "$description [$unit]";
    }
    
    bool validate() const override {
        // Wolfram-validated constant
        return true;
    }
};

"@
}

function Generate-SystemClass {
    param(
        [string]$Name,
        [object]$Entity
    )
    
    $className = "Wolfram_$(Get-SafeCppIdentifier $Name)"
    $category = Get-CategoryEnum "system"
    $description = if ($Entity.description) { $Entity.description } else { $Name }
    
    # Escape special characters
    $description = $description -replace '\\', '\\\\'
    $description = $description -replace '"', '\"'
    
    # Extract properties for system modeling
    $mass = if ($Entity.mass) { $Entity.mass } else { "params.M_ns" }
    $radius = if ($Entity.radius) { $Entity.radius } else { "params.R_ns" }
    
    return @"
// Wolfram Knowledgebase: $Name
class $className : public PhysicsTerm {
public:
    $className() : PhysicsTerm("$Name", $category) {}
    
    double calculate(const SystemParams& params) const override {
        // Astrophysical system modeling from Wolfram data
        const double M = $mass;
        const double R = $radius;
        
        // Basic gravitational potential energy
        const double G = 6.67430e-11;  // m^3 kg^-1 s^-2
        return -G * M * M / R;
    }
    
    std::string getDescription() const override {
        return "$description";
    }
    
    bool validate() const override {
        return true;
    }
};

"@
}

function Generate-RegistrationFunction {
    param(
        [string]$SourceFile,
        [array]$ClassNames
    )
    
    $funcName = "registerWolframTerms_$(Get-SafeCppIdentifier $SourceFile)"
    
    $registrations = $ClassNames | ForEach-Object {
        "        registry.registerTerm(std::make_unique<$_>());"
    }
    
    return @"

// ===================================================================
// Registration function for $SourceFile Wolfram terms
// ===================================================================
void $funcName(PhysicsTermRegistry& registry) {
$($registrations -join "`n")
    
    std::cout << "Registered " << $($ClassNames.Count) << " Wolfram terms from $SourceFile" << std::endl;
}
"@
}

# ===================================================================
# Main Generation Logic
# ===================================================================

Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Phase 4: C++ Class Generation" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host ""

$totalClasses = 0
$fileStats = @{}

# Group entities by source file
$entityGroups = @{}

# Process constants
if ($validatedData.validatedConstants) {
    foreach ($const in $validatedData.validatedConstants) {
        $sourceFile = if ($const.sourceFile) { $const.sourceFile } else { "general" }
        
        if (-not $entityGroups.ContainsKey($sourceFile)) {
            $entityGroups[$sourceFile] = @{
                constants = @()
                systems = @()
            }
        }
        
        $entityGroups[$sourceFile].constants += $const
    }
}

# Process systems
if ($validatedData.validatedSystems) {
    foreach ($system in $validatedData.validatedSystems) {
        $sourceFile = if ($system.sourceFile) { $system.sourceFile } else { "general" }
        
        if (-not $entityGroups.ContainsKey($sourceFile)) {
            $entityGroups[$sourceFile] = @{
                constants = @()
                systems = @()
            }
        }
        
        $entityGroups[$sourceFile].systems += $system
    }
}

Write-Host "Entity groups found: $($entityGroups.Keys.Count)" -ForegroundColor Green
Write-Host ""

# Generate C++ files for each source group
foreach ($sourceFile in $entityGroups.Keys | Sort-Object) {
    $group = $entityGroups[$sourceFile]
    $constantCount = $group.constants.Count
    $systemCount = $group.systems.Count
    $totalCount = $constantCount + $systemCount
    
    if ($totalCount -eq 0) { continue }
    
    Write-Host "Processing $sourceFile ($constantCount constants, $systemCount systems)..." -ForegroundColor Yellow
    
    # Generate output filename
    $outputFile = if ($sourceFile -eq "general") {
        "wolfram_general.cpp"
    } else {
        "${sourceFile}_wolfram.cpp"
    }
    
    $outputPath = Join-Path $OutputDir $outputFile
    
    # Build C++ file content
    $cppContent = @"
// ===================================================================
// Auto-generated Wolfram Physics Terms from $sourceFile
// Generated by: wolfram_extraction_phase4.ps1
// Date: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
// ===================================================================
// Source: $sourceFile
// Constants: $constantCount
// Systems: $systemCount
// Total Terms: $totalCount
// ===================================================================

#include <memory>
#include <string>
#include <iostream>

// Forward declarations (assume these exist in main codebase)
class PhysicsTerm;
class PhysicsTermRegistry;
struct SystemParams;
enum class PhysicsTermCategory;

// ===================================================================
// Constant Definitions
// ===================================================================

"@
    
    # Generate constant classes
    $classNames = @()
    foreach ($const in $group.constants) {
        $className = "Wolfram_$(Get-SafeCppIdentifier $const.name)"
        $classNames += $className
        $cppContent += Generate-ConstantClass -Name $const.name -Entity $const
    }
    
    $cppContent += @"

// ===================================================================
// Astrophysical System Definitions
// ===================================================================

"@
    
    # Generate system classes
    foreach ($system in $group.systems) {
        $className = "Wolfram_$(Get-SafeCppIdentifier $system.name)"
        $classNames += $className
        $cppContent += Generate-SystemClass -Name $system.name -Entity $system
    }
    
    # Generate registration function
    $cppContent += Generate-RegistrationFunction -SourceFile $sourceFile -ClassNames $classNames
    
    # Save to file
    if ($DryRun) {
        Write-Host "  [DRY RUN] Would write: $outputPath" -ForegroundColor Gray
    } else {
        Set-Content -Path $outputPath -Value $cppContent -Encoding UTF8
        Write-Host "  ✓ Generated: $outputPath" -ForegroundColor Green
    }
    
    $fileStats[$sourceFile] = @{
        file = $outputFile
        constants = $constantCount
        systems = $systemCount
        total = $totalCount
        classes = $classNames
    }
    
    $totalClasses += $totalCount
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Generation Summary" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Total C++ files generated: $($fileStats.Keys.Count)" -ForegroundColor Green
Write-Host "Total PhysicsTerm classes: $totalClasses" -ForegroundColor Green
Write-Host ""

# Generate summary report
$summaryPath = Join-Path $OutputDir "generation_summary.json"
$summaryData = @{
    timestamp = Get-Date -Format "o"
    totalFiles = $fileStats.Keys.Count
    totalClasses = $totalClasses
    fileDetails = $fileStats
    sourceGroups = @($entityGroups.Keys | Sort-Object)
}

$summaryData | ConvertTo-Json -Depth 10 | Set-Content $summaryPath
Write-Host "Summary saved to: $summaryPath" -ForegroundColor Cyan

# Generate master registration header
$masterHeaderPath = Join-Path $OutputDir "wolfram_master_registration.h"

if (-not $DryRun) {
    # Build content with placeholders to avoid PowerShell parsing issues
    $content = @'
// ===================================================================
// Master Registration Header for All Wolfram Physics Terms
// Auto-generated by: wolfram_extraction_phase4.ps1
// Date: {TIMESTAMP}
// ===================================================================

#ifndef WOLFRAM_MASTER_REGISTRATION_H
#define WOLFRAM_MASTER_REGISTRATION_H

#include <iostream>

// Forward declaration
class PhysicsTermRegistry;

// Registration function declarations
{FUNCTION_DECLARATIONS}

// ===================================================================
// Master registration function
// Call this to register ALL Wolfram terms
// ===================================================================
inline void registerAllExtractedWolframTerms(PhysicsTermRegistry{AMP} registry) {
    std::cout {LT}{LT} "Registering Wolfram Physics Terms..." {LT}{LT} std::endl;
    
{FUNCTION_CALLS}
    
    std::cout {LT}{LT} "All Wolfram terms registered successfully." {LT}{LT} std::endl;
}

#endif // WOLFRAM_MASTER_REGISTRATION_H
'@
    
    # Generate function declarations
    $declarations = @()
    foreach ($sourceFile in $fileStats.Keys | Sort-Object) {
        $funcName = "registerWolframTerms_$(Get-SafeCppIdentifier $sourceFile)"
        $declarations += "extern void $funcName(PhysicsTermRegistry{AMP} registry);"
    }
    
    # Generate function calls
    $calls = @()
    foreach ($sourceFile in $fileStats.Keys | Sort-Object) {
        $funcName = "registerWolframTerms_$(Get-SafeCppIdentifier $sourceFile)"
        $calls += "    $funcName(registry);"
    }
    
    # Replace placeholders
    $content = $content -replace '\{TIMESTAMP\}', (Get-Date -Format 'yyyy-MM-dd HH:mm:ss')
    $content = $content -replace '\{FUNCTION_DECLARATIONS\}', ($declarations -join "`n")
    $content = $content -replace '\{FUNCTION_CALLS\}', ($calls -join "`n")
    $content = $content -replace '\{AMP\}', '&'
    $content = $content -replace '\{LT\}', '<'
    
    Set-Content -Path $masterHeaderPath -Value $content -Encoding UTF8
    Write-Host "Master registration header: $masterHeaderPath" -ForegroundColor Green
}

Write-Host ""
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "Next Steps" -ForegroundColor Cyan
Write-Host "=====================================" -ForegroundColor Cyan
Write-Host "1. Review generated files in: $OutputDir" -ForegroundColor Yellow
Write-Host "2. Add files to CMakeLists.txt target_sources" -ForegroundColor Yellow
Write-Host "3. Include wolfram_master_registration.h in MAIN_1_CoAnQi.cpp" -ForegroundColor Yellow
Write-Host "4. Call registerAllExtractedWolframTerms in main function" -ForegroundColor Yellow
Write-Host "5. Rebuild and verify 6597/6597 registration success" -ForegroundColor Yellow
Write-Host ""
