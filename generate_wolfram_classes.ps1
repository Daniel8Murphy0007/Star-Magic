# Generate C++ PhysicsTerm class definitions from Wolfram entities
# Input: wolfram_entities_parsed.json
# Output: wolfram_physics_classes.cpp (class definitions + registration calls)

$data = Get-Content "wolfram_entities_parsed.json" -Raw | ConvertFrom-Json

$classDefinitions = @()
$registrationCalls = @()
$totalClasses = 0

Write-Host "=== GENERATING C++ CLASSES ===" -ForegroundColor Cyan

# Helper function to sanitize names for C++
function Get-SafeCppName {
    param($name)
    $safe = $name -replace '[^a-zA-Z0-9_]', '_'
    $safe = $safe -replace '^(\d)', '_$1'  # Can't start with digit
    return $safe
}

# Generate PhysicalConstant classes
Write-Host "Generating PhysicalConstant classes..." -ForegroundColor Yellow
foreach ($name in $data.entities.PhysicalConstant) {
    $className = "$(Get-SafeCppName $name)ConstantTerm"
    $totalClasses++
    
    $classDefinitions += @"
class $className : public PhysicsTerm {
private:
    double constantValue;
public:
    $className() : constantValue(0.0) {
        setMetadata("source", "Wolfram-PhysicalConstant");
        setMetadata("wolfram_entity", "$name");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Physical constant - time-independent value
        // TODO: Query Wolfram for actual value via QuantityMagnitude
        return constantValue;
    }
    
    std::string getName() const override { return "$name"; }
    std::string getDescription() const override { 
        return "Wolfram Physical Constant: $name"; 
    }
};

"@

    $registrationCalls += "    core.registerPhysicsTerm(""$name"", std::make_unique<$className>(), ""Wolfram-PhysicalConstant"");"
}

# Generate Particle classes
Write-Host "Generating Particle classes..." -ForegroundColor Yellow
foreach ($name in $data.entities.Particle) {
    $className = "$(Get-SafeCppName $name)ParticleTerm"
    $totalClasses++
    
    $classDefinitions += @"
class $className : public PhysicsTerm {
private:
    double restMassEnergy;
public:
    $className() : restMassEnergy(0.0) {
        setMetadata("source", "Wolfram-Particle");
        setMetadata("wolfram_entity", "$name");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // E = m*c^2 for rest mass energy
        // TODO: Query Wolfram for particle mass
        const double c = 2.998e8;  // m/s
        return restMassEnergy * c * c;
    }
    
    std::string getName() const override { return "$name"; }
    std::string getDescription() const override { 
        return "Wolfram Particle: $name"; 
    }
};

"@

    $registrationCalls += "    core.registerPhysicsTerm(""$name"", std::make_unique<$className>(), ""Wolfram-Particle"");"
}

# Generate Isotope classes
Write-Host "Generating Isotope classes..." -ForegroundColor Yellow
foreach ($name in $data.entities.Isotope) {
    $className = "$(Get-SafeCppName $name)IsotopeTerm"
    $totalClasses++
    
    $classDefinitions += @"
class $className : public PhysicsTerm {
private:
    int atomicNumber;
    int massNumber;
public:
    $className() : atomicNumber(0), massNumber(0) {
        setMetadata("source", "Wolfram-Isotope");
        setMetadata("wolfram_entity", "$name");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Nuclear binding energy approximation (semi-empirical mass formula)
        // TODO: Query Wolfram for actual binding energy
        if (massNumber == 0) return 0.0;
        
        const double a_v = 15.75;   // Volume term (MeV)
        const double a_s = 17.8;    // Surface term (MeV)
        double bindingEnergy = a_v * massNumber - a_s * std::pow(massNumber, 2.0/3.0);
        return bindingEnergy * 1e6 * 1.602e-19;  // Convert MeV to Joules
    }
    
    std::string getName() const override { return "$name"; }
    std::string getDescription() const override { 
        return "Wolfram Isotope: $name"; 
    }
};

"@

    $registrationCalls += "    core.registerPhysicsTerm(""$name"", std::make_unique<$className>(), ""Wolfram-Isotope"");"
}

# Generate PhysicalQuantity classes
Write-Host "Generating PhysicalQuantity classes..." -ForegroundColor Yellow
foreach ($name in $data.entities.PhysicalQuantity) {
    $className = "$(Get-SafeCppName $name)QuantityTerm"
    $totalClasses++
    
    $classDefinitions += @"
class $className : public PhysicsTerm {
private:
    double quantityValue;
public:
    $className() : quantityValue(0.0) {
        setMetadata("source", "Wolfram-PhysicalQuantity");
        setMetadata("wolfram_entity", "$name");
    }
    
    double compute(double t, const std::map<std::string, double>& params) const override {
        // Physical quantity measurement
        // TODO: Query Wolfram for quantity definition
        return quantityValue;
    }
    
    std::string getName() const override { return "$name"; }
    std::string getDescription() const override { 
        return "Wolfram Physical Quantity: $name"; 
    }
};

"@

    $registrationCalls += "    core.registerPhysicsTerm(""$name"", std::make_unique<$className>(), ""Wolfram-PhysicalQuantity"");"
}

# Output results
Write-Host "`nGenerated $totalClasses class definitions" -ForegroundColor Green
Write-Host "Generated $($registrationCalls.Count) registration calls" -ForegroundColor Green

# Create output file
$outputContent = @"
// ============================================================================
// WOLFRAM PHYSICS TERMS - AUTO-GENERATED FROM WOLFRAM KNOWLEDGEBASE
// Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
// Total Classes: $totalClasses
// Source: wolfram_physics_terms_FULL.txt (5,703 entities)
// ============================================================================

// NOTE: These classes provide structural integration of Wolfram entities.
// TODO: Replace placeholder compute() methods with actual Wolfram API queries
// using WolframLibrary or WSTP for real-time physics data retrieval.

#include <cmath>
#include <string>
#include <map>
#include <memory>

// ============================================================================
// CLASS DEFINITIONS
// ============================================================================

$($classDefinitions -join "`n")

// ============================================================================
// REGISTRATION FUNCTION
// ============================================================================

void registerAllWolframPhysicsTerms(CalculatorCore& core) {
$($registrationCalls -join "`n")
}

// ============================================================================
// END OF AUTO-GENERATED WOLFRAM PHYSICS TERMS
// Total: $totalClasses classes registered
// ============================================================================
"@

$outputContent | Out-File "wolfram_physics_classes.cpp" -Encoding UTF8

Write-Host "`nOutput file: wolfram_physics_classes.cpp" -ForegroundColor Cyan
Write-Host "File size: $('{0:N2}' -f ((Get-Item 'wolfram_physics_classes.cpp').Length / 1MB)) MB" -ForegroundColor Yellow
