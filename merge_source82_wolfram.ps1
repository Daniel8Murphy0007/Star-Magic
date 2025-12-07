# PowerShell script to merge source82_wolfram.cpp files
# C++20/MSVC 19.44+ compliant

Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host "SOURCE82_WOLFRAM.CPP MERGE UTILITY (PowerShell)" -ForegroundColor Cyan
Write-Host "C++20/MSVC 19.44+" -ForegroundColor Cyan
Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host ""

# Read both files
Write-Host "[1/6] Reading SMBH workspace version..." -ForegroundColor Yellow
$smbh_content = Get-Content "source82_wolfram.cpp" -Raw
Write-Host "      ✓ Loaded $($smbh_content.Length) bytes" -ForegroundColor Green

Write-Host "[2/6] Reading Virgo Cluster version..." -ForegroundColor Yellow
$virgo_content = Get-Content "source82_wolfram_VIRGO_EXTRACTION.cpp" -Raw -Encoding UTF8
Write-Host "      ✓ Loaded $($virgo_content.Length) bytes" -ForegroundColor Green

# Extract Virgo classes (skip SMBHMSigmaRelationTerm - keep SMBH version)
Write-Host "[3/6] Extracting Virgo classes (skipping duplicate SMBHMSigmaRelationTerm)..." -ForegroundColor Yellow

# Find the start and end of each Virgo class
$virgo_classes = @()
$class_names = @(
    'VirgoClusterMassTerm',
    'VirgoClusterIntraclusterMediumTerm',
    'VirgoClusterGravitationalPotentialTerm',
    'VirgoClusterDarkMatterTerm',
    'VirgoClusterM87JetTerm',
    'VirgoClusterTidalStrippingTerm',
    'VirgoClusterVirialTerm',
    'VirgoClusterXRayLuminosityTerm',
    'VirgoClusterVelocityDispersionTerm'
    # Skip SMBHMSigmaRelationTerm - duplicate
)

foreach ($class_name in $class_names) {
    # Find class definition with comments
    $pattern = "(?s)(//\s*=+\s*`n//\s*CLASS\s+\d+:\s*$class_name.*?`n//\s*=+.*?)class\s+$class_name\s*\{.*?\};`n"
    $match = [regex]::Match($virgo_content, $pattern)
    
    if ($match.Success) {
        $virgo_classes += $match.Value
        Write-Host "      ✓ Extracted: $class_name" -ForegroundColor Gray
    }
    else {
        Write-Host "      ⚠ Not found: $class_name" -ForegroundColor Yellow
    }
}

Write-Host "      Total Virgo classes extracted: $($virgo_classes.Count)/9" -ForegroundColor Green

# Create merged file header
Write-Host "[4/6] Generating unified C++20-compliant header..." -ForegroundColor Yellow

$merged_content = @'
// ============================================================================
// UNIFIED SOURCE82_WOLFRAM.CPP - DUAL-SCALE PHYSICS
// ============================================================================
// Merged: December 4, 2025
// Part 1: SMBH Local-Scale Physics (15 classes, 10^6-10^9 M_sun)
// Part 2: Virgo Cluster Cosmological-Scale Physics (9 unique classes, 10^15 M_sun)
// Total Classes: 24 (15 SMBH + 9 Virgo, 1 duplicate removed)
// Compiler: C++20 / MSVC 19.44.35207+
// Copyright: Daniel T. Murphy
// ============================================================================

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Physical constants (C++20 constexpr)
constexpr double G_CONST = 6.6743e-11;      // Gravitational constant (m³/kg·s²)
constexpr double M_SUN = 1.989e30;          // Solar mass (kg)
constexpr double MPC_TO_M = 3.086e22;       // Megaparsec to meters
constexpr double KPC_TO_M = 3.086e19;       // Kiloparsec to meters
constexpr double KEV_TO_J = 1.602e-16;      // keV to Joules
constexpr double K_BOLTZ = 1.381e-23;       // Boltzmann constant (J/K)
constexpr double M_PROTON = 1.673e-27;      // Proton mass (kg)
constexpr double C_LIGHT = 2.998e8;         // Speed of light (m/s)
constexpr double YEAR_TO_S = 3.156e7;       // Year to seconds

// ============================================================================
// BASE PHYSICS TERM INTERFACE (C++20)
// ============================================================================

class PhysicsTerm {
public:
    virtual ~PhysicsTerm() = default;  // C++20 defaulted virtual destructor
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const { return true; }
};

// ============================================================================
// PART 1: SMBH LOCAL-SCALE PHYSICS (15 CLASSES)
// Scale: 10^6 - 10^9 M_sun
// Physics: Black hole M-σ relation, vacuum energy, quantum coupling, UQFF layers
// ============================================================================

'@

# Extract SMBH classes from workspace version (skip header/comments before first class)
Write-Host "[5/6] Extracting SMBH classes from workspace..." -ForegroundColor Yellow
$smbh_classes_start = $smbh_content.IndexOf('class SMBHDynamicVacuumTerm')
$smbh_classes_end = $smbh_content.LastIndexOf('// End of source82_wolfram.cpp')
if ($smbh_classes_end -lt 0) {
    $smbh_classes_end = $smbh_content.Length
}

$smbh_classes_section = $smbh_content.Substring($smbh_classes_start, $smbh_classes_end - $smbh_classes_start).Trim()

# Combine all sections
$merged_content += $smbh_classes_section
$merged_content += @'


// ============================================================================
// PART 2: VIRGO CLUSTER COSMOLOGICAL-SCALE PHYSICS (9 UNIQUE CLASSES)
// Scale: 10^15 M_sun cluster mass
// Physics: Galaxy cluster dynamics, ICM, dark matter, X-ray emission
// Distance: 16.5 Mpc, Virial radius: 2.2 Mpc, Velocity dispersion: 700 km/s
// ============================================================================

'@

# Add Virgo classes
foreach ($virgo_class in $virgo_classes) {
    $merged_content += $virgo_class + "`n`n"
}

# Add registration function footer
$merged_content += @'

// ============================================================================
// UNIFIED REGISTRATION FUNCTION (C++20 std::make_unique)
// ============================================================================

/*
void registerWolframTerms_source82(PhysicsTermRegistry& registry) {
    // Part 1: SMBH Local-Scale Terms (15)
    registry.registerPhysicsTerm("SMBH_DynamicVacuum", std::make_unique<SMBHDynamicVacuumTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_QuantumCoupling", std::make_unique<SMBHQuantumCouplingTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_MSigmaRelation", std::make_unique<SMBHMSigmaRelationTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_BulgeGravity", std::make_unique<SMBHBulgeGravityTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Ug1", std::make_unique<SMBHUg1Term>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Ug2", std::make_unique<SMBHUg2Term>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Ug3", std::make_unique<SMBHUg3Term>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Ug4", std::make_unique<SMBHUg4Term>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_ReactorEfficiency", std::make_unique<SMBHReactorEfficiencyTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_PseudoMonopole", std::make_unique<SMBHPseudoMonopoleTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_RedshiftCorrection", std::make_unique<SMBHRedshiftCorrectionTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Ui", std::make_unique<SMBHUiTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_Um", std::make_unique<SMBHUmTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_OmegaS", std::make_unique<SMBHOmegaSGalacticTerm>(), "wolfram-smbh");
    registry.registerPhysicsTerm("SMBH_CosmicTime", std::make_unique<SMBHCosmicTimeTerm>(), "wolfram-smbh");
    
    // Part 2: Virgo Cluster Cosmological-Scale Terms (9)
    registry.registerPhysicsTerm("Virgo_ClusterMass", std::make_unique<VirgoClusterMassTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_ICM", std::make_unique<VirgoClusterIntraclusterMediumTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_Potential", std::make_unique<VirgoClusterGravitationalPotentialTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_DarkMatter", std::make_unique<VirgoClusterDarkMatterTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_M87Jet", std::make_unique<VirgoClusterM87JetTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_TidalStripping", std::make_unique<VirgoClusterTidalStrippingTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_Virial", std::make_unique<VirgoClusterVirialTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_XRay", std::make_unique<VirgoClusterXRayLuminosityTerm>(), "wolfram-virgo");
    registry.registerPhysicsTerm("Virgo_VelocityDispersion", std::make_unique<VirgoClusterVelocityDispersionTerm>(), "wolfram-virgo");
    
    std::cout << "Registered 24 Wolfram terms from source82_wolfram.cpp (15 SMBH + 9 Virgo)" << std::endl;
}
*/

// ============================================================================
// SUMMARY: 24 UNIFIED PHYSICS TERM CLASSES
// ============================================================================
// SMBH LOCAL-SCALE (15 classes):
//   1. SMBHDynamicVacuumTerm
//   2. SMBHQuantumCouplingTerm
//   3. SMBHMSigmaRelationTerm
//   4. SMBHBulgeGravityTerm
//   5. SMBHUg1Term
//   6. SMBHUg2Term
//   7. SMBHUg3Term
//   8. SMBHUg4Term
//   9. SMBHReactorEfficiencyTerm
//   10. SMBHPseudoMonopoleTerm
//   11. SMBHRedshiftCorrectionTerm
//   12. SMBHUiTerm
//   13. SMBHUmTerm
//   14. SMBHOmegaSGalacticTerm
//   15. SMBHCosmicTimeTerm
//
// VIRGO CLUSTER COSMOLOGICAL-SCALE (9 classes):
//   16. VirgoClusterMassTerm
//   17. VirgoClusterIntraclusterMediumTerm
//   18. VirgoClusterGravitationalPotentialTerm
//   19. VirgoClusterDarkMatterTerm
//   20. VirgoClusterM87JetTerm
//   21. VirgoClusterTidalStrippingTerm
//   22. VirgoClusterVirialTerm
//   23. VirgoClusterXRayLuminosityTerm
//   24. VirgoClusterVelocityDispersionTerm
//
// Duplicate removed: SMBHMSigmaRelationTerm (kept SMBH version from workspace)
// ============================================================================
// Integration: Include in MAIN_1_CoAnQi.cpp and uncomment registration function
// Build: MSVC 19.44+ with /std:c++20
// ============================================================================
'@

# Write merged file
Write-Host "[6/6] Writing merged source82_wolfram_MERGED.cpp..." -ForegroundColor Yellow
$merged_content | Out-File -FilePath "source82_wolfram_MERGED.cpp" -Encoding UTF8 -NoNewline
Write-Host "      ✓ Written $($merged_content.Length) bytes" -ForegroundColor Green

Write-Host ""
Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host "MERGE COMPLETE" -ForegroundColor Green
Write-Host "=" * 80 -ForegroundColor Cyan
Write-Host "Total classes: 24 (15 SMBH + 9 Virgo)" -ForegroundColor White
Write-Host "Output: source82_wolfram_MERGED.cpp" -ForegroundColor Yellow
Write-Host ""
Write-Host "Next steps:" -ForegroundColor Cyan
Write-Host "  1. Review source82_wolfram_MERGED.cpp" -ForegroundColor White
Write-Host "  2. Rename to source82_wolfram.cpp" -ForegroundColor White
Write-Host "  3. Build test: cmake --build build_msvc --config Release" -ForegroundColor White
Write-Host ""
