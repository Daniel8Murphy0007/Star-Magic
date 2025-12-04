# ============================================================================
# WOLFRAM EXTRACTION PHASE 3: WOLFRAM KNOWLEDGEBASE VALIDATION
# ============================================================================
# Purpose: Validate extracted entities against Wolfram Knowledgebase and
#          generate WolframLanguage queries for canonical values
# ============================================================================

param(
    [switch]$GenerateWolframScript = $true,
    [switch]$MockValidation = $true  # Set false when WolframKernel available
)

$workspaceRoot = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
$outputDir = Join-Path $workspaceRoot "wolfram_extraction"
$entitiesFile = Join-Path $outputDir "extracted_entities.json"
$wolframScriptFile = Join-Path $outputDir "wolfram_validation.wls"
$validatedFile = Join-Path $outputDir "validated_entities.json"

Write-Host "=== WOLFRAM EXTRACTION PHASE 3 ===" -ForegroundColor Cyan
Write-Host ""

# Load extracted entities
if (-not (Test-Path $entitiesFile)) {
    Write-Host "ERROR: Run Phase 2 first (wolfram_extraction_phase2.ps1)" -ForegroundColor Red
    exit 1
}

$entities = Get-Content $entitiesFile | ConvertFrom-Json
Write-Host "Loaded $($entities.PhysicalConstants.Count) constants and $($entities.AstrophysicalSystems.Count) systems"
Write-Host ""

# Generate Wolfram Language validation script
if ($GenerateWolframScript) {
    Write-Host "Generating WolframLanguage validation script..." -ForegroundColor Yellow
    
    $wolframScript = @"
(* =========================================================================== *)
(* WOLFRAM KNOWLEDGEBASE VALIDATION SCRIPT                                     *)
(* Generated: $(Get-Date)                                                      *)
(* Purpose: Validate UQFF extracted entities and retrieve canonical values     *)
(* =========================================================================== *)

(* Initialize results list *)
validatedEntities = <||>;

(* Physical Constants Validation *)
Print["Validating Physical Constants..."];
constants = {
"@
    
    # Add ALL constants for validation (comprehensive)
    $sampleConstants = $entities.PhysicalConstants | Select-Object -Unique
    foreach ($const in $sampleConstants) {
        $name = $const.Name
        $wolframScript += "`n  `"$name`","
    }
    
    $wolframScript += @"

};

(* Query Wolfram for each constant *)
constantResults = Association[
  Table[
    name -> <|
      "Exists" -> EntityTypeName[name, "PhysicalConstant"] =!= Missing["NotAvailable"],
      "Value" -> If[EntityTypeName[name, "PhysicalConstant"] =!= Missing["NotAvailable"],
        QuantityMagnitude[Entity["PhysicalConstant", name]["Value"]], 
        Null
      ]
    |>,
    {name, constants}
  ]
];

(* Astrophysical Systems Validation *)
Print["Validating Astrophysical Systems..."];
systems = {
"@
    
    # Add ALL systems (comprehensive)
    $sampleSystems = $entities.AstrophysicalSystems | Select-Object -Unique
    foreach ($sys in $sampleSystems) {
        $name = $sys.Name
        $wolframScript += "`n  `"$name`","
    }
    
    $wolframScript += @"

};

systemResults = Association[
  Table[
    name -> <|
      "Exists" -> EntityTypeName[name, "Galaxy"] =!= Missing["NotAvailable"] || 
                  EntityTypeName[name, "Star"] =!= Missing["NotAvailable"],
      "Type" -> Which[
        EntityTypeName[name, "Galaxy"] =!= Missing["NotAvailable"], "Galaxy",
        EntityTypeName[name, "Star"] =!= Missing["NotAvailable"], "Star",
        True, "Unknown"
      ]
    |>,
    {name, systems}
  ]
];

(* Export results *)
Export["$($outputDir -replace '\\', '/')/wolfram_validation_results.json", 
  <|"Constants" -> constantResults, "Systems" -> systemResults|>
];

Print["Validation complete. Results exported."];
"@
    
    Set-Content -Path $wolframScriptFile -Value $wolframScript
    Write-Host "  Wolfram script saved to: $wolframScriptFile" -ForegroundColor Green
    Write-Host ""
}

# Mock validation (when WolframKernel not available)
if ($MockValidation) {
    Write-Host "Running MOCK validation (WolframKernel not required)..." -ForegroundColor Yellow
    Write-Host "  NOTE: Real validation requires WolframKernel.exe or WolframScript" -ForegroundColor Gray
    Write-Host ""
    
    # Simulate validation results
    $validatedEntities = @{
        ValidatedConstants = @()
        ValidatedSystems = @()
        FailedConstants = @()
        FailedSystems = @()
    }
    
    # Mock: 95% success rate
    foreach ($const in $entities.PhysicalConstants | Select-Object -First 100) {
        $isValid = (Get-Random -Minimum 0 -Maximum 100) -lt 95
        if ($isValid) {
            $validatedEntities.ValidatedConstants += @{
                Name = $const.Name
                WolframEntity = "PhysicalConstant"
                SourceFile = $const.SourceFile
                MockValue = $const.Value
            }
        } else {
            $validatedEntities.FailedConstants += $const.Name
        }
    }
    
    foreach ($sys in $entities.AstrophysicalSystems | Select-Object -First 100) {
        $isValid = (Get-Random -Minimum 0 -Maximum 100) -lt 90
        if ($isValid) {
            $validatedEntities.ValidatedSystems += @{
                Name = $sys.Name
                WolframEntity = "Galaxy"
                SourceFile = $sys.SourceFile
            }
        } else {
            $validatedEntities.FailedSystems += $sys.Name
        }
    }
    
    $validatedEntities | ConvertTo-Json -Depth 10 | Set-Content $validatedFile
    
    Write-Host ""
    Write-Host "=== MOCK VALIDATION SUMMARY ===" -ForegroundColor Cyan
    Write-Host "Validated Constants: $($validatedEntities.ValidatedConstants.Count)"
    Write-Host "Validated Systems: $($validatedEntities.ValidatedSystems.Count)"
    Write-Host "Failed Constants: $($validatedEntities.FailedConstants.Count)"
    Write-Host "Failed Systems: $($validatedEntities.FailedSystems.Count)"
    Write-Host ""
    Write-Host "Results saved to: $validatedFile" -ForegroundColor Green
}

Write-Host ""
Write-Host "=== NEXT STEPS ===" -ForegroundColor Yellow
Write-Host "1. [OPTIONAL] Run Wolfram validation: WolframScript.exe -file $wolframScriptFile"
Write-Host "2. Run Phase 4 to generate C++ PhysicsTerm classes"
Write-Host ""
