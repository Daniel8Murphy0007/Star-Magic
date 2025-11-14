# Smart Extraction Script - Phase 1 Batch Processing
# Purpose: Extract 157 modules into organized Core/ structure
# Strategy: Process by priority and complexity, not just category

param(
    [Parameter(Mandatory=$false)]
    [string]$Phase = "1",  # 1=Foundation, 2=Simple, 3=Moderate, 4=Complex
    
    [Parameter(Mandatory=$false)]
    [switch]$DryRun,
    
    [Parameter(Mandatory=$false)]
    [int]$BatchSize = 20
)

# Load analysis results
$analysisPath = Get-ChildItem -Filter "module_analysis_*.json" | Sort-Object LastWriteTime -Descending | Select-Object -First 1
if (-not $analysisPath) {
    Write-Error "No analysis file found. Run: .\batch_extract.ps1 -Mode analyze"
    exit 1
}

$analysis = Get-Content $analysisPath.FullName | ConvertFrom-Json
Write-Host "Loaded analysis: $($analysis.TotalModules) modules" -ForegroundColor Green

# Define extraction phases
$script:extractionPlan = @{
    "Phase1_Foundation" = @{
        Description = "Critical foundation modules (already extracted in Weeks 1-2)"
        Modules = @(10, 4)  # SystemCatalogue, UQFFModule4
        Priority = 1
        Status = "COMPLETE"
    }
    "Phase2_SimpleEnhanced" = @{
        Description = "Simple enhanced modules (<500 lines, already have framework)"
        Modules = @()  # Will be populated
        Priority = 2
        Status = "PENDING"
    }
    "Phase3_ModerateEnhanced" = @{
        Description = "Moderate enhanced modules (500-2000 lines)"
        Modules = @()
        Priority = 3
        Status = "PENDING"
    }
    "Phase4_ComplexEnhanced" = @{
        Description = "Complex enhanced modules (2000-5000 lines)"
        Modules = @()
        Priority = 4
        Status = "PENDING"
    }
    "Phase5_NonEnhanced" = @{
        Description = "Modules without enhancement framework (need manual review)"
        Modules = @()
        Priority = 5
        Status = "PENDING"
    }
}

# Categorize modules into phases
foreach ($module in $analysis.Modules) {
    $num = $module.Number
    
    # Skip already processed
    if ($num -in @(10, 4)) { continue }
    
    if ($module.HasEnhancedFramework) {
        switch ($module.EstimatedComplexity) {
            "Simple" { $script:extractionPlan.Phase2_SimpleEnhanced.Modules += $num }
            "Moderate" { $script:extractionPlan.Phase3_ModerateEnhanced.Modules += $num }
            "Complex" { $script:extractionPlan.Phase4_ComplexEnhanced.Modules += $num }
            "Very Complex" { $script:extractionPlan.Phase4_ComplexEnhanced.Modules += $num }
        }
    } else {
        $script:extractionPlan.Phase5_NonEnhanced.Modules += $num
    }
}

# Display extraction plan
Write-Host "`n========== SMART EXTRACTION PLAN ==========" -ForegroundColor Cyan
foreach ($phaseName in $script:extractionPlan.Keys | Sort-Object) {
    $p = $script:extractionPlan[$phaseName]
    $count = $p.Modules.Count
    Write-Host "$phaseName [$($p.Status)]" -ForegroundColor Yellow
    Write-Host "  Description: $($p.Description)" -ForegroundColor Gray
    Write-Host "  Modules: $count" -ForegroundColor Gray
    if ($count -gt 0 -and $count -le 20) {
        Write-Host "  Numbers: $($p.Modules -join ', ')" -ForegroundColor DarkGray
    }
}
Write-Host "==========================================`n" -ForegroundColor Cyan

# Function to extract a single module intelligently
function Extract-ModuleIntelligent {
    param(
        [int]$ModuleNumber,
        [string]$OutputDir
    )
    
    # Find module file
    $lowerPath = "source$ModuleNumber.cpp"
    $upperPath = "Source$ModuleNumber.cpp"
    $sourcePath = if (Test-Path $lowerPath) { $lowerPath } 
                  elseif (Test-Path $upperPath) { $upperPath }
                  else { return $false }
    
    Write-Host "Processing Module $ModuleNumber..." -ForegroundColor Green
    
    if ($DryRun) {
        Write-Host "  [DRY RUN] Would extract $sourcePath" -ForegroundColor Yellow
        return $true
    }
    
    # Read source
    $content = Get-Content $sourcePath -Raw
    
    # Find all class definitions
    $allMatches = [regex]::Matches($content, 'class\s+([A-Za-z_]\w*)')
    $classNames = $allMatches | ForEach-Object { $_.Groups[1].Value }
    
    # Filter out base framework classes
    $baseClasses = @('PhysicsTerm', 'DynamicVacuumTerm', 'QuantumCouplingTerm')
    $mainClasses = $classNames | Where-Object { $_ -notin $baseClasses }
    
    # Prioritize UQFF/Module classes, then take first non-base class
    $className = $null
    if ($mainClasses.Count -gt 0) {
        # Try to find best match
        $className = $mainClasses | Where-Object { $_ -match 'UQFF.*Module' } | Select-Object -First 1
        if (-not $className) {
            $className = $mainClasses | Where-Object { $_ -match 'Module' } | Select-Object -First 1
        }
        if (-not $className) {
            $className = $mainClasses | Where-Object { $_ -match 'Solver|Engine' } | Select-Object -First 1
        }
        if (-not $className) {
            # Take first non-base class
            $className = $mainClasses[0]
        }
    }
    
    if ($className) {
        Write-Host "  Found class: $className" -ForegroundColor Cyan
        
        # Create Core directory if needed
        if (-not (Test-Path $OutputDir)) {
            New-Item -ItemType Directory -Path $OutputDir -Force | Out-Null
        }
        
        # For now, copy entire module to Core/Modules/
        $targetPath = Join-Path $OutputDir "$className.cpp"
        Copy-Item $sourcePath $targetPath -Force
        Write-Host "  Copied to: $targetPath" -ForegroundColor Green
        
        # TODO: Smart parsing and header/implementation split
        # This is a placeholder - full extraction would:
        # 1. Parse class definitions
        # 2. Extract to .hpp/.cpp pair
        # 3. Update includes
        # 4. Generate CMakeLists.txt
        
        return $true
    } else {
        Write-Host "  WARNING: No recognizable class found" -ForegroundColor Red
        return $false
    }
}

# Execute extraction based on phase
function Execute-Phase {
    param([string]$PhaseName)
    
    $phaseData = $script:extractionPlan[$PhaseName]
    
    if ($null -eq $phaseData) {
        Write-Error "Phase not found in extraction plan: $PhaseName"
        return
    }
    
    $modules = $phaseData.Modules
    
    if ($modules.Count -eq 0) {
        Write-Host "No modules in $PhaseName" -ForegroundColor Yellow
        return
    }
    
    Write-Host "`nExecuting $PhaseName - $($modules.Count) modules" -ForegroundColor Cyan
    Write-Host "Description: $($phaseData.Description)`n" -ForegroundColor Gray
    
    # Create output directory
    $outputDir = Join-Path $PSScriptRoot "Core\Modules"
    
    $success = 0
    $failed = 0
    
    # Process in batches
    for ($i = 0; $i -lt $modules.Count; $i += $BatchSize) {
        $batch = $modules[$i..([Math]::Min($i + $BatchSize - 1, $modules.Count - 1))]
        Write-Host "Batch $([Math]::Floor($i / $BatchSize) + 1): Processing modules $($batch[0])-$($batch[-1])" -ForegroundColor Yellow
        
        foreach ($moduleNum in $batch) {
            if (Extract-ModuleIntelligent -ModuleNumber $moduleNum -OutputDir $outputDir) {
                $success++
            } else {
                $failed++
            }
        }
        
        Write-Host "  Batch complete: $success succeeded, $failed failed`n" -ForegroundColor Gray
    }
    
    Write-Host "$PhaseName COMPLETE: $success/$($modules.Count) modules extracted`n" -ForegroundColor Green
}

# Main execution
$phaseMap = @{
    "1" = "Phase1_Foundation"
    "2" = "Phase2_SimpleEnhanced"
    "3" = "Phase3_ModerateEnhanced"
    "4" = "Phase4_ComplexEnhanced"
    "5" = "Phase5_NonEnhanced"
}

switch ($Phase) {
    "1" { Execute-Phase $phaseMap["1"] }
    "2" { Execute-Phase $phaseMap["2"] }
    "3" { Execute-Phase $phaseMap["3"] }
    "4" { Execute-Phase $phaseMap["4"] }
    "5" { Execute-Phase $phaseMap["5"] }
    "all" {
        foreach ($p in @("2", "3", "4", "5")) {
            Execute-Phase $phaseMap[$p]
        }
    }
    default {
        Write-Error "Unknown phase: $Phase. Valid options: 1, 2, 3, 4, 5, all"
    }
}

Write-Host "`n✅ Extraction script complete!" -ForegroundColor Green
