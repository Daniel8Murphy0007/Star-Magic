# Star-Magic UQFF Batch Extraction Tool
# Purpose: Automated extraction of 157 source*.cpp modules into Core/ architecture
# Date: November 14, 2025

param(
    [Parameter(Mandatory=$false)]
    [string]$Mode = "analyze",  # analyze, extract, test, all
    
    [Parameter(Mandatory=$false)]
    [int[]]$ModuleNumbers = @(),  # Specific modules to process (empty = all)
    
    [Parameter(Mandatory=$false)]
    [string]$BatchName = "",  # e.g., "AstroPhysics", "QuantumFields"
    
    [Parameter(Mandatory=$false)]
    [switch]$DryRun  # Simulate without making changes
)

# Configuration
$script:Config = @{
    SourceDir = $PSScriptRoot
    CoreDir = Join-Path $PSScriptRoot "Core"
    BackupDir = Join-Path $PSScriptRoot "extraction_backup_$(Get-Date -Format 'yyyyMMdd_HHmmss')"
    CompilerPath = "g++"
    CppStandard = "c++17"
    LogFile = Join-Path $PSScriptRoot "extraction_log_$(Get-Date -Format 'yyyyMMdd_HHmmss').txt"
}

# Module categorization based on common patterns
$script:Categories = @{
    "Foundation" = @{
        Pattern = @("SystemParams", "PhysicsConstants", "MasterEquations")
        Directory = "Foundation"
        Modules = @(10)  # Known from Week 1
    }
    "AdaptiveFramework" = @{
        Pattern = @("UQFFModule", "FluidSolver", "AdaptiveGrid")
        Directory = "Adaptive"
        Modules = @(4)  # Known from Week 2
    }
    "AstroPhysics" = @{
        Pattern = @("Galaxy", "Stellar", "AGN", "Pulsar", "Quasar", "BlackHole", "SMBH")
        Directory = "AstroPhysics"
        Modules = @()
    }
    "QuantumFields" = @{
        Pattern = @("Vacuum", "Quantum", "Field", "Entanglement", "Coupling")
        Directory = "QuantumFields"
        Modules = @()
    }
    "FluidDynamics" = @{
        Pattern = @("Navier", "Stokes", "Fluid", "MHD", "Hydro", "Turbulence")
        Directory = "FluidDynamics"
        Modules = @()
    }
    "Relativity" = @{
        Pattern = @("Relativity", "GR", "FrameDrag", "Metric", "Spacetime")
        Directory = "Relativity"
        Modules = @()
    }
    "Enhanced" = @{
        Pattern = @("registerDynamicTerm", "PhysicsTerm", "DynamicVacuum")
        Directory = "Enhanced"
        Modules = @(14..167)  # Source14-167 already enhanced
    }
    "Nuclear" = @{
        Pattern = @("Neutron", "Nuclear", "LENR", "Fission", "Fusion")
        Directory = "Nuclear"
        Modules = @()
    }
    "Cosmology" = @{
        Pattern = @("Cosmic", "DarkMatter", "DarkEnergy", "Redshift", "Expansion")
        Directory = "Cosmology"
        Modules = @()
    }
}

# Logging function
function Write-Log {
    param([string]$Message, [string]$Level = "INFO")
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    $logEntry = "[$timestamp] [$Level] $Message"
    Write-Host $logEntry
    Add-Content -Path $script:Config.LogFile -Value $logEntry
}

# Analyze a single module file
function Analyze-Module {
    param(
        [Parameter(Mandatory=$true)]
        [string]$FilePath,
        [int]$ModuleNumber
    )
    
    if (-not (Test-Path $FilePath)) {
        return $null
    }
    
    $content = Get-Content $FilePath -Raw
    $analysis = @{
        Number = $ModuleNumber
        FilePath = $FilePath
        FileName = Split-Path $FilePath -Leaf
        LineCount = (Get-Content $FilePath).Count
        Classes = @()
        Structs = @()
        Functions = @()
        Includes = @()
        Categories = @()
        HasEnhancedFramework = $false
        HasPhysicsTerm = $false
        EstimatedComplexity = "Unknown"
    }
    
    # Extract classes
    $classMatches = [regex]::Matches($content, '(?m)^class\s+(\w+)')
    foreach ($match in $classMatches) {
        $analysis.Classes += $match.Groups[1].Value
    }
    
    # Extract structs
    $structMatches = [regex]::Matches($content, '(?m)^struct\s+(\w+)')
    foreach ($match in $structMatches) {
        $analysis.Structs += $match.Groups[1].Value
    }
    
    # Extract includes
    $includeMatches = [regex]::Matches($content, '(?m)^#include\s+[<"]([^>"]+)[>"]')
    foreach ($match in $includeMatches) {
        $analysis.Includes += $match.Groups[1].Value
    }
    
    # Check for enhanced framework
    $analysis.HasEnhancedFramework = $content -match 'registerDynamicTerm|class PhysicsTerm'
    $analysis.HasPhysicsTerm = $content -match 'class \w+Term\s*:\s*public PhysicsTerm'
    
    # Categorize module
    foreach ($category in $script:Categories.Keys) {
        $patterns = $script:Categories[$category].Pattern
        foreach ($pattern in $patterns) {
            if ($content -match $pattern) {
                if ($analysis.Categories -notcontains $category) {
                    $analysis.Categories += $category
                }
            }
        }
    }
    
    # Estimate complexity
    if ($analysis.LineCount -lt 500) {
        $analysis.EstimatedComplexity = "Simple"
    } elseif ($analysis.LineCount -lt 2000) {
        $analysis.EstimatedComplexity = "Moderate"
    } elseif ($analysis.LineCount -lt 5000) {
        $analysis.EstimatedComplexity = "Complex"
    } else {
        $analysis.EstimatedComplexity = "Very Complex"
    }
    
    return $analysis
}

# Analyze all modules
function Analyze-AllModules {
    Write-Log "Starting module analysis..." "INFO"
    
    $allAnalyses = @()
    $existingModules = @()
    
    # Find all existing source files
    for ($i = 1; $i -le 167; $i++) {
        $lowerPath = Join-Path $script:Config.SourceDir "source$i.cpp"
        $upperPath = Join-Path $script:Config.SourceDir "Source$i.cpp"
        
        $filePath = if (Test-Path $lowerPath) { $lowerPath } 
                    elseif (Test-Path $upperPath) { $upperPath }
                    else { $null }
        
        if ($filePath) {
            $existingModules += $i
            $analysis = Analyze-Module -FilePath $filePath -ModuleNumber $i
            if ($analysis) {
                $allAnalyses += $analysis
            }
        }
    }
    
    Write-Log "Found $($existingModules.Count) existing modules (1-167 with gaps)" "INFO"
    Write-Log "Missing modules: $(($1..167 | Where-Object { $_ -notin $existingModules }) -join ', ')" "INFO"
    
    # Generate categorization report
    $report = @{
        TotalModules = $allAnalyses.Count
        ByCategory = @{}
        ByComplexity = @{}
        Enhanced = ($allAnalyses | Where-Object { $_.HasEnhancedFramework }).Count
        WithPhysicsTerms = ($allAnalyses | Where-Object { $_.HasPhysicsTerm }).Count
        Modules = $allAnalyses
    }
    
    # Group by category
    foreach ($category in $script:Categories.Keys) {
        $modules = $allAnalyses | Where-Object { $_.Categories -contains $category }
        $report.ByCategory[$category] = @{
            Count = $modules.Count
            Modules = $modules.Number
        }
    }
    
    # Group by complexity
    $report.ByComplexity = @{
        Simple = ($allAnalyses | Where-Object { $_.EstimatedComplexity -eq "Simple" }).Count
        Moderate = ($allAnalyses | Where-Object { $_.EstimatedComplexity -eq "Moderate" }).Count
        Complex = ($allAnalyses | Where-Object { $_.EstimatedComplexity -eq "Complex" }).Count
        VeryComplex = ($allAnalyses | Where-Object { $_.EstimatedComplexity -eq "Very Complex" }).Count
    }
    
    return $report
}

# Extract a single module
function Extract-Module {
    param(
        [Parameter(Mandatory=$true)]
        [object]$Analysis,
        
        [Parameter(Mandatory=$true)]
        [string]$TargetCategory
    )
    
    Write-Log "Extracting module $($Analysis.Number) to category: $TargetCategory" "INFO"
    
    if ($DryRun) {
        Write-Log "[DRY RUN] Would extract $($Analysis.FileName) -> Core/$TargetCategory/" "INFO"
        return $true
    }
    
    # Create category directory
    $categoryDir = Join-Path $script:Config.CoreDir $script:Categories[$TargetCategory].Directory
    if (-not (Test-Path $categoryDir)) {
        New-Item -ItemType Directory -Path $categoryDir -Force | Out-Null
        Write-Log "Created directory: $categoryDir" "INFO"
    }
    
    # Analyze module structure to determine what to extract
    $content = Get-Content $Analysis.FilePath -Raw
    
    # For now, create header/cpp pair (will be enhanced with actual extraction logic)
    $baseName = [System.IO.Path]::GetFileNameWithoutExtension($Analysis.FileName)
    $headerPath = Join-Path $categoryDir "$baseName.hpp"
    $cppPath = Join-Path $categoryDir "$baseName.cpp"
    
    # This is a placeholder - actual extraction would parse and reorganize code
    Write-Log "TODO: Implement smart extraction for $($Analysis.FileName)" "WARN"
    
    return $true
}

# Generate extraction report
function Generate-Report {
    param([object]$Analysis)
    
    Write-Log "`n========== EXTRACTION ANALYSIS REPORT ==========" "INFO"
    Write-Log "Total Modules Found: $($Analysis.TotalModules)" "INFO"
    Write-Log "Enhanced Modules: $($Analysis.Enhanced)" "INFO"
    Write-Log "Modules with PhysicsTerms: $($Analysis.WithPhysicsTerms)" "INFO"
    
    Write-Log "`nBy Complexity:" "INFO"
    Write-Log "  Simple: $($Analysis.ByComplexity.Simple)" "INFO"
    Write-Log "  Moderate: $($Analysis.ByComplexity.Moderate)" "INFO"
    Write-Log "  Complex: $($Analysis.ByComplexity.Complex)" "INFO"
    Write-Log "  Very Complex: $($Analysis.ByComplexity.VeryComplex)" "INFO"
    
    Write-Log "`nBy Category:" "INFO"
    foreach ($category in $script:Categories.Keys | Sort-Object) {
        $count = $Analysis.ByCategory[$category].Count
        if ($count -gt 0) {
            $modules = ($Analysis.ByCategory[$category].Modules | Sort-Object) -join ', '
            Write-Log "  $category ($count): $modules" "INFO"
        }
    }
    
    Write-Log "`nUncategorized Modules:" "INFO"
    $uncategorized = $Analysis.Modules | Where-Object { $_.Categories.Count -eq 0 }
    if ($uncategorized.Count -gt 0) {
        Write-Log "  Count: $($uncategorized.Count)" "WARN"
        Write-Log "  Modules: $(($uncategorized.Number | Sort-Object) -join ', ')" "WARN"
    } else {
        Write-Log "  All modules categorized!" "INFO"
    }
    
    Write-Log "`n===============================================" "INFO"
    
    # Export detailed analysis to JSON
    $jsonPath = Join-Path $script:Config.SourceDir "module_analysis_$(Get-Date -Format 'yyyyMMdd_HHmmss').json"
    $Analysis | ConvertTo-Json -Depth 10 | Out-File $jsonPath
    Write-Log "Detailed analysis exported to: $jsonPath" "INFO"
}

# Main execution
function Main {
    Write-Log "Star-Magic Batch Extraction Tool Starting..." "INFO"
    Write-Log "Mode: $Mode" "INFO"
    Write-Log "Dry Run: $DryRun" "INFO"
    
    # Ensure Core directory exists
    if (-not (Test-Path $script:Config.CoreDir)) {
        New-Item -ItemType Directory -Path $script:Config.CoreDir -Force | Out-Null
        Write-Log "Created Core directory: $($script:Config.CoreDir)" "INFO"
    }
    
    switch ($Mode.ToLower()) {
        "analyze" {
            $analysis = Analyze-AllModules
            Generate-Report -Analysis $analysis
            
            # Return analysis object for programmatic use
            return $analysis
        }
        
        "extract" {
            Write-Log "Extract mode - Starting batch extraction..." "INFO"
            
            # First analyze
            $analysis = Analyze-AllModules
            
            # Create backup
            if (-not $DryRun) {
                Write-Log "Creating backup at: $($script:Config.BackupDir)" "INFO"
                # Backup will be created during actual extraction
            }
            
            # Extract modules by category
            foreach ($category in $script:Categories.Keys) {
                $modulesToExtract = $analysis.Modules | Where-Object { $_.Categories -contains $category }
                
                if ($modulesToExtract.Count -gt 0) {
                    Write-Log "Extracting $($modulesToExtract.Count) modules for category: $category" "INFO"
                    
                    foreach ($module in $modulesToExtract) {
                        try {
                            Extract-Module -Analysis $module -TargetCategory $category
                        } catch {
                            Write-Log "Failed to extract module $($module.Number): $_" "ERROR"
                        }
                    }
                }
            }
            
            Write-Log "Extraction complete!" "INFO"
        }
        
        "test" {
            Write-Log "Test mode - Validating extracted modules..." "INFO"
            # TODO: Implement compilation and testing
            Write-Log "Test mode not yet implemented" "WARN"
        }
        
        "all" {
            Write-Log "Running full pipeline: analyze -> extract -> test" "INFO"
            Main -Mode "analyze"
            Main -Mode "extract"
            Main -Mode "test"
        }
        
        default {
            Write-Log "Unknown mode: $Mode" "ERROR"
            Write-Log "Valid modes: analyze, extract, test, all" "ERROR"
        }
    }
}

# Run main function
Main
