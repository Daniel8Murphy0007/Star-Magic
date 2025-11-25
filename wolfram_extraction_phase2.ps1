# ============================================================================
# WOLFRAM EXTRACTION PHASE 2: ENTITY EXTRACTION ENGINE
# ============================================================================
# Purpose: Extract physics entities from source files and prepare for Wolfram
#          Knowledgebase validation
# ============================================================================

param(
    [string[]]$SourceFiles = @(),  # Specific files to process, or all if empty
    [int]$TopN = 20,               # Process top N files by extraction potential
    [switch]$ProcessAll = $false   # Process all 163 files
)

$workspaceRoot = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
$outputDir = Join-Path $workspaceRoot "wolfram_extraction"
$entitiesFile = Join-Path $outputDir "extracted_entities.json"

Write-Host "=== WOLFRAM EXTRACTION PHASE 2 ===" -ForegroundColor Cyan
Write-Host ""

# Load Phase 1 analysis
$analysisPath = Join-Path $outputDir "phase1_analysis_report.csv"
if (-not (Test-Path $analysisPath)) {
    Write-Host "ERROR: Run Phase 1 first (wolfram_extraction_phase1.ps1)" -ForegroundColor Red
    exit 1
}

$analysis = Import-Csv $analysisPath | Sort-Object { [int]$_.TotalExtractable } -Descending

# Determine files to process
if ($ProcessAll) {
    $filesToProcess = $analysis
    Write-Host "Processing ALL $($filesToProcess.Count) source files..." -ForegroundColor Yellow
} elseif ($SourceFiles.Count -gt 0) {
    $filesToProcess = $analysis | Where-Object { $_.FileName -in $SourceFiles }
    Write-Host "Processing $($filesToProcess.Count) specified files..." -ForegroundColor Yellow
} else {
    $filesToProcess = $analysis | Select-Object -First $TopN
    Write-Host "Processing top $TopN files by extraction potential..." -ForegroundColor Yellow
}

Write-Host ""

# Entity extraction patterns
$patterns = @{
    PhysicalConstant = @(
        'const\s+double\s+(\w+)\s*=\s*([0-9.eE+-]+)'  # const double NAME = VALUE
        '(G|c|h|k_B|m_e|m_p)\s*=\s*([0-9.eE+-]+)'     # Common physics constants
    )
    Particle = @(
        'class\s+(\w*(?:Proton|Neutron|Electron|Muon|Quark|Meson|Baryon|Lepton)\w*)'
        '(Proton|Neutron|Electron|Muon|Tau|Quark|Gluon|Photon|W|Z)(?:Mass|Charge|Spin)'
    )
    Quantity = @(
        '(mass|charge|spin|energy|momentum|frequency|wavelength|temperature)\s+(\w+)'
        'double\s+(\w*(?:Mass|Energy|Freq|Temp|Press|Dens)\w*)'
    )
    AstrophysicalSystem = @(
        '(NGC|M\d+|IC|Abell|SMBH|SGR|Magnetar|Galaxy|Nebula|Cluster)[\w\d-]+'
        'class\s+(\w*(?:System|Galaxy|Star|Nebula|Cluster)\w*)'
    )
}

# Storage for extracted entities
$allEntities = @{
    PhysicalConstants = @()
    Particles = @()
    Quantities = @()
    AstrophysicalSystems = @()
    SourceFileMap = @{}
}

$fileCount = 0
foreach ($fileInfo in $filesToProcess) {
    $fileCount++
    $filePath = Join-Path $workspaceRoot $fileInfo.FileName
    
    Write-Host "[$fileCount/$($filesToProcess.Count)] Processing $($fileInfo.FileName)..." -ForegroundColor Gray
    
    if (-not (Test-Path $filePath)) {
        Write-Host "  WARNING: File not found, skipping" -ForegroundColor Yellow
        continue
    }
    
    $content = Get-Content $filePath -Raw
    $entities = @{
        PhysicalConstants = @()
        Particles = @()
        Quantities = @()
        Systems = @()
    }
    
    # Extract Physical Constants
    foreach ($pattern in $patterns.PhysicalConstant) {
        $matches = [regex]::Matches($content, $pattern)
        foreach ($match in $matches) {
            if ($match.Groups.Count -ge 3) {
                $entities.PhysicalConstants += @{
                    Name = $match.Groups[1].Value
                    Value = $match.Groups[2].Value
                    SourceFile = $fileInfo.FileName
                }
            }
        }
    }
    
    # Extract Particles
    foreach ($pattern in $patterns.Particle) {
        $matches = [regex]::Matches($content, $pattern)
        foreach ($match in $matches) {
            $particleName = $match.Groups[1].Value
            if ($particleName -and $particleName.Length -gt 2) {
                $entities.Particles += @{
                    Name = $particleName
                    SourceFile = $fileInfo.FileName
                }
            }
        }
    }
    
    # Extract Astrophysical Systems
    foreach ($pattern in $patterns.AstrophysicalSystem) {
        $matches = [regex]::Matches($content, $pattern)
        foreach ($match in $matches) {
            $systemName = $match.Groups[1].Value
            if ($systemName -and $systemName.Length -gt 2) {
                $entities.Systems += @{
                    Name = $systemName
                    SourceFile = $fileInfo.FileName
                }
            }
        }
    }
    
    # Store entities
    $allEntities.PhysicalConstants += $entities.PhysicalConstants
    $allEntities.Particles += $entities.Particles
    $allEntities.AstrophysicalSystems += $entities.Systems
    $allEntities.SourceFileMap[$fileInfo.FileName] = $entities
    
    $totalFound = $entities.PhysicalConstants.Count + $entities.Particles.Count + $entities.Systems.Count
    Write-Host "  Found: $totalFound entities ($($entities.PhysicalConstants.Count) constants, $($entities.Particles.Count) particles, $($entities.Systems.Count) systems)"
}

# Save results
Write-Host ""
Write-Host "Saving extraction results..." -ForegroundColor Yellow
$allEntities | ConvertTo-Json -Depth 10 | Set-Content $entitiesFile

Write-Host ""
Write-Host "=== EXTRACTION SUMMARY ===" -ForegroundColor Cyan
Write-Host "Files Processed: $fileCount"
Write-Host "Physical Constants: $($allEntities.PhysicalConstants.Count)"
Write-Host "Particles: $($allEntities.Particles.Count)"
Write-Host "Astrophysical Systems: $($allEntities.AstrophysicalSystems.Count)"
Write-Host ""
Write-Host "Results saved to: $entitiesFile" -ForegroundColor Green
Write-Host ""
Write-Host "Next: Run Phase 3 to validate entities with Wolfram Knowledgebase"
Write-Host ""
