# ============================================================================
# WOLFRAM EXTRACTION PHASE 1: SOURCE FILE INVENTORY & ANALYSIS
# ============================================================================
# Purpose: Scan all source*.cpp files, extract physics entities that can be 
#          mapped to Wolfram Knowledgebase (constants, particles, quantities)
# Author: CoAnQi UQFF Framework
# Date: November 24, 2025
# ============================================================================

param(
    [switch]$DryRun = $false,
    [int]$BatchSize = 10
)

$workspaceRoot = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
$outputDir = Join-Path $workspaceRoot "wolfram_extraction"

Write-Host "=== WOLFRAM EXTRACTION PHASE 1 ===" -ForegroundColor Cyan
Write-Host "Workspace: $workspaceRoot"
Write-Host "Output: $outputDir"
Write-Host ""

# Create output directory
if (-not (Test-Path $outputDir)) {
    New-Item -ItemType Directory -Path $outputDir | Out-Null
    Write-Host "Created output directory: $outputDir" -ForegroundColor Green
}

# Step 1: Inventory all source files
Write-Host "[Step 1] Inventorying source*.cpp files..." -ForegroundColor Yellow
$sourceFiles = Get-ChildItem -Path $workspaceRoot -Filter "source*.cpp" | 
    Where-Object { $_.Name -match '^source\d+\.cpp$' } | 
    Sort-Object { [int]($_.Name -replace '[^\d]', '') }

Write-Host "Found $($sourceFiles.Count) source files (source1.cpp to source173.cpp)"
Write-Host ""

# Step 2: Analyze file sizes and complexity
Write-Host "[Step 2] Analyzing file characteristics..." -ForegroundColor Yellow
$analysisResults = @()

foreach ($file in $sourceFiles) {
    $content = Get-Content $file.FullName -Raw
    $lineCount = ($content -split "`n").Count
    $sizeKB = [math]::Round($file.Length / 1KB, 2)
    
    # Count potential Wolfram-extractable patterns
    $constants = ([regex]::Matches($content, 'const\s+double\s+\w+')).Count
    $particles = ([regex]::Matches($content, '(Proton|Neutron|Electron|Muon|Quark|Meson|Baryon|Lepton)')).Count
    $frequencies = ([regex]::Matches($content, '(frequency|wavelength|energy)\s*=')).Count
    $masses = ([regex]::Matches($content, 'mass\s*=')).Count
    
    $totalExtractable = $constants + $particles + $frequencies + $masses
    
    $analysisResults += [PSCustomObject]@{
        FileName = $file.Name
        FileNumber = [int]($file.Name -replace '[^\d]', '')
        SizeKB = $sizeKB
        Lines = $lineCount
        Constants = $constants
        Particles = $particles
        Frequencies = $frequencies
        Masses = $masses
        TotalExtractable = $totalExtractable
        Priority = if ($totalExtractable -gt 50) { "HIGH" } elseif ($totalExtractable -gt 10) { "MEDIUM" } else { "LOW" }
    }
}

# Step 3: Generate summary report
Write-Host "[Step 3] Generating analysis report..." -ForegroundColor Yellow

$reportPath = Join-Path $outputDir "phase1_analysis_report.csv"
$analysisResults | Export-Csv -Path $reportPath -NoTypeInformation

Write-Host ""
Write-Host "=== ANALYSIS SUMMARY ===" -ForegroundColor Cyan
Write-Host "Total Files: $($analysisResults.Count)"
Write-Host "HIGH Priority (>50 entities): $($analysisResults | Where-Object Priority -eq 'HIGH' | Measure-Object | Select-Object -ExpandProperty Count)"
Write-Host "MEDIUM Priority (10-50): $($analysisResults | Where-Object Priority -eq 'MEDIUM' | Measure-Object | Select-Object -ExpandProperty Count)"
Write-Host "LOW Priority (<10): $($analysisResults | Where-Object Priority -eq 'LOW' | Measure-Object | Select-Object -ExpandProperty Count)"
Write-Host ""
Write-Host "Total Extractable Entities: $(($analysisResults | Measure-Object -Property TotalExtractable -Sum).Sum)"
Write-Host ""

# Show top 20 files by extraction potential
Write-Host "=== TOP 20 FILES BY EXTRACTION POTENTIAL ===" -ForegroundColor Cyan
$analysisResults | Sort-Object TotalExtractable -Descending | Select-Object -First 20 | 
    Format-Table FileName, SizeKB, Lines, TotalExtractable, Priority -AutoSize

Write-Host ""
Write-Host "Report saved to: $reportPath" -ForegroundColor Green
Write-Host ""

# Step 4: Identify special Wolfram-related files (174-178)
Write-Host "[Step 4] Checking Wolfram-specific source files..." -ForegroundColor Yellow
$wolframFiles = @(174, 175, 176, 177, 178)
$wolframAnalysis = $analysisResults | Where-Object { $_.FileNumber -in $wolframFiles }

if ($wolframAnalysis) {
    Write-Host "Wolfram-specific files (174-178) found:" -ForegroundColor Green
    $wolframAnalysis | Format-Table FileName, TotalExtractable, Priority -AutoSize
} else {
    Write-Host "WARNING: Wolfram-specific files (174-178) not all present!" -ForegroundColor Red
}

Write-Host ""
Write-Host "=== PHASE 1 COMPLETE ===" -ForegroundColor Green
Write-Host "Next: Run Phase 2 to extract entities from high-priority files"
Write-Host ""
