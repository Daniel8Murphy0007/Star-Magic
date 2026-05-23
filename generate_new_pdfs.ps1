#!/usr/bin/env pwsh
<#
.SYNOPSIS
Generate PDFs for newly converted .md files (PAPER_1196, 1199-1214)
Uses pandoc with pdflatex engine (arXiv-approved workflow)
#>

# Configuration
$pandoc_path = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$pdflatex_path = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$md_folder = "whitepapers"
$pdf_folder = "pdf"

# Verify pandoc exists
if (-not (Test-Path $pandoc_path)) {
    Write-Host "[ERROR] Pandoc not found at $pandoc_path" -ForegroundColor Red
    exit 1
}

# Verify pdflatex exists
if (-not (Test-Path $pdflatex_path)) {
    Write-Host "[ERROR] pdflatex not found at $pdflatex_path" -ForegroundColor Red
    exit 1
}

# Create pdf folder if it doesn't exist
if (-not (Test-Path $pdf_folder)) {
    New-Item -ItemType Directory -Path $pdf_folder -Force | Out-Null
    Write-Host "[CREATED] $pdf_folder folder"
}

# Files to process (30 newly created .md files)
$files_to_process = @(
    "PAPER_1196_UQFF_Plasma_Fusion_Unified_Proof_Set.md",
    "PAPER_1199_UQFF_Information_Math_Unified_Proof_Set.md",
    "PAPER_1200_UQFF_GR_Precision_Unified_Proof_Set.md",
    "PAPER_1201_UQFF_Materials_Photonics_Unified_Proof_Set.md",
    "PAPER_1202_UQFF_Chemistry_Spectroscopy_Unified_Proof_Set.md",
    "PAPER_1203_UQFF_Nuclear_Physics_Unified_Proof_Set.md",
    "PAPER_1204_UQFF_Fluid_Dynamics_Unified_Proof_Set.md",
    "PAPER_1205_UQFF_Geometry_Topology_Unified_Proof_Set.md",
    "PAPER_1206_UQFF_Solar_System_Unified_Proof_Set.md",
    "PAPER_1207_UQFF_Biology_Allometry_Unified_Proof_Set.md",
    "PAPER_1208_UQFF_Transcendentals_Unified_Proof_Set.md",
    "PAPER_1209_UQFF_Particle_Physics_Unified_Proof_Set.md",
    "PAPER_1209AA_UQFF_Chemistry_Unified_Proof_Set.md",
    "PAPER_1209BB_UQFF_Biology_Unified_Proof_Set.md",
    "PAPER_1209CC_UQFF_Geophysics_Unified_Proof_Set.md",
    "PAPER_1209DD_UQFF_Electromagnetism_Unified_Proof_Set.md",
    "PAPER_1209EE_UQFF_Quantum_Thermo_Unified_Proof_Set.md",
    "PAPER_1209FF_UQFF_Math_Constants_Unified_Proof_Set.md",
    "PAPER_1209GG_UQFF_Cosmological_Constants_Unified_Proof_Set.md",
    "PAPER_1209HH_UQFF_Particle_Masses_Unified_Proof_Set.md",
    "PAPER_1209II_UQFF_Nuclear_Binding_Energies_Unified_Proof_Set.md",
    "PAPER_1209JJ_UQFF_Geophysics_Unified_Proof_Set.md",
    "PAPER_1209KK_Tier_KK_Solar_System.md",
    "PAPER_1209X_UQFF_Climate_Atmosphere_Unified_Proof_Set.md",
    "PAPER_1209Y_UQFF_Engineering_Unified_Proof_Set.md",
    "PAPER_1209Z_UQFF_Astronomical_Units_Unified_Proof_Set.md",
    "PAPER_1210_UQFF_Lagrangian_Bridge_172_Closures.md",
    "PAPER_1212_UQFF_Cosmological_Constant_Closure.md",
    "PAPER_1213_UQFF_Page_Curve_Closure.md",
    "PAPER_1214_UQFF_Habitable_Zone_Universal_Buoyancy.md"
)

$success_count = 0
$failed_count = 0
$failed_files = @()

Write-Host "Starting PDF generation for 30 files..." -ForegroundColor Cyan
Write-Host "Using: pdflatex (arXiv-approved, UQFF canonical)" -ForegroundColor Gray
Write-Host ""

foreach ($file in $files_to_process) {
    $md_path = Join-Path $md_folder $file
    $pdf_name = $file -replace '\.md$', '.pdf'
    $pdf_path = Join-Path $pdf_folder $pdf_name
    
    # Skip if .md doesn't exist
    if (-not (Test-Path $md_path)) {
        Write-Host "[SKIP] $file (source .md not found)" -ForegroundColor Yellow
        $failed_count++
        $failed_files += $file
        continue
    }
    
    try {
        # Generate PDF using pandoc + pdflatex
        & $pandoc_path $md_path `
            --pdf-engine=$pdflatex_path `
            -V geometry:margin=1in `
            -V fontsize=11pt `
            -V colorlinks=true `
            --highlight-style=tango `
            -o $pdf_path
        
        if (Test-Path $pdf_path) {
            $size_kb = (Get-Item $pdf_path).Length / 1KB
            Write-Host "[OK] $pdf_name ($([Math]::Round($size_kb, 1)) KB)" -ForegroundColor Green
            $success_count++
        } else {
            Write-Host "[FAIL] $file -- $pdf_name (generation failed, no output)" -ForegroundColor Red
            $failed_count++
            $failed_files += $file
        }
    } catch {
        Write-Host "[ERROR] $file -- $pdf_name (Error: $_)" -ForegroundColor Red
        $failed_count++
        $failed_files += $file
    }
}

Write-Host ""
Write-Host ("=" * 70)
Write-Host "PDF GENERATION SUMMARY" -ForegroundColor Cyan
Write-Host ("=" * 70)
Write-Host "[SUCCESS] $success_count PDFs generated"
Write-Host "[FAILED]  $failed_count files"

if ($failed_files.Count -gt 0) {
    Write-Host ""
    Write-Host "Failed files:" -ForegroundColor Yellow
    foreach ($f in $failed_files) {
        Write-Host "  - $f"
    }
}

Write-Host ""
Write-Host "[OUTPUT] Output folder: $((Get-Item $pdf_folder).FullName)"
Write-Host "[COUNT] Total PDFs in pdf/ folder: $((Get-ChildItem "$pdf_folder\*.pdf" -ErrorAction SilentlyContinue | Measure-Object).Count)"
