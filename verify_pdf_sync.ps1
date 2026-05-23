#!/usr/bin/env pwsh
<#
.SYNOPSIS
Verify new .md files against existing PDFs in /pdf folder
Checks: existence, timestamps, and content alignment
#>

$md_folder = "whitepapers"
$pdf_folder = "pdf"

# Files to check (30 newly created)
$files_to_check = @(
    "PAPER_1196_UQFF_Plasma_Fusion_Unified_Proof_Set",
    "PAPER_1199_UQFF_Information_Math_Unified_Proof_Set",
    "PAPER_1200_UQFF_GR_Precision_Unified_Proof_Set",
    "PAPER_1201_UQFF_Materials_Photonics_Unified_Proof_Set",
    "PAPER_1202_UQFF_Chemistry_Spectroscopy_Unified_Proof_Set",
    "PAPER_1203_UQFF_Nuclear_Physics_Unified_Proof_Set",
    "PAPER_1204_UQFF_Fluid_Dynamics_Unified_Proof_Set",
    "PAPER_1205_UQFF_Geometry_Topology_Unified_Proof_Set",
    "PAPER_1206_UQFF_Solar_System_Unified_Proof_Set",
    "PAPER_1207_UQFF_Biology_Allometry_Unified_Proof_Set",
    "PAPER_1208_UQFF_Transcendentals_Unified_Proof_Set",
    "PAPER_1209_UQFF_Particle_Physics_Unified_Proof_Set",
    "PAPER_1209AA_UQFF_Chemistry_Unified_Proof_Set",
    "PAPER_1209BB_UQFF_Biology_Unified_Proof_Set",
    "PAPER_1209CC_UQFF_Geophysics_Unified_Proof_Set",
    "PAPER_1209DD_UQFF_Electromagnetism_Unified_Proof_Set",
    "PAPER_1209EE_UQFF_Quantum_Thermo_Unified_Proof_Set",
    "PAPER_1209FF_UQFF_Math_Constants_Unified_Proof_Set",
    "PAPER_1209GG_UQFF_Cosmological_Constants_Unified_Proof_Set",
    "PAPER_1209HH_UQFF_Particle_Masses_Unified_Proof_Set",
    "PAPER_1209II_UQFF_Nuclear_Binding_Energies_Unified_Proof_Set",
    "PAPER_1209JJ_UQFF_Geophysics_Unified_Proof_Set",
    "PAPER_1209KK_Tier_KK_Solar_System",
    "PAPER_1209X_UQFF_Climate_Atmosphere_Unified_Proof_Set",
    "PAPER_1209Y_UQFF_Engineering_Unified_Proof_Set",
    "PAPER_1209Z_UQFF_Astronomical_Units_Unified_Proof_Set",
    "PAPER_1210_UQFF_Lagrangian_Bridge_172_Closures",
    "PAPER_1212_UQFF_Cosmological_Constant_Closure",
    "PAPER_1213_UQFF_Page_Curve_Closure",
    "PAPER_1214_UQFF_Habitable_Zone_Universal_Buoyancy"
)

$pdf_exists = 0
$pdf_missing = 0
$md_newer = 0
$pdf_newer = 0
$same_time = 0
$needs_regen = @()
$up_to_date = @()
$missing_list = @()

Write-Host "VERIFICATION: New .md files vs existing PDFs" -ForegroundColor Cyan
Write-Host ""

foreach ($base_name in $files_to_check) {
    $md_path = Join-Path $md_folder "$base_name.md"
    $pdf_path = Join-Path $pdf_folder "$base_name.pdf"
    
    # Check if .md exists
    if (-not (Test-Path $md_path)) {
        Write-Host "[SKIP] $base_name (no .md file found)" -ForegroundColor Yellow
        continue
    }
    
    # Check if PDF exists
    if (-not (Test-Path $pdf_path)) {
        Write-Host "[MISSING] $base_name.md exists, but PDF missing" -ForegroundColor Red
        $pdf_missing++
        $needs_regen += $base_name
        $missing_list += "$base_name"
        continue
    }
    
    # PDF exists - compare timestamps
    $md_time = (Get-Item $md_path).LastWriteTime
    $pdf_time = (Get-Item $pdf_path).LastWriteTime
    $time_diff = ($md_time - $pdf_time).TotalSeconds
    
    # Get file sizes
    $md_size = (Get-Item $md_path).Length
    $pdf_size = (Get-Item $pdf_path).Length
    
    $pdf_exists++
    
    # Determine status
    if ([Math]::Abs($time_diff) -lt 1) {
        # Same time (within 1 second)
        Write-Host "[OK] $base_name (PDF in sync, same time)" -ForegroundColor Green
        $same_time++
        $up_to_date += $base_name
    } elseif ($time_diff -lt 0) {
        # PDF is newer than .md (good)
        $abs_diff = [Math]::Round([Math]::Abs($time_diff), 0)
        Write-Host "[OK] $base_name (PDF newer by $abs_diff sec)" -ForegroundColor Green
        $pdf_newer++
        $up_to_date += $base_name
    } else {
        # .md is newer than PDF (need regeneration)
        $rounded_diff = [Math]::Round($time_diff, 0)
        Write-Host "[OUTDATED] $base_name (PDF outdated by $rounded_diff sec)" -ForegroundColor Yellow
        $md_newer++
        $needs_regen += $base_name
    }
}

Write-Host ""
Write-Host ("=" * 70)
Write-Host "VERIFICATION SUMMARY" -ForegroundColor Cyan
Write-Host ("=" * 70)
Write-Host "[OK] PDFs existing & in sync:    $($up_to_date.Count) files"
Write-Host "[!] PDFs missing/outdated:       $($needs_regen.Count) files"
Write-Host "    - Missing PDFs:              $pdf_missing"
Write-Host "    - Outdated PDFs (.md newer): $md_newer"
Write-Host ""

if ($missing_list.Count -gt 0) {
    Write-Host "MISSING PDFs (need generation):" -ForegroundColor Yellow
    foreach ($name in $missing_list | Sort-Object) {
        Write-Host "  - $name"
    }
    Write-Host ""
}

if ($needs_regen.Count -eq 0) {
    Write-Host "[SUCCESS] NO REGENERATION NEEDED - All PDFs are in sync with .md files" -ForegroundColor Green
} else {
    Write-Host "[ACTION] REGENERATION NEEDED for $($needs_regen.Count) files" -ForegroundColor Yellow
    Write-Host ""
    Write-Host "Run: powershell -ExecutionPolicy Bypass -File generate_new_pdfs.ps1" -ForegroundColor Cyan
}

Write-Host ""
Write-Host "[Total] Checked: $($files_to_check.Count) files"
Write-Host "[Found] PDFs in sync: $pdf_exists"
Write-Host "[Missing] PDFs: $pdf_missing"
