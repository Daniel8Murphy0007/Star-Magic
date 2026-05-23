#!/usr/bin/env pwsh
<#
.SYNOPSIS
Generate PDFs for SCm .md files using EXACT SAME pipeline as PAPER files
Uses pandoc + pdflatex with identical settings
#>

# Configuration (IDENTICAL to PAPER pipeline)
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

# 6 SCm files - EXACT SAME FORMAT as PAPER files
$files_to_process = @(
    "SCm_Holmlid_KER_Validation.md",
    "SCm_Holmlid_Parkhomov_PonsFleischmann_Upgrade.md",
    "SCm_Holmlid_Rossi_Parkhomov_Validation.md",
    "SCm_Mizuno_LENR_Transmutation.md",
    "SCm_PonsFleischmann_Derivation.md",
    "SCm_Rossi_ECat_Variants_Unified.md"
)

$success_count = 0
$failed_count = 0
$failed_files = @()

Write-Host "Generating PDFs for 6 SCm .md files..." -ForegroundColor Cyan
Write-Host "Using EXACT SAME pipeline as 30 PAPER files:" -ForegroundColor Gray
Write-Host "  pandoc + pdflatex with standard settings" -ForegroundColor Gray
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
        # Generate PDF using EXACT SAME pandoc + pdflatex command as PAPER files
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
Write-Host "[SUCCESS] $success_count SCm PDFs generated"
Write-Host "[FAILED]  $failed_count files"

if ($failed_files.Count -gt 0) {
    Write-Host ""
    Write-Host "Failed files:" -ForegroundColor Yellow
    foreach ($f in $failed_files) {
        Write-Host "  - $f"
    }
}

Write-Host ""
Write-Host "[OUTPUT] PDFs now in: $((Get-Item $pdf_folder).FullName)"
Write-Host "[COUNT] Total PDFs in pdf/ folder: $((Get-ChildItem "$pdf_folder\*.pdf" -ErrorAction SilentlyContinue | Measure-Object).Count)"

if ($success_count -eq 6) {
    Write-Host ""
    Write-Host "[COMPLETE] All 6 SCm files: .md ✓ .pdf ✓ UNIFORM!" -ForegroundColor Green
}
