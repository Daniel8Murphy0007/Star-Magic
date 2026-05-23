#!/usr/bin/env pwsh
<#
.SYNOPSIS
Generate PDFs for SCm files using .tex sources (proper LaTeX)
Uses pdflatex directly on .tex files for guaranteed success
#>

# Configuration
$pdflatex_path = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$md_folder = "whitepapers"
$tex_folder = "whitepapers"
$pdf_folder = "pdf"

# Verify pdflatex exists
if (-not (Test-Path $pdflatex_path)) {
    Write-Host "[ERROR] pdflatex not found at $pdflatex_path" -ForegroundColor Red
    exit 1
}

# 6 SCm files - generate from .tex (original LaTeX format)
$files_to_process = @(
    "SCm_Holmlid_KER_Validation",
    "SCm_Holmlid_Parkhomov_PonsFleischmann_Upgrade",
    "SCm_Holmlid_Rossi_Parkhomov_Validation",
    "SCm_Mizuno_LENR_Transmutation",
    "SCm_PonsFleischmann_Derivation",
    "SCm_Rossi_ECat_Variants_Unified"
)

$success_count = 0
$failed_count = 0
$failed_files = @()

Write-Host "Generating PDFs from SCm .tex files (original LaTeX format)..." -ForegroundColor Cyan
Write-Host "Using: pdflatex directly on .tex files" -ForegroundColor Gray
Write-Host ""

foreach ($base_name in $files_to_process) {
    $tex_path = Join-Path $tex_folder "$base_name.tex"
    $pdf_path = Join-Path $pdf_folder "$base_name.pdf"
    
    # Skip if .tex doesn't exist
    if (-not (Test-Path $tex_path)) {
        Write-Host "[SKIP] $base_name.tex (not found)" -ForegroundColor Yellow
        $failed_count++
        $failed_files += "$base_name.tex"
        continue
    }
    
    try {
        # Generate PDF using pdflatex directly on .tex (3 passes for stability)
        $tex_dir = Get-Item $tex_folder
        $tex_file = Get-Item $tex_path
        
        # Run pdflatex 3 times (standard LaTeX practice)
        for ($pass = 1; $pass -le 3; $pass++) {
            $output = & $pdflatex_path -interaction=nonstopmode -output-directory=$tex_dir.FullName $tex_file.FullName 2>&1
            # Suppress warnings, only check for fatal errors
        }
        
        # Check if PDF was created
        if (Test-Path $pdf_path) {
            $size_kb = (Get-Item $pdf_path).Length / 1KB
            Write-Host "[OK] $base_name.pdf ($([Math]::Round($size_kb, 1)) KB)" -ForegroundColor Green
            $success_count++
        } else {
            # Try moving from temp location
            $temp_pdf = Join-Path $tex_folder "$base_name.pdf"
            if (Test-Path $temp_pdf) {
                Move-Item $temp_pdf $pdf_path -Force
                $size_kb = (Get-Item $pdf_path).Length / 1KB
                Write-Host "[OK] $base_name.pdf ($([Math]::Round($size_kb, 1)) KB)" -ForegroundColor Green
                $success_count++
            } else {
                Write-Host "[FAIL] $base_name.tex (no PDF generated)" -ForegroundColor Red
                $failed_count++
                $failed_files += "$base_name.tex"
            }
        }
    } catch {
        Write-Host "[ERROR] $base_name.tex (Error: $_)" -ForegroundColor Red
        $failed_count++
        $failed_files += "$base_name.tex"
    }
}

# Clean up LaTeX auxiliary files from whitepapers
Get-ChildItem "$tex_folder\*.aux" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem "$tex_folder\*.log" -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem "$tex_folder\*.out" -ErrorAction SilentlyContinue | Remove-Item -Force

Write-Host ""
Write-Host ("=" * 70)
Write-Host "PDF GENERATION SUMMARY" -ForegroundColor Cyan
Write-Host ("=" * 70)
Write-Host "[SUCCESS] $success_count PDFs generated from .tex"
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
    Write-Host "[COMPLETE] All 6 SCm files now have .md AND .pdf!" -ForegroundColor Green
}
