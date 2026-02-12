# Create arXiv submission package for UQFF Production Manuscript
# PowerShell version for Windows
# Usage: .\create_arxiv_package.ps1

$ErrorActionPreference = "Stop"

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "arXiv Submission Package Creator" -ForegroundColor Cyan
Write-Host "UQFF Production Manuscript" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Configuration
$PackageName = "uqff_arxiv_submission"
$Manuscript = "uqff_production_arxiv.tex"
$FiguresDir = "figures"
$AncDir = "anc"

# Clean previous build
Write-Host "[1/6] Cleaning previous builds..." -ForegroundColor Yellow
if (Test-Path $PackageName) {
    Remove-Item -Recurse -Force $PackageName
}
if (Test-Path "$PackageName.zip") {
    Remove-Item -Force "$PackageName.zip"
}
Get-ChildItem -Path manuscript -Filter "uqff_production_arxiv.*" | 
    Where-Object { $_.Extension -match '^\.(aux|log|out)$' } | 
    Remove-Item -Force -ErrorAction SilentlyContinue
Write-Host "✓ Cleaned" -ForegroundColor Green

# Create package directory structure
Write-Host "[2/6] Creating package structure..." -ForegroundColor Yellow
New-Item -ItemType Directory -Path "$PackageName" -Force | Out-Null
New-Item -ItemType Directory -Path "$PackageName\$FiguresDir" -Force | Out-Null
New-Item -ItemType Directory -Path "$PackageName\$AncDir" -Force | Out-Null
Write-Host "✓ Directories created" -ForegroundColor Green

# Generate figures
Write-Host "[3/6] Generating figures..." -ForegroundColor Yellow
Push-Location manuscript
try {
    python generate_figures.py
    if ($LASTEXITCODE -ne 0) {
        throw "Figure generation failed"
    }
} catch {
    Write-Host "✗ Figure generation failed: $_" -ForegroundColor Red
    Pop-Location
    exit 1
}
Pop-Location
Write-Host "✓ Figures generated" -ForegroundColor Green

# Compile manuscript (twice for references)
Write-Host "[4/6] Compiling manuscript..." -ForegroundColor Yellow
Push-Location manuscript
try {
    pdflatex -interaction=nonstopmode $Manuscript | Out-Null
    pdflatex -interaction=nonstopmode $Manuscript | Out-Null
    if (-not (Test-Path "uqff_production_arxiv.pdf")) {
        throw "PDF not generated"
    }
} catch {
    Write-Host "✗ Compilation failed: $_" -ForegroundColor Red
    Pop-Location
    exit 1
}
Pop-Location
Write-Host "✓ Manuscript compiled (31 pages)" -ForegroundColor Green

# Copy files to package
Write-Host "[5/6] Copying files to package..." -ForegroundColor Yellow
Copy-Item "manuscript\$Manuscript" "$PackageName\"
Copy-Item "manuscript\$FiguresDir\*.pdf" "$PackageName\$FiguresDir\"
Copy-Item "manuscript\generate_figures.py" "$PackageName\$AncDir\"
Copy-Item "manuscript\README_ARXIV.md" "$PackageName\$AncDir\"
Copy-Item "manuscript\SUBMISSION_CHECKLIST.md" "$PackageName\$AncDir\"
Write-Host "✓ Files copied" -ForegroundColor Green

# Create ZIP archive
Write-Host "[6/6] Creating submission archive..." -ForegroundColor Yellow
Compress-Archive -Path $PackageName -DestinationPath "$PackageName.zip" -Force
$Size = (Get-Item "$PackageName.zip").Length / 1MB
$SizeStr = "{0:N2} MB" -f $Size
Write-Host "✓ Archive created: $PackageName.zip ($SizeStr)" -ForegroundColor Green

# Summary
Write-Host ""
Write-Host "========================================" -ForegroundColor Cyan
Write-Host "✓ PACKAGE READY FOR SUBMISSION" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Contents:" -ForegroundColor White
Write-Host "  - Main manuscript: $Manuscript"
Write-Host "  - Figures: 5 PDF files"
Write-Host "  - Ancillary: 3 supplementary files"
Write-Host ""
Write-Host "Archive: $PackageName.zip ($SizeStr)" -ForegroundColor Yellow
Write-Host ""
Write-Host "Next steps:" -ForegroundColor White
Write-Host "  1. Review: manuscript\uqff_production_arxiv.pdf"
Write-Host "  2. Upload: $PackageName.zip to arxiv.org/submit"
Write-Host "  3. Categories: gr-qc (primary), astro-ph.CO, hep-ph, physics.comp-ph"
Write-Host ""
Write-Host "Submission URL: https://arxiv.org/submit" -ForegroundColor Cyan
Write-Host ""

# Optional: Open PDF for review
$Response = Read-Host "Open PDF for review? [y/N]"
if ($Response -match '^[Yy]$') {
    Start-Process "manuscript\uqff_production_arxiv.pdf"
}

Write-Host "Done!" -ForegroundColor Green
