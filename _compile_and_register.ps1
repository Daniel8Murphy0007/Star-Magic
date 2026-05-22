<#
_compile_and_register.ps1 — Post-edit pipeline for v5.78 paper authoring.

Steps:
  1. Two-pass pdflatex (direct, NO PANDOC) on whitepapers/<paper>.tex
  2. Move PDF to pdf/<paper>.pdf
  3. Harvest `% CLOSURE ::` banners → master_closures.csv (idempotent)
  4. Re-run _profile_master_ledger.py → ledger CSVs
  5. git add (NO auto-commit — review then commit manually)

Usage:
    .\_compile_and_register.ps1 PAPER_0500_lambda_qcd
    .\_compile_and_register.ps1 PAPER_0500_lambda_qcd -Force   # overwrite existing closures
    .\_compile_and_register.ps1 PAPER_0500_lambda_qcd -SkipProfile

Path layout:
    Source : whitepapers\<paper>.tex
    Output : pdf\<paper>.pdf
#>
[CmdletBinding()]
param(
    [Parameter(Mandatory=$true, Position=0)]
    [string]$PaperSlug,
    [switch]$Force,
    [switch]$SkipProfile,
    [switch]$NoGitAdd
)

$ErrorActionPreference = 'Stop'
$RepoRoot   = Split-Path -Parent $MyInvocation.MyCommand.Path
$PdfLatex   = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$Python     = Join-Path $RepoRoot ".venv_py314_backup\Scripts\python.exe"
$WhitepDir  = Join-Path $RepoRoot "whitepapers"
$PdfDir     = Join-Path $RepoRoot "pdf"

if (-not (Test-Path $PdfLatex)) { throw "pdflatex not found: $PdfLatex" }
if (-not (Test-Path $Python))   { throw "python venv not found: $Python" }

# Strip .tex if user supplied it.
$slug = $PaperSlug -replace '\.tex$',''
$texPath = Join-Path $WhitepDir "$slug.tex"
if (-not (Test-Path $texPath)) { throw "Source not found: $texPath" }

Write-Host "==> 1/4 pdflatex (two-pass): $slug.tex" -ForegroundColor Cyan
Push-Location $WhitepDir
try {
    for ($pass = 1; $pass -le 2; $pass++) {
        & $PdfLatex -interaction=nonstopmode -halt-on-error "$slug.tex" | Out-Null
        if ($LASTEXITCODE -ne 0) {
            Pop-Location
            throw "pdflatex pass $pass failed (exit=$LASTEXITCODE). Check $slug.log"
        }
    }
    # Move PDF, clean aux/log/out.
    $srcPdf = Join-Path $WhitepDir "$slug.pdf"
    $dstPdf = Join-Path $PdfDir    "$slug.pdf"
    if (-not (Test-Path $PdfDir)) { New-Item -ItemType Directory -Path $PdfDir | Out-Null }
    Move-Item -Force $srcPdf $dstPdf
    foreach ($ext in '.aux','.log','.out','.toc') {
        $f = Join-Path $WhitepDir "$slug$ext"
        if (Test-Path $f) { Remove-Item $f -Force }
    }
    Write-Host "    PDF -> pdf\$slug.pdf" -ForegroundColor Green
} finally {
    Pop-Location
}

Write-Host "==> 2/4 harvest banners -> master_closures.csv" -ForegroundColor Cyan
$harvestArgs = @("$slug.tex")
if ($Force) { $harvestArgs += "--force" }
& $Python (Join-Path $RepoRoot "_harvest_template_banners.py") @harvestArgs
if ($LASTEXITCODE -ne 0) { throw "harvest failed (exit=$LASTEXITCODE)" }

if (-not $SkipProfile) {
    Write-Host "==> 3/4 re-profile ledger" -ForegroundColor Cyan
    & $Python (Join-Path $RepoRoot "_profile_master_ledger.py")
    if ($LASTEXITCODE -ne 0) { throw "profiler failed (exit=$LASTEXITCODE)" }
} else {
    Write-Host "==> 3/4 SKIPPED (--SkipProfile)" -ForegroundColor Yellow
}

if (-not $NoGitAdd) {
    Write-Host "==> 4/4 git add (no commit)" -ForegroundColor Cyan
    Push-Location $RepoRoot
    try {
        git add "whitepapers/$slug.tex" "pdf/$slug.pdf" master_closures.csv `
                MASTER_LEDGER_BY_CATEGORY.csv MASTER_LEDGER_BY_STATUS.csv `
                MASTER_LEDGER_BY_SCRIPT.csv LEDGER_VS_PRIMITIVES_XREF.csv 2>$null
        git status --short
    } finally { Pop-Location }
} else {
    Write-Host "==> 4/4 SKIPPED (--NoGitAdd)" -ForegroundColor Yellow
}

Write-Host "`nDone. Review with 'git diff --cached' then commit manually." -ForegroundColor Green
