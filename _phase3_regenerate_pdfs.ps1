# Phase 3: Regenerate all 1216+ PDFs using pdflatex
#
# PDF generation workflow:
#  - Source: whitepapers/*.md (markdown files with YAML frontmatter)
#  - Engine: pdflatex (arXiv-approved, NOT xelatex)
#  - Output: pdf/*.pdf (canonical user-facing PDFs)
#  - Pipeline: pandoc .md -> pdflatex two-pass -> pdf/

$pandoc = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$xelatex = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
$srcDir = "whitepapers"
$outDir = "pdf"

# Verify tools exist
if (!(Test-Path $pandoc)) {
    Write-Host "ERROR: pandoc not found at $pandoc"
    exit 1
}
if (!(Test-Path $xelatex)) {
    Write-Host "ERROR: xelatex not found at $xelatex"
    exit 1
}

# Create output directory
if (!(Test-Path $outDir)) {
    New-Item -ItemType Directory -Path $outDir | Out-Null
}

Write-Host "[Phase 3] Starting PDF regeneration"
Write-Host "  pandoc: $pandoc"
Write-Host "  xelatex: $xelatex"
Write-Host ""

# Get all PAPER_*.md files
$papers = @(Get-ChildItem -Path $srcDir -Filter "PAPER_*.md" -File | Sort-Object Name)
$total = $papers.Count

Write-Host "[Phase 3] Found $total papers to convert"
Write-Host ""

$ok = 0
$fail = 0
$failList = @()
$startTime = Get-Date

$i = 0
foreach ($paper in $papers) {
    $stem = $paper.BaseName
    $pdf = Join-Path $outDir "$stem.pdf"
    $src = $paper.FullName
    
    # Progress indicator every 50 papers
    if (($i % 50) -eq 0) {
        Write-Host "[Phase 3] Processing $i/$total..."
    }
    $i++
    
    # Convert MD to PDF: Pandoc with xelatex (Unicode math support)
    & $pandoc $src -o $pdf --pdf-engine=$xelatex 2>&1 | Out-Null
    
    $success = $false
    if (Test-Path $pdf) {
        $success = $true
    }
    
    if ($success) {
        $ok++
    } else {
        Write-Host "  FAIL: $stem"
        $fail++
        $failList += $stem
    }
}

$elapsed = (Get-Date) - $startTime
$rate = [math]::Round($total / $elapsed.TotalSeconds, 1)

Write-Host ""
Write-Host "=========================================="
Write-Host "[Phase 3] PDF Regeneration Complete"
Write-Host "=========================================="
Write-Host "Total papers: $total"
Write-Host "Success: $ok"
Write-Host "Failed: $fail"
Write-Host "Rate: $rate papers/sec"
Write-Host "Elapsed: $($elapsed.TotalMinutes.ToString('F1')) min"
Write-Host ""

if ($fail -gt 0) {
    Write-Host "Failed papers ($fail):"
    $failList | ForEach-Object { Write-Host "  $_" }
    Write-Host ""
    Write-Host "Next steps:"
    Write-Host "  1. Check pandoc/pdflatex installation"
    Write-Host "  2. Manually convert failed papers or investigate LaTeX errors"
}

Write-Host ""
Write-Host "Verify outputs:"
Write-Host "  cd pdf; (Get-ChildItem *.pdf | Measure-Object).Count"
