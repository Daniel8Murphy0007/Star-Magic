#!/usr/bin/env pwsh
# _generate_all_pdfs.ps1  — Regenerate ALL whitepapers -> PDFs
# Generated: April 23, 2026

$pandoc  = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$xe      = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
$srcDir  = "whitepapers"
$outDir  = "pdf"

if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

$papers   = @(Get-ChildItem "$srcDir\*.md" | Sort-Object Name)
$total    = $papers.Count
$ok       = 0
$fail     = 0
$skipped  = 0
$failList = @()

Write-Host "=== PDF Generation: $total whitepapers ===" -ForegroundColor Cyan
Write-Host ""

foreach ($src in $papers) {
    $stem = $src.BaseName
    $pdf  = Join-Path $outDir "$stem.pdf"
    $idx  = $ok + $fail + $skipped + 1

    Write-Host -NoNewline "[$idx/$total] $($src.Name) ... "

    $errOut = & $pandoc $src.FullName -o $pdf `
        "--pdf-engine=$xe" `
        -V "geometry:margin=1in" `
        -V "fontsize=11pt" `
        -V "colorlinks=true" `
        "--highlight-style=tango" `
        2>&1

    if ($LASTEXITCODE -eq 0) {
        $sz = (Get-Item $pdf).Length
        Write-Host "OK ($([math]::Round($sz/1024,0)) KB)" -ForegroundColor Green
        $ok++
    } else {
        Write-Host "FAIL" -ForegroundColor Red
        if ($errOut) { Write-Host "  ERR: $errOut" -ForegroundColor DarkRed }
        $fail++
        $failList += $src.Name
    }
}

Write-Host ""
Write-Host "================================================" -ForegroundColor Cyan
Write-Host "  RESULT: $ok OK   $fail FAILED   $skipped SKIPPED" -ForegroundColor Cyan
Write-Host "  Total:  $total" -ForegroundColor Cyan
Write-Host "================================================" -ForegroundColor Cyan

if ($failList.Count -gt 0) {
    Write-Host "`nFailed papers:" -ForegroundColor Yellow
    $failList | ForEach-Object { Write-Host "  - $_" -ForegroundColor Yellow }
    $failList | Out-File "_pdf_failures.txt" -Encoding utf8
    Write-Host "  (saved to _pdf_failures.txt)"
}
