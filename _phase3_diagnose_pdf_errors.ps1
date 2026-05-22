# Phase 3 Diagnostic: Detailed PDF Generation Error Analysis
# Purpose: Identify root cause of 882/1216 PDF generation failures
# Strategy: Test subset (PAPER_001-050) with full error logging enabled

param(
    [int]$StartPaper = 1,
    [int]$EndPaper = 50,
    [switch]$CaptureAllLogs = $false
)

# Setup paths
$pandoc = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$pdflatex = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$markdown_dir = "whitepapers"
$pdf_dir = "pdf"
$log_dir = "pdf_diagnostic_logs"

# Create log directory
if (-not (Test-Path $log_dir)) { New-Item -ItemType Directory -Path $log_dir | Out-Null }

Write-Host "====================================================" -ForegroundColor Cyan
Write-Host "[Phase 3 Diagnostic] PDF Generation Error Analysis" -ForegroundColor Cyan
Write-Host "====================================================" -ForegroundColor Cyan
Write-Host ""
Write-Host "Testing papers: PAPER_$($StartPaper.ToString('000'))-PAPER_$($EndPaper.ToString('000'))"
Write-Host "Log output: $log_dir\"
Write-Host "Total papers to test: $($EndPaper - $StartPaper + 1)"
Write-Host ""

# Verify prerequisites
if (-not (Test-Path $pandoc)) {
    Write-Host "ERROR: pandoc not found at $pandoc" -ForegroundColor Red
    exit 1
}
if (-not (Test-Path $pdflatex)) {
    Write-Host "ERROR: pdflatex not found at $pdflatex" -ForegroundColor Red
    exit 1
}

# Find markdown files matching pattern
$papers = @()
$i = $StartPaper
while ($i -le $EndPaper) {
    $num_str = $i.ToString('000')
    $pattern = "PAPER_${num_str}_*"
    $matching = @(Get-ChildItem -Path $markdown_dir -Filter "$pattern.md" -ErrorAction SilentlyContinue)
    
    if ($matching.Count -gt 0) {
        $papers += $matching[0].FullName
    }
    $i++
}

Write-Host "Found: $($papers.Count) markdown files to convert`n"

$success = 0
$failed = 0
$errors = @()

# Process each paper with detailed logging
foreach ($src in $papers) {
    $filename = Split-Path -Leaf $src
    $basename = [System.IO.Path]::GetFileNameWithoutExtension($filename)
    $pdf = Join-Path $pdf_dir "$basename.pdf"
    $log_file = Join-Path $log_dir "$basename.log"
    $yaml_log = Join-Path $log_dir "$basename.yaml.log"
    $pdflatex_log = Join-Path $log_dir "$basename.pdflatex.log"
    
    Write-Host -NoNewline "Processing: $filename ... "
    
    # Capture ALL output (stdout + stderr)
    $full_output = @()
    $pandoc_error = $null
    
    try {
        # Run pandoc with environment variable to capture pdflatex logs
        $env:TEXMFOUTPUT = (Resolve-Path $log_dir).Path
        
        # Execute with detailed error capture
        $proc = Start-Process -FilePath $pandoc `
            -ArgumentList $src, "-o", $pdf, `
                "--pdf-engine=$pdflatex", `
                '-V', '"geometry:margin=1in"', `
                '-V', '"fontsize=11pt"', `
                '-V', '"colorlinks=true"', `
                '--syntax-highlighting=none' `
            -NoNewWindow -RedirectStandardOutput (Join-Path $log_dir "$basename.stdout.log") `
            -RedirectStandardError (Join-Path $log_dir "$basename.stderr.log") `
            -PassThru -Wait
        
        # Check if PDF was created
        if (Test-Path $pdf) {
            $size = (Get-Item $pdf).Length / 1KB
            Write-Host "✓ SUCCESS ($([math]::Round($size, 1)) KB)" -ForegroundColor Green
            $success++
        } else {
            # PDF not created - read error logs
            $stderr = Get-Content (Join-Path $log_dir "$basename.stderr.log") -ErrorAction SilentlyContinue
            $stdout = Get-Content (Join-Path $log_dir "$basename.stdout.log") -ErrorAction SilentlyContinue
            
            Write-Host "✗ FAILED" -ForegroundColor Red
            $failed++
            
            $error_msg = ""
            if ($stderr) { $error_msg += "STDERR: $([string]::Join(' | ', $stderr[0..2]))" }
            if ($stdout) { $error_msg += " | STDOUT: $([string]::Join(' | ', $stdout[0..2]))" }
            
            $errors += @{
                Paper = $basename
                Error = $error_msg
                StderrFile = Join-Path $log_dir "$basename.stderr.log"
                StdoutFile = Join-Path $log_dir "$basename.stdout.log"
            }
        }
        
    } catch {
        Write-Host "✗ EXCEPTION" -ForegroundColor Red
        $failed++
        $errors += @{
            Paper = $basename
            Error = $_.Exception.Message
            StderrFile = $null
            StdoutFile = $null
        }
    }
}

# Print diagnostic summary
Write-Host ""
Write-Host "====================================================" -ForegroundColor Cyan
Write-Host "[Phase 3 Diagnostic] Results" -ForegroundColor Cyan
Write-Host "====================================================" -ForegroundColor Cyan
Write-Host "Total tested: $($papers.Count)"
Write-Host "Success: $success ($([math]::Round($success/$($papers.Count)*100, 1))%)" -ForegroundColor Green
Write-Host "Failed: $failed ($([math]::Round($failed/$($papers.Count)*100, 1))%)" -ForegroundColor Red
Write-Host ""

# Group errors by pattern
if ($errors.Count -gt 0) {
    Write-Host "Error Patterns:" -ForegroundColor Yellow
    Write-Host "=" * 80
    
    foreach ($error in $errors | Select-Object -First 10) {
        Write-Host ""
        Write-Host "Paper: $($error.Paper)" -ForegroundColor Cyan
        Write-Host "Error: $($error.Error)" -ForegroundColor Red
        
        # Try to extract meaningful error from stderr
        if ($error.StderrFile -and (Test-Path $error.StderrFile)) {
            $stderr_content = Get-Content $error.StderrFile | Select-Object -First 10
            Write-Host "Details:" -ForegroundColor Yellow
            foreach ($line in $stderr_content) {
                Write-Host "  $line"
            }
        }
    }
    
    if ($errors.Count -gt 10) {
        Write-Host ""
        Write-Host "... and $($errors.Count - 10) more failures" -ForegroundColor Red
    }
}

Write-Host ""
Write-Host "Log files saved to: $(Resolve-Path $log_dir)" -ForegroundColor Yellow
Write-Host "Next steps:"
Write-Host "  1. Review error patterns in $log_dir"
Write-Host "  2. Check PAPER_001.stderr.log through PAPER_050.stderr.log for common issues"
Write-Host "  3. Look for patterns: YAML parse errors? LaTeX compilation? Missing fonts? Memory?"
Write-Host "  4. If root cause identified, run full Phase 3 with fix applied"
Write-Host ""
