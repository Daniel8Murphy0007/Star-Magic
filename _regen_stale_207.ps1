$ErrorActionPreference = 'Continue'
$pandoc = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$pdflatex = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$header = "_pdf_unicode_header.tex"
$list = Get-Content "_stale_pdfs_to_regen.txt" | Where-Object { $_.Trim() -ne "" }
$total = $list.Count
$ok = 0
$fail = 0
$failed = @()
$i = 0
$sw = [Diagnostics.Stopwatch]::StartNew()
foreach ($md in $list) {
    $i++
    $name = [System.IO.Path]::GetFileNameWithoutExtension($md)
    $pdf  = "pdf\$name.pdf"
    $err = & $pandoc $md `
        --pdf-engine=$pdflatex `
        -H $header `
        -V geometry:margin=1in `
        -V fontsize=11pt `
        -V colorlinks=true `
        --highlight-style=tango `
        -o $pdf 2>&1
    $exit = $LASTEXITCODE
    if ($exit -eq 0 -and (Test-Path $pdf)) {
        $ok++
        if ($i % 10 -eq 0 -or $i -eq $total) {
            Write-Host ("[{0}/{1}] OK {2}  (ok={3} fail={4} elapsed={5:N0}s)" -f $i,$total,$name,$ok,$fail,$sw.Elapsed.TotalSeconds)
        }
    } else {
        $fail++
        $failed += $name
        $errMsg = ($err | Out-String).Trim()
        $errLine = ($errMsg -split "`n" | Where-Object { $_ -match "^l\.|^! " } | Select-Object -First 2) -join " | "
        Write-Host ("[{0}/{1}] FAIL {2}  ::  {3}" -f $i,$total,$name,$errLine)
    }
}
$sw.Stop()
""
"==============================="
"Regenerated OK: $ok / $total"
"Failed:         $fail"
"Elapsed:        $([math]::Round($sw.Elapsed.TotalMinutes,1)) min"
$failed | Set-Content -Encoding utf8 "_regen_stale_207_failures.txt"
"Failure list -> _regen_stale_207_failures.txt"
