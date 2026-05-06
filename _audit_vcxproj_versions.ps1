$vcxprojs = Get-ChildItem -Path "." -Filter "*.vcxproj" -File |
    Where-Object { $_.Name -notmatch "ZERO_CHECK|ALL_BUILD|INSTALL|VCTargets|CompilerId" }

$toolsets = @{}
$langs = @{}
$stdOpts = @{}
$mismatched = @()

foreach ($f in $vcxprojs) {
    $c = Get-Content $f.FullName -Raw
    $ts = if ($c -match '<PlatformToolset>([^<]+)') { $Matches[1].Trim() } else { 'NONE' }
    $lg = if ($c -match '<LanguageStandard>([^<]+)') { $Matches[1].Trim() } else { 'NONE' }
    $so = if ($c -match '/std:c\+\+(\w+)') { "/std:c++$($Matches[1])" } else { 'NONE' }

    if (-not $toolsets[$ts]) { $toolsets[$ts] = [System.Collections.Generic.List[string]]::new() }
    $toolsets[$ts].Add($f.Name)
    if (-not $langs[$lg]) { $langs[$lg] = 0 } ; $langs[$lg]++
    if (-not $stdOpts[$so]) { $stdOpts[$so] = 0 } ; $stdOpts[$so]++

    # Flag any that are NOT v143 + stdcpp20
    if ($ts -ne 'v143' -or ($lg -ne 'stdcpp20' -and $lg -ne 'NONE')) {
        $mismatched += "$($f.Name)  [toolset=$ts, lang=$lg, /std=$so]"
    }
}

Write-Host "=== PlatformToolset distribution ==="
$toolsets.GetEnumerator() | Sort-Object Key | % { "  $($_.Value.Count)x  $($_.Key)  --  Examples: $($_.Value[0..2] -join ', ')" }

Write-Host "`n=== LanguageStandard distribution ==="
$langs.GetEnumerator() | Sort-Object Key | % { "  $($_.Value)x  $($_.Key)" }

Write-Host "`n=== /std: flag distribution ==="
$stdOpts.GetEnumerator() | Sort-Object Key | % { "  $($_.Value)x  $($_.Key)" }

Write-Host "`n=== Files NOT on v143 / stdcpp20 ($($mismatched.Count) total) ==="
$mismatched | % { "  $_" }
