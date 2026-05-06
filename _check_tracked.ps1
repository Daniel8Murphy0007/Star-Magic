$files = @("Source2.vcxproj", "build_check/Source2.vcxproj")
foreach ($path in $files) {
    $c = Get-Content $path -Raw
    $ts = if ($c -match '<PlatformToolset>([^<]+)') { $Matches[1] } else { 'NONE' }
    $lg = if ($c -match '<LanguageStandard>([^<]+)') { $Matches[1] } else { 'NONE' }
    Write-Host "${path}: toolset=$ts  lang=$lg"
}
