Get-ChildItem -Recurse -Filter "*.vcxproj" -File |
    Where-Object { $_.DirectoryName -ne (Get-Location).Path -and $_.Name -notmatch "ZERO_CHECK|ALL_BUILD|INSTALL|VCTargets|CompilerId|RUN_TESTS" } |
    ForEach-Object {
        $c = Get-Content $_.FullName -Raw
        $ts = if ($c -match '<PlatformToolset>([^<]+)') { $Matches[1].Trim() } else { 'NONE' }
        $lg = if ($c -match '<LanguageStandard>([^<]+)') { $Matches[1].Trim() } else { 'NONE' }
        if ($ts -ne 'v143' -or ($lg -ne 'stdcpp20' -and $lg -ne 'NONE')) {
            "MISMATCH: $($_.Name)  [$($_.Directory.Name)]  toolset=$ts  lang=$lg"
        } else {
            "OK:       $($_.Name)  [$($_.Directory.Name)]  toolset=$ts  lang=$lg"
        }
    }
