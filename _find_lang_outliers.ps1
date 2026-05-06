Get-ChildItem -Path "." -Filter "*.vcxproj" -File |
    Where-Object { $_.Name -notmatch "ZERO_CHECK|ALL_BUILD|INSTALL|VCTargets|CompilerId" } |
    ForEach-Object {
        $c = Get-Content $_.FullName -Raw
        $lg = if ($c -match '<LanguageStandard>([^<]+)') { $Matches[1].Trim() } else { 'NONE' }
        if ($lg -ne 'stdcpp20') { "$($_.Name)  =>  lang=$lg" }
    }
