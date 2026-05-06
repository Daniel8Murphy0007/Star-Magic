# Fix all C++ toolset/standard mismatches across the entire project
# Canonical: PlatformToolset=v143 (cl.exe 19.44.35223), LanguageStandard=stdcpp20
# Run from repo root: powershell -ExecutionPolicy Bypass -File _fix_vcxproj_versions.ps1

$fixedCount = 0

# ---- 1. Fix all build/ directory files: v145 → v143 ----
Write-Host "=== Fixing build/ directory: v145 → v143 ==="
Get-ChildItem -Path "build" -Filter "*.vcxproj" -File -Recurse -ErrorAction SilentlyContinue |
    Where-Object { $_.Name -notmatch "VCTargets|CompilerId" } |
    ForEach-Object {
        $c = Get-Content $_.FullName -Raw
        if ($c -match '<PlatformToolset>v145</PlatformToolset>') {
            $c = $c -replace '<PlatformToolset>v145</PlatformToolset>', '<PlatformToolset>v143</PlatformToolset>'
            Set-Content -Path $_.FullName -Value $c -NoNewline
            Write-Host "  FIXED: $($_.Name)"
            $fixedCount++
        }
    }

# ---- 2. Fix test_core_integration.vcxproj: stdcpp17 → stdcpp20 ----
Write-Host "`n=== Fixing test_core_integration.vcxproj: stdcpp17 → stdcpp20 ==="
$f = "test_core_integration.vcxproj"
if (Test-Path $f) {
    $c = Get-Content $f -Raw
    if ($c -match '<LanguageStandard>stdcpp17</LanguageStandard>') {
        $c = $c -replace '<LanguageStandard>stdcpp17</LanguageStandard>', '<LanguageStandard>stdcpp20</LanguageStandard>'
        Set-Content -Path $f -Value $c -NoNewline
        Write-Host "  FIXED: $f"
        $fixedCount++
    } else {
        Write-Host "  Already OK or not found: $f"
    }
}

# ---- 3. Fix RUN_TESTS.vcxproj: add missing LanguageStandard tag ----
Write-Host "`n=== Fixing RUN_TESTS.vcxproj: add LanguageStandard=stdcpp20 ==="
$f = "RUN_TESTS.vcxproj"
if (Test-Path $f) {
    $c = Get-Content $f -Raw
    if ($c -notmatch '<LanguageStandard>') {
        # Insert after <WarningLevel> or <Optimization> in the first ClCompile ItemDefinitionGroup
        if ($c -match '(<WarningLevel>[^<]+</WarningLevel>)') {
            $insertAfter = $Matches[0]
            $c = $c.Replace($insertAfter, "$insertAfter`n      <LanguageStandard>stdcpp20</LanguageStandard>")
            Set-Content -Path $f -Value $c -NoNewline
            Write-Host "  FIXED: $f (added after WarningLevel)"
            $fixedCount++
        } else {
            Write-Host "  WARNING: Could not find insertion point in $f — check manually"
        }
    } else {
        Write-Host "  Already has LanguageStandard: $f"
    }
}

Write-Host "`n=== DONE. Fixed $fixedCount files. ==="
