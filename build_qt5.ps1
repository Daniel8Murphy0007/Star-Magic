# Build Qt5 with Perl in PATH
$env:PATH = "C:\Strawberry\perl\bin;$env:PATH"
$env:PERL = "C:\Strawberry\perl\bin\perl.exe"

# Verify Perl accessible
& perl --version
if ($LASTEXITCODE -ne 0) {
    Write-Error "Perl not found in PATH"
    exit 1
}

# Run vcpkg install
& C:\vcpkg\vcpkg.exe install qt5-base:x64-mingw-dynamic --recurse
