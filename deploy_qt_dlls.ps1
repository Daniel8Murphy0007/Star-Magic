# Deploy Qt DLLs for MAIN_1_CoAnQi.exe
# Created: December 1, 2025
# Purpose: Copy all required Qt6 DLLs and plugins for Grok API (source178)

param(
    [string]$BuildDir = "build_msvc\Release",
    [string]$Qt6Path = "C:\Qt\6.10.0\msvc2022_64"
)

Write-Host "Deploying Qt6 DLLs to $BuildDir..." -ForegroundColor Cyan

# Core Qt6 DLLs required by source178_grok_api.cpp (QNetworkAccessManager)
$requiredDlls = @(
    'Qt6Core.dll',
    'Qt6Network.dll',
    'Qt6Widgets.dll',
    'Qt6Gui.dll'
)

# Copy DLLs
foreach ($dll in $requiredDlls) {
    $src = Join-Path "$Qt6Path\bin" $dll
    if (Test-Path $src) {
        Copy-Item $src $BuildDir -Force
        Write-Host "  ✓ Copied $dll" -ForegroundColor Green
    } else {
        Write-Host "  ✗ Missing: $dll" -ForegroundColor Red
    }
}

# Create platforms directory for Qt platform plugins
$platformsDir = Join-Path $BuildDir "platforms"
New-Item -ItemType Directory -Path $platformsDir -Force | Out-Null

# Copy Windows platform plugin (required for Qt to initialize on Windows)
$qwindowsDll = Join-Path "$Qt6Path\plugins\platforms" "qwindows.dll"
if (Test-Path $qwindowsDll) {
    Copy-Item $qwindowsDll $platformsDir -Force
    Write-Host "  ✓ Copied qwindows.dll to platforms\" -ForegroundColor Green
} else {
    Write-Host "  ✗ Missing: qwindows.dll" -ForegroundColor Red
}

# Optional: Copy TLS plugins for HTTPS (Grok API uses HTTPS)
$networkDir = Join-Path $BuildDir "networkinformation"
New-Item -ItemType Directory -Path $networkDir -Force | Out-Null

$tlsPlugins = @(
    'qcertonlybackend.dll',
    'qschannelbackend.dll'
)

foreach ($plugin in $tlsPlugins) {
    $src = Join-Path "$Qt6Path\plugins\tls" $plugin
    if (Test-Path $src) {
        Copy-Item $src (Join-Path $BuildDir "tls") -Force -ErrorAction SilentlyContinue
    }
}

Write-Host "`nDeployment complete!" -ForegroundColor Cyan
Write-Host "Run: .\$BuildDir\MAIN_1_CoAnQi.exe" -ForegroundColor Yellow
