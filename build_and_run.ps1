# Star-Magic UQFF Build and Run Script (Release-MaxCompress)
# Last Updated: November 19, 2025
# Supports: Visual Studio 2022 (primary) and MinGW (alternative)

param(
    [Parameter(Mandatory=$false)]
    [ValidateSet("vs", "mingw", "wolfram")]
    [string]$BuildType = "vs",
    
    [Parameter(Mandatory=$false)]
    [switch]$Clean,
    
    [Parameter(Mandatory=$false)]
    [switch]$Run
)

Write-Host "`n=== Star-Magic UQFF Build System ===" -ForegroundColor Cyan
Write-Host "Build Type: $BuildType" -ForegroundColor Yellow

# Configuration
$VSBuildDir = "build_msvc"
$MinGWBuildDir = "build"
$WolframEnginePath = "C:\Program Files\Wolfram Research\Wolfram Engine\14.1\SystemFiles\Libraries\Windows-x86-64"

# Clean build directories if requested
if ($Clean) {
    Write-Host "`n[CLEAN] Removing build directories..." -ForegroundColor Yellow
    
    if ($BuildType -eq "vs" -or $BuildType -eq "wolfram") {
        Remove-Item -Recurse -Force $VSBuildDir -ErrorAction SilentlyContinue
        Write-Host "  ✓ Removed $VSBuildDir" -ForegroundColor Green
    }
    
    if ($BuildType -eq "mingw") {
        Remove-Item -Recurse -Force $MinGWBuildDir -ErrorAction SilentlyContinue
        Write-Host "  ✓ Removed $MinGWBuildDir" -ForegroundColor Green
    }
}

# Build based on type
switch ($BuildType) {
    "vs" {
        Write-Host "`n[CMAKE] Configuring Visual Studio 2022 (Release-MaxCompress)..." -ForegroundColor Cyan
        cmake -S . -B $VSBuildDir -G "Visual Studio 17 2022" -A x64
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "`n[BUILD] Compiling with Release-MaxCompress optimizations..." -ForegroundColor Cyan
            cmake --build $VSBuildDir --config Release --target MAIN_1_CoAnQi
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host "`n✓ Build successful!" -ForegroundColor Green
                $ExePath = ".\$VSBuildDir\Release\MAIN_1_CoAnQi.exe"
                
                if (Test-Path $ExePath) {
                    $Size = (Get-Item $ExePath).Length / 1MB
                    Write-Host "  Executable: $ExePath" -ForegroundColor White
                    Write-Host "  Size: $([math]::Round($Size, 2)) MB" -ForegroundColor White
                    
                    if ($Run) {
                        Write-Host "`n[RUN] Launching MAIN_1_CoAnQi..." -ForegroundColor Cyan
                        & $ExePath
                    }
                }
            } else {
                Write-Host "`n✗ Build failed!" -ForegroundColor Red
                exit 1
            }
        } else {
            Write-Host "`n✗ CMake configuration failed!" -ForegroundColor Red
            exit 1
        }
    }
    
    "mingw" {
        Write-Host "`n[CMAKE] Configuring MinGW Makefiles..." -ForegroundColor Cyan
        cmake -S . -B $MinGWBuildDir -G "MinGW Makefiles"
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "`n[BUILD] Compiling with MinGW..." -ForegroundColor Cyan
            cmake --build $MinGWBuildDir --target MAIN_1_CoAnQi
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host "`n✓ Build successful!" -ForegroundColor Green
                $ExePath = ".\$MinGWBuildDir\MAIN_1_CoAnQi.exe"
                
                if (Test-Path $ExePath) {
                    $Size = (Get-Item $ExePath).Length / 1MB
                    Write-Host "  Executable: $ExePath" -ForegroundColor White
                    Write-Host "  Size: $([math]::Round($Size, 2)) MB" -ForegroundColor White
                    
                    if ($Run) {
                        Write-Host "`n[RUN] Launching MAIN_1_CoAnQi..." -ForegroundColor Cyan
                        & $ExePath
                    }
                }
            } else {
                Write-Host "`n✗ Build failed!" -ForegroundColor Red
                exit 1
            }
        } else {
            Write-Host "`n✗ CMake configuration failed!" -ForegroundColor Red
            exit 1
        }
    }
    
    "wolfram" {
        Write-Host "`n[CMAKE] Configuring Visual Studio 2022 with Wolfram WSTP..." -ForegroundColor Cyan
        cmake -S . -B $VSBuildDir -G "Visual Studio 17 2022" -A x64 -DUSE_EMBEDDED_WOLFRAM=ON
        
        if ($LASTEXITCODE -eq 0) {
            Write-Host "`n[BUILD] Compiling with Wolfram integration..." -ForegroundColor Cyan
            cmake --build $VSBuildDir --config Release --target MAIN_1_CoAnQi
            
            if ($LASTEXITCODE -eq 0) {
                Write-Host "`n✓ Build successful!" -ForegroundColor Green
                $ExePath = ".\$VSBuildDir\Release\MAIN_1_CoAnQi.exe"
                
                if (Test-Path $ExePath) {
                    $Size = (Get-Item $ExePath).Length / 1MB
                    Write-Host "  Executable: $ExePath" -ForegroundColor White
                    Write-Host "  Size: $([math]::Round($Size, 2)) MB" -ForegroundColor White
                    
                    if ($Run) {
                        Write-Host "`n[RUN] Launching MAIN_1_CoAnQi with Wolfram Engine..." -ForegroundColor Cyan
                        
                        # Add Wolfram DLLs to PATH
                        $env:PATH = "$WolframEnginePath;" + $env:PATH
                        Write-Host "  Added Wolfram Engine to PATH" -ForegroundColor Gray
                        
                        & $ExePath
                    }
                }
            } else {
                Write-Host "`n✗ Build failed!" -ForegroundColor Red
                exit 1
            }
        } else {
            Write-Host "`n✗ CMake configuration failed!" -ForegroundColor Red
            exit 1
        }
    }
}

Write-Host "`n=== Build Complete ===" -ForegroundColor Cyan
Write-Host ""

# Usage examples
Write-Host "Usage Examples:" -ForegroundColor Yellow
Write-Host "  .\build_and_run.ps1 -BuildType vs -Clean -Run          # Clean build with VS2022 and run" -ForegroundColor Gray
Write-Host "  .\build_and_run.ps1 -BuildType mingw                   # Build with MinGW" -ForegroundColor Gray
Write-Host "  .\build_and_run.ps1 -BuildType wolfram -Run            # Build with Wolfram support and run" -ForegroundColor Gray
Write-Host ""
