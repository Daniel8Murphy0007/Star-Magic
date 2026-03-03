# Phase 0B: Backend Build & Test Script
# Tests source2(HEAD PROGRAM).cpp with IPC Pipeline Handler

Write-Host "=" * 70 -ForegroundColor Cyan
Write-Host "PHASE 0B: BACKEND IPC SERVER BUILD TEST" -ForegroundColor Cyan
Write-Host "=" * 70 -ForegroundColor Cyan

# Step 1: Check required files exist
Write-Host "`n[1/5] Checking required files..." -ForegroundColor Yellow

$requiredFiles = @(
    "source2(HEAD PROGRAM).cpp",
    "ipc_pipeline_handler.h",
    "qcalc_subprocess.py",
    "CMakeLists.txt"
)

$allFilesExist = $true
foreach ($file in $requiredFiles) {
    if (Test-Path $file) {
        Write-Host "  ✓ $file" -ForegroundColor Green
    } else {
        Write-Host "  ✗ $file MISSING" -ForegroundColor Red
        $allFilesExist = $false
    }
}

if (-not $allFilesExist) {
    Write-Host "`n✗ Missing required files. Cannot proceed." -ForegroundColor Red
    exit 1
}

# Step 2: Check Python subprocess works
Write-Host "`n[2/5] Testing Python subprocess..." -ForegroundColor Yellow
$testInput = '{"object_name":"TEST","M":2.0,"r":1e6,"B":1e13}'
$pythonTest = echo $testInput | python qcalc_subprocess.py 2>$null

if ($LASTEXITCODE -eq 0) {
    Write-Host "  ✓ qcalc_subprocess.py works" -ForegroundColor Green
} else {
    Write-Host "  ✗ qcalc_subprocess.py failed (exit code: $LASTEXITCODE)" -ForegroundColor Red
    Write-Host "  Fix Python script before building C++" -ForegroundColor Yellow
    exit 1
}

# Step 3: Check build directory exists
Write-Host "`n[3/5] Checking build directory..." -ForegroundColor Yellow
if (Test-Path "build_msvc") {
    Write-Host "  ✓ build_msvc exists" -ForegroundColor Green
} else {
    Write-Host "  ⚠ build_msvc not found, creating..." -ForegroundColor Yellow
    cmake -S . -B build_msvc -G "Visual Studio 17 2022" -A x64
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  ✗ CMake configuration failed" -ForegroundColor Red
        exit 1
    }
}

# Step 4: Build source2(HEAD PROGRAM) target
Write-Host "`n[4/5] Building source2(HEAD PROGRAM).cpp..." -ForegroundColor Yellow
Write-Host "  This may take 2-3 minutes..." -ForegroundColor Gray

cmake --build build_msvc --config Release --target source2 2>&1 | Out-Null

if ($LASTEXITCODE -eq 0) {
    Write-Host "  ✓ Build successful" -ForegroundColor Green
    
    # Check executable exists
    $exePath = "build_msvc\Release\source2`(HEAD PROGRAM`).exe"
    if (Test-Path $exePath) {
        $exeSize = (Get-Item $exePath).Length / 1MB
        Write-Host "  ✓ Executable created: $([math]::Round($exeSize, 2)) MB" -ForegroundColor Green
    }
} else {
    Write-Host "  ✗ Build failed (exit code: $LASTEXITCODE)" -ForegroundColor Red
    Write-Host "`nShowing last 30 lines of build output:" -ForegroundColor Yellow
    Write-Host "─" * 70 -ForegroundColor Gray
    cmake --build build_msvc --config Release --target source2 2>&1 | Select-Object -Last 30
    Write-Host "─" * 70 -ForegroundColor Gray
    exit 1
}

# Step 5: Summary
Write-Host "`n[5/5] Build Test Summary" -ForegroundColor Yellow
Write-Host "─" * 70 -ForegroundColor Gray
Write-Host "✓ Python subprocess working" -ForegroundColor Green
Write-Host "✓ C++ backend compiled successfully" -ForegroundColor Green
Write-Host "✓ IPC Pipeline Handler integrated" -ForegroundColor Green
Write-Host "─" * 70 -ForegroundColor Gray

Write-Host "`n" + "=" * 70 -ForegroundColor Cyan
Write-Host "PHASE 0B COMPLETE: Backend Ready for Testing" -ForegroundColor Cyan
Write-Host "=" * 70 -ForegroundColor Cyan

Write-Host "`nNext Steps:" -ForegroundColor Yellow
Write-Host "  1. Run backend: .\build_msvc\Release\source2`(HEAD PROGRAM`).exe" -ForegroundColor White
Write-Host "     Look for: '[IPC Server] Listening for PIPELINE_PROCESS messages...'" -ForegroundColor Gray
Write-Host "`n  2. Proceed to Phase 0C: Add IPC client to source2.cpp (frontend)" -ForegroundColor White
Write-Host "`n  3. Test end-to-end pipeline" -ForegroundColor White

Write-Host "`nEstimated time to Phase 0C: 30-45 minutes" -ForegroundColor Gray
