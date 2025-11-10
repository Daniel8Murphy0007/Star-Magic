#!/usr/bin/env pwsh
# UQFF Computational Engine Launcher v2.1
# Star-Magic Enhanced Edition with 106 Systems
# Sources 4-6 Integration: F_U_Bi_i Framework & Self-Expanding Capabilities

Write-Host "=== UQFF Computational Engine v2.1 - Enhanced Edition ===" -ForegroundColor Cyan
Write-Host "106 Celestial Systems | 7 Source Modules (91-94, 4-6)" -ForegroundColor Yellow
Write-Host "Features: 26-Layer Gravity, F_U_Bi_i Integration, Self-Expanding Framework" -ForegroundColor Yellow
Write-Host ""

# Execute index.js
node index.js @args

# Check exit code
if ($LASTEXITCODE -eq 0) {
    Write-Host "`n✓ Computation complete" -ForegroundColor Green
    Write-Host "23,790 lines | 1.12 MB | v2.1.0" -ForegroundColor Gray
}
else {
    Write-Host "`n✗ Error during execution (Exit code: $LASTEXITCODE)" -ForegroundColor Red
}
