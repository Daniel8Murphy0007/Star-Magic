# Add Wolfram DLL to PATH
$env:PATH = "C:\Program Files\Wolfram Research\Wolfram Engine\14.1\SystemFiles\Libraries\Windows-x86-64;" + $env:PATH

# Run MAIN_1_CoAnQi with option 9 (Wolfram), then 10 (Exit)
Write-Host "=== Testing Wolfram Integration in MAIN_1_CoAnQi ===" -ForegroundColor Cyan
Write-Host "Selecting Option 9 (Wolfram symbolic bridge)..." -ForegroundColor Yellow

$input = @"
9
10
"@

$input | .\build_msvc\Release\MAIN_1_CoAnQi.exe 2>&1 | Tee-Object -FilePath wolfram_output.log

Write-Host "`n=== Searching for Wolfram-related output ===" -ForegroundColor Cyan
Get-Content wolfram_output.log | Select-String -Pattern "Wolfram|kernel|Bridge|WSTP|mathlink|symbolic" -Context 1,1
