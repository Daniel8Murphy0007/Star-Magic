# Add MinGW to System PATH
# This script must be run as Administrator

$mingwPath = "C:\MinGW\bin"

# Get current system PATH
$currentPath = [Environment]::GetEnvironmentVariable("Path", "Machine")

# Check if MinGW is already in PATH
if ($currentPath -notlike "*$mingwPath*") {
    Write-Host "Adding $mingwPath to system PATH..." -ForegroundColor Green
    
    # Add MinGW to PATH
    $newPath = $currentPath + ";$mingwPath"
    [Environment]::SetEnvironmentVariable("Path", $newPath, "Machine")
    
    Write-Host "✓ Successfully added MinGW to system PATH!" -ForegroundColor Green
    Write-Host ""
    Write-Host "IMPORTANT: You need to restart VS Code for the changes to take effect." -ForegroundColor Yellow
    Write-Host ""
    Write-Host "To verify after restart, run: g++ --version" -ForegroundColor Cyan
} else {
    Write-Host "✓ MinGW is already in system PATH" -ForegroundColor Green
}

Write-Host ""
Write-Host "Press any key to exit..."
$null = $Host.UI.RawUI.ReadKey("NoEcho,IncludeKeyDown")
