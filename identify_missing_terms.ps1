# Create a script to identify the 14 missing Wolfram terms
# by comparing attempted registrations with actual registered terms

$wolframFile = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\wolfram_physics_classes.cpp"
$content = Get-Content $wolframFile -Raw

# Extract all REG("name"... calls
$regPattern = 'REG\("([^"]+)"'
$allAttempted = [regex]::Matches($content, $regPattern) | ForEach-Object { $_.Groups[1].Value }

Write-Host "Total attempted registrations: $($allAttempted.Count)"
Write-Host ""
Write-Host "Creating verification program..."

# Create a small C++ program that will run and tell us which ones registered
$cppCode = @"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

int main() {
    // This list should match all attempted registrations
    std::vector<std::string> attempted = {
"@

foreach ($name in $allAttempted) {
    $escaped = $name -replace '\\', '\\\\' -replace '"', '\"'
    $cppCode += "`n        `"$escaped`","
}

$cppCode += @"

    };
    
    std::cout << "Total attempted: " << attempted.size() << std::endl;
    return 0;
}
"@

Set-Content -Path "verify_registrations.cpp" -Value $cppCode
Write-Host "Created verify_registrations.cpp with $($allAttempted.Count) attempted registrations"
Write-Host ""
Write-Host "Sample of attempted registrations:"
$allAttempted | Select-Object -First 10
