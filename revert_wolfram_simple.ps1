# Temporarily remove error handling to see if the program runs
$filePath = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\wolfram_physics_classes.cpp"
$content = Get-Content $filePath -Raw

# Convert SAFE_REGISTER back to core.registerPhysicsTerm
$pattern = 'SAFE_REGISTER\("([^"]+)",\s*([^,]+),\s*"([^"]+)"\);'
$replacement = 'core.registerPhysicsTerm("$1", std::make_unique<$2>(), "$3");'

$newContent = $content -replace $pattern, $replacement

# Remove the macro definition and debug output
$newContent = $newContent -replace '(?s)std::cout << "DEBUG.*?std::cout\.flush\(\);', ''
$newContent = $newContent -replace '(?s)#define SAFE_REGISTER.*?\}', ''

# Remove the summary logging that uses undefined variables
$newContent = $newContent -replace '(?s)std::cout << "Wolfram Registration Summary.*?error_log\.close\(\);', ''

# Simplify the function start
$newContent = $newContent -replace 'void registerAllWolframPhysicsTerms\(CalculatorCore& core\) \{[^c]*core', 'void registerAllWolframPhysicsTerms(CalculatorCore& core) {
    core'

Set-Content -Path $filePath -Value $newContent -NoNewline

Write-Host "Removed error handling from wolfram_physics_classes.cpp"
Write-Host "File now back to simple core.registerPhysicsTerm calls"
