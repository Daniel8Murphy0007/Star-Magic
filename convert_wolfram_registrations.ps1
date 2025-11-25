# Convert all core.registerPhysicsTerm calls to SAFE_REGISTER macro calls
$filePath = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\wolfram_physics_classes.cpp"
$content = Get-Content $filePath -Raw

# Pattern to match: core.registerPhysicsTerm("Name", std::make_unique<NameTerm>(), "Category");
# Replace with: SAFE_REGISTER("Name", NameTerm, "Category");

$pattern = 'core\.registerPhysicsTerm\("([^"]+)",\s*std::make_unique<([^>]+)>\(\),\s*"([^"]+)"\);'
$replacement = 'SAFE_REGISTER("$1", $2, "$3");'

$newContent = $content -replace $pattern, $replacement

# Save the file
Set-Content -Path $filePath -Value $newContent -NoNewline

Write-Host "Converted all registration calls to SAFE_REGISTER macro"
Write-Host "Original file size: $($content.Length) bytes"
Write-Host "New file size: $($newContent.Length) bytes"

# Count conversions
$originalCount = ([regex]::Matches($content, 'core\.registerPhysicsTerm')).Count
$newCount = ([regex]::Matches($newContent, 'SAFE_REGISTER')).Count
Write-Host "Converted $originalCount core.registerPhysicsTerm calls to $newCount SAFE_REGISTER calls"
