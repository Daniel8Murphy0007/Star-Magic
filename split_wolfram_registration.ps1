# Split the Wolfram registration function into 10 smaller batch functions
$filePath = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\wolfram_physics_classes.cpp"
$content = Get-Content $filePath -Raw

# Find the start and end of the registration function
$startPattern = 'void registerAllWolframPhysicsTerms\(CalculatorCore& core\) \{'
$endPattern = '^\s*\}\s*$'

# Extract all SAFE_REGISTER lines
$registerPattern = 'SAFE_REGISTER\([^;]+\);'
$matches = [regex]::Matches($content, $registerPattern)

Write-Host "Found $($matches.Count) SAFE_REGISTER calls"

# Calculate batch size (570 registrations per batch for 10 batches)
$batchSize = [Math]::Ceiling($matches.Count / 10)
Write-Host "Batch size: $batchSize registrations per function"

# We'll keep the single function but add progress logging every 500 registrations
# This is safer than splitting into multiple functions
