# Parse wolfram_physics_terms_FULL.txt into categorized entity lists
# Output: JSON file with entity types and names for C++ class generation

$content = Get-Content "wolfram_physics_terms_FULL.txt" -Raw

# Extract all Entity["Type", "Name"] patterns
$entityPattern = 'Entity\["([^"]+)",\s*"([^"]+)"\]'
$entityWithParamsPattern = 'Entity\["([^"]+)",\s*\{([^}]+)\}\]'

$entities = @{
    "PhysicalConstant" = @()
    "PhysicalQuantity" = @()
    "Particle" = @()
    "Isotope" = @()
    "Star" = @()
    "Galaxy" = @()
    "Exoplanet" = @()
    "PlanetaryMoon" = @()
}

# Parse simple entities: Entity["Type", "Name"]
[regex]::Matches($content, $entityPattern) | ForEach-Object {
    $type = $_.Groups[1].Value
    $name = $_.Groups[2].Value
    
    if ($entities.ContainsKey($type)) {
        if ($entities[$type] -notcontains $name) {
            $entities[$type] += $name
        }
    }
}

# Parse complex entities: Entity["Type", {...}]
[regex]::Matches($content, $entityWithParamsPattern) | ForEach-Object {
    $type = $_.Groups[1].Value
    $params = $_.Groups[2].Value
    
    if ($entities.ContainsKey($type)) {
        # Use full parameter set as unique identifier
        $uniqueName = $params -replace '"', '' -replace ',', '_' -replace '\s+', ''
        if ($entities[$type] -notcontains $uniqueName) {
            $entities[$type] += $uniqueName
        }
    }
}

# Generate summary
Write-Host "`n=== WOLFRAM ENTITY PARSING ===" -ForegroundColor Cyan
$totalParsed = 0
foreach ($type in $entities.Keys | Sort-Object) {
    $count = $entities[$type].Count
    $totalParsed += $count
    if ($count -gt 0) {
        Write-Host "$type : $count entities" -ForegroundColor Green
    }
}
Write-Host "`nTotal parsed: $totalParsed" -ForegroundColor Yellow

# Save to JSON
$output = @{
    "entities" = $entities
    "total_count" = $totalParsed
    "timestamp" = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
}

$output | ConvertTo-Json -Depth 10 | Out-File "wolfram_entities_parsed.json" -Encoding UTF8
Write-Host "Saved to: wolfram_entities_parsed.json" -ForegroundColor Cyan

# Generate initial statistics
Write-Host "`n=== C++ CLASS GENERATION PLAN ===" -ForegroundColor Cyan
Write-Host "Existing classes in MAIN_1_CoAnQi.cpp: 774" -ForegroundColor White
Write-Host "New Wolfram classes to generate: $totalParsed" -ForegroundColor Yellow
Write-Host "Total after integration: $($774 + $totalParsed)" -ForegroundColor Green
