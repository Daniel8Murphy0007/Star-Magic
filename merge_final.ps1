# Simple file merger for source82_wolfram.cpp
# Combines SMBH (15 classes) + Virgo (9 classes, excluding duplicate)

$smbhFile = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\source82_wolfram.cpp"
$virgoFile = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\source82_wolfram_VIRGO_EXTRACTION.cpp"
$outputFile = "c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\source82_wolfram_MERGED.cpp"

# Read both files
$smbh = Get-Content $smbhFile -Raw
$virgo = Get-Content $virgoFile -Encoding UTF8 -Raw

# Extract SMBH content (everything from workspace version)
$smbhContent = $smbh

# Extract Virgo classes only (skip header, stop before duplicate SMBHMSigmaRelationTerm)
$virgoStart = $virgo.IndexOf('class VirgoClusterMassTerm')
$virgoEnd = $virgo.LastIndexOf('class SMBHMSigmaRelationTerm')
$virgoClasses = $virgo.Substring($virgoStart, $virgoEnd - $virgoStart).TrimEnd()

# Merge: SMBH + separator + Virgo classes
$separator = @"

// ========================================================================
// VIRGO CLUSTER COSMOLOGICAL TERMS (Added from orphaned branch integration)
// 9 classes: VirgoClusterMass, VirgoClusterICM, VirgoClusterPotential,
//            VirgoClusterDarkMatter, VirgoClusterM87Jet, VirgoTidalStripping,
//            VirgoClusterVirial, VirgoXRay, VirgoVelocityDispersion
// Source: origin/copilot/create-source82-wolfram (Nov 26, 2025)
// ========================================================================

"@

$merged = $smbhContent.TrimEnd() + $separator + $virgoClasses

# Write merged file
$merged | Out-File -FilePath $outputFile -Encoding UTF8 -NoNewline

Write-Host "✓ Merged file created: $outputFile" -ForegroundColor Green
Write-Host "  - SMBH classes: 15" -ForegroundColor Cyan
Write-Host "  - Virgo classes: 9" -ForegroundColor Cyan
Write-Host "  - Total classes: 24" -ForegroundColor Yellow
Write-Host "  - Duplicate SMBHMSigmaRelationTerm excluded (kept workspace version)" -ForegroundColor Gray
