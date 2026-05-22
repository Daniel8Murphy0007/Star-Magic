# Final YAML fix - use PowerShell string replacement

$papers = Get-ChildItem -Path "whitepapers" -Filter "PAPER_*.md" -File | Sort-Object Name

Write-Host "[Fix YAML Final] Processing $($papers.Count) files...`n"

$fixed = 0

foreach ($paper in $papers) {
    $content = Get-Content -Path $paper.FullName -Raw -Encoding UTF8
    
    # Pattern 1: title: "...$\...$..." → title: "..."
    # Remove $\<command> and $ from title line
    $original = $content
    
    # Find and fix title line (very specific)
    if ($content -match '(?m)^title:\s*"([^"]*\$[^"]*)"') {
        # Extract title and sanitize it
        if ($content -match '(?m)^title:\s*"([^"]*)"') {
            $titleVal = $matches[1]
            
            # Remove $ and \\command
            $cleanTitle = $titleVal -replace '\$', '' -replace '\\[a-zA-Z]+', '' -replace '\s+', ' '
            $cleanTitle = $cleanTitle.Trim()
            
            # Replace in content
            $content = $content -replace '(?m)^title:\s*"[^"]*"', "title: `"$cleanTitle`""
            $fixed++
        }
    }
    
    # Write back if changed
    if ($content -ne $original) {
        Set-Content -Path $paper.FullName -Value $content -Encoding UTF8 -Force
        if ($fixed % 100 -eq 0) {
            Write-Host "  Fixed $fixed/$($papers.Count)..."
        }
    }
}

Write-Host "`n[Fix YAML Final] Completed: Fixed $fixed/$($papers.Count) files"
