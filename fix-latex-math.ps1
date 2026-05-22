param(
    [string]$LatexFile
)

# Read the LaTeX file
$content = Get-Content -Path $LatexFile -Raw

# Replace \(...\) with $...$ for pdflatex compatibility
$fixed = $content -replace '\\\\?\(', '$'
$fixed = $fixed -replace '\\\\?\)', '$'

# Write back
Set-Content -Path $LatexFile -Value $fixed

Write-Host "Fixed math delimiters in $LatexFile"
