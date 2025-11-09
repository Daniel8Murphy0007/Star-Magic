# Comprehensive C++ File Inspection Script
# Star-Magic UQFF Codebase Analysis

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "  Star-Magic C++ File Inspection" -ForegroundColor Cyan
Write-Host "========================================" -ForegroundColor Cyan
Write-Host ""

# Get all .cpp files excluding backups
$cppFiles = Get-ChildItem -Path . -Filter *.cpp -File | Where-Object { $_.DirectoryName -notmatch 'module_backups' }

Write-Host "SUMMARY" -ForegroundColor Yellow
Write-Host "-------" -ForegroundColor Yellow
Write-Host "Total C++ files: $($cppFiles.Count)" -ForegroundColor White
Write-Host ""

# Categorize files
$categories = @{
    'MAIN' = @()
    'Source1-50' = @()
    'Source51-100' = @()
    'Source101-162' = @()
    'Enhanced' = @()
    'Test' = @()
    'Other' = @()
}

foreach ($file in $cppFiles) {
    $name = $file.Name
    
    if ($name -match '^MAIN_1\.cpp$') {
        $categories['MAIN'] += $file
    }
    elseif ($name -match '^[Ss]ource([1-9]|[1-4][0-9]|50)\.cpp$') {
        $categories['Source1-50'] += $file
    }
    elseif ($name -match '^[Ss]ource(5[1-9]|[6-9][0-9]|100)\.cpp$') {
        $categories['Source51-100'] += $file
    }
    elseif ($name -match '^[Ss]ource(10[1-9]|1[1-5][0-9]|16[0-2])\.cpp$') {
        $categories['Source101-162'] += $file
    }
    elseif ($name -match '_[Ee]nhanced\.cpp$|_baseline_backup\.cpp$') {
        $categories['Enhanced'] += $file
    }
    elseif ($name -match '^test_|^verify_') {
        $categories['Test'] += $file
    }
    else {
        $categories['Other'] += $file
    }
}

Write-Host "CATEGORIES" -ForegroundColor Yellow
Write-Host "----------" -ForegroundColor Yellow
foreach ($cat in $categories.Keys | Sort-Object) {
    $count = $categories[$cat].Count
    if ($count -gt 0) {
        Write-Host "${cat}: $count files" -ForegroundColor Green
    }
}
Write-Host ""

# Analyze file headers to identify module types
Write-Host "ANALYZING FILE CONTENTS..." -ForegroundColor Yellow
Write-Host ""

$analysis = @()

foreach ($file in $cppFiles | Sort-Object Name) {
    $content = Get-Content $file.FullName -TotalCount 100 -ErrorAction SilentlyContinue
    if (-not $content) { continue }
    
    $info = [PSCustomObject]@{
        File = $file.Name
        Size = $file.Length
        Lines = (Get-Content $file.FullName | Measure-Object -Line).Lines
        HasMain = $false
        HasClass = $false
        ClassName = ""
        HasComplex = $false
        HasUQFF = $false
        Includes = 0
        Comments = ""
    }
    
    # Scan first 100 lines
    foreach ($line in $content) {
        if ($line -match '^\s*int\s+main\s*\(') { $info.HasMain = $true }
        if ($line -match '^\s*class\s+(\w+)') { 
            $info.HasClass = $true 
            if (-not $info.ClassName) {
                $info.ClassName = $matches[1]
            }
        }
        if ($line -match '#include\s*<complex>|std::complex') { $info.HasComplex = $true }
        if ($line -match 'UQFF|Unified.*Quantum.*Field|F_U_Bi_i') { $info.HasUQFF = $true }
        if ($line -match '^\s*#include') { $info.Includes++ }
        if ($line -match '^\s*//\s*(.+)$' -and -not $info.Comments) {
            $info.Comments = $matches[1].Substring(0, [Math]::Min(50, $matches[1].Length))
        }
    }
    
    $analysis += $info
}

# Display detailed analysis
Write-Host "DETAILED FILE ANALYSIS" -ForegroundColor Yellow
Write-Host "----------------------" -ForegroundColor Yellow
Write-Host ""

Write-Host "Main Programs (with main() function):" -ForegroundColor Cyan
$mainPrograms = $analysis | Where-Object { $_.HasMain } | Sort-Object File
foreach ($prog in $mainPrograms) {
    Write-Host "  $($prog.File) - $($prog.Lines) lines, Class: $($prog.ClassName)" -ForegroundColor White
}
Write-Host "  Total: $($mainPrograms.Count)" -ForegroundColor Green
Write-Host ""

Write-Host "Class-Based Modules (no main):" -ForegroundColor Cyan
$classModules = $analysis | Where-Object { $_.HasClass -and -not $_.HasMain } | Sort-Object File
Write-Host "  Total: $($classModules.Count)" -ForegroundColor Green
if ($classModules.Count -gt 0 -and $classModules.Count -le 20) {
    foreach ($mod in $classModules) {
        Write-Host "  $($mod.File) - Class: $($mod.ClassName)" -ForegroundColor White
    }
}
Write-Host ""

Write-Host "UQFF Physics Modules:" -ForegroundColor Cyan
$uqffModules = $analysis | Where-Object { $_.HasUQFF } | Sort-Object File
Write-Host "  Total: $($uqffModules.Count)" -ForegroundColor Green
Write-Host ""

Write-Host "Complex Number Support:" -ForegroundColor Cyan
$complexModules = $analysis | Where-Object { $_.HasComplex } | Sort-Object File
Write-Host "  Total: $($complexModules.Count)" -ForegroundColor Green
Write-Host ""

# Size distribution
Write-Host "SIZE DISTRIBUTION" -ForegroundColor Yellow
Write-Host "-----------------" -ForegroundColor Yellow
$small = ($analysis | Where-Object { $_.Size -lt 20000 }).Count
$medium = ($analysis | Where-Object { $_.Size -ge 20000 -and $_.Size -lt 50000 }).Count
$large = ($analysis | Where-Object { $_.Size -ge 50000 -and $_.Size -lt 100000 }).Count
$xlarge = ($analysis | Where-Object { $_.Size -ge 100000 }).Count

Write-Host "  Small (<20KB):    $small files" -ForegroundColor White
Write-Host "  Medium (20-50KB): $medium files" -ForegroundColor White
Write-Host "  Large (50-100KB): $large files" -ForegroundColor White
Write-Host "  X-Large (>100KB): $xlarge files" -ForegroundColor White
Write-Host ""

# Largest files
Write-Host "LARGEST FILES" -ForegroundColor Yellow
Write-Host "-------------" -ForegroundColor Yellow
$largest = $analysis | Sort-Object Size -Descending | Select-Object -First 10
foreach ($file in $largest) {
    $sizeKB = [Math]::Round($file.Size / 1KB, 1)
    Write-Host "  $($file.File): ${sizeKB}KB ($($file.Lines) lines)" -ForegroundColor White
}
Write-Host ""

# Export to CSV for detailed review
$csvPath = "cpp_file_analysis.csv"
$analysis | Export-Csv -Path $csvPath -NoTypeInformation
Write-Host "Detailed analysis exported to: $csvPath" -ForegroundColor Green
Write-Host ""

# Recently modified files
Write-Host "RECENTLY MODIFIED (Last 7 days)" -ForegroundColor Yellow
Write-Host "--------------------------------" -ForegroundColor Yellow
$recent = $cppFiles | Where-Object { $_.LastWriteTime -gt (Get-Date).AddDays(-7) } | Sort-Object LastWriteTime -Descending
foreach ($file in $recent) {
    Write-Host "  $($file.Name) - $($file.LastWriteTime)" -ForegroundColor Magenta
}
Write-Host "  Total: $($recent.Count)" -ForegroundColor Green
Write-Host ""

Write-Host "========================================" -ForegroundColor Cyan
Write-Host "Analysis Complete!" -ForegroundColor Green
Write-Host "========================================" -ForegroundColor Cyan
