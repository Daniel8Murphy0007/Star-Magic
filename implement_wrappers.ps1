# implement_wrappers.ps1
# Generates functional PhysicsTerm wrapper implementations that call actual module methods

Write-Host "`n=== Generating Functional PhysicsTerm Wrapper Implementations ===" -ForegroundColor Cyan
Write-Host "Reading module mappings..." -ForegroundColor Yellow

# Load the JSON mapping
$methodMapping = Get-Content "original_module_methods.json" | ConvertFrom-Json
$catalog = Import-Csv "module_catalog.csv"

# Create wrapper implementations
$implementations = @()
$implementations += "// ============================================"
$implementations += "// FUNCTIONAL PhysicsTerm Wrapper Classes"
$implementations += "// Generated: $(Get-Date -Format 'yyyy-MM-dd HH:mm:ss')"
$implementations += "// Each wrapper calls its embedded module's compute method"
$implementations += "// ============================================`n"

$count = 0
foreach ($row in $catalog) {
    $sourceNum = $row.Source
    $className = $row.ClassName
    $sourceKey = "Source$sourceNum"
    
    # Get method info from JSON
    $methodInfo = $methodMapping.$sourceKey
    
    if ($methodInfo) {
        $originalClass = $methodInfo.Class
        $computeMethod = $methodInfo.Method
        
        $wrapperName = "${className}Term"
        
        $implementations += "// Wrapper for $className (Source$sourceNum)"
        $implementations += "class $wrapperName : public PhysicsTerm {"
        $implementations += "private:"
        $implementations += "    $originalClass instance;"
        $implementations += "public:"
        $implementations += "    ${wrapperName}() {}"
        $implementations += "    "
        $implementations += "    double compute(double t, const std::map<std::string, double>& params) const override {"
        $implementations += "        // Call embedded module's actual compute method"
        $implementations += "        return instance.$computeMethod(t);"
        $implementations += "    }"
        $implementations += "    "
        $implementations += "    std::string getName() const override {"
        $implementations += "        return `"$className`";"
        $implementations += "    }"
        $implementations += "    "
        $implementations += "    std::string getDescription() const override {"
        $desc = "Source${sourceNum}: ${originalClass}.${computeMethod}()"
        $implementations += "        return `"$desc`";"
        $implementations += "    }"
        $implementations += "};"
        $implementations += ""
        
        $count++
        if ($count % 10 -eq 0) {
            Write-Host "  Generated $count wrappers..." -ForegroundColor Gray
        }
    }
    else {
        Write-Host "  WARNING: No method mapping for Source$sourceNum ($className)" -ForegroundColor Red
    }
}

# Save to file
$outFile = "functional_physics_terms.cpp"
$implementations | Out-File $outFile -Encoding UTF8

Write-Host "`n=== Generation Complete ===" -ForegroundColor Green
Write-Host "  Total wrappers: $count"
Write-Host "  Output file: $outFile"
Write-Host "  File size: $([math]::Round((Get-Item $outFile).Length/1KB, 1)) KB"
Write-Host "`nEach wrapper now calls its embedded module's actual compute method!" -ForegroundColor Cyan
