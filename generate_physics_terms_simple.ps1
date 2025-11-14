# Simple Physics Term Generator
# Generates PhysicsTerm wrappers for Source14-162 modules

$outputFile = "generated_physics_terms.cpp"
$catalog = Import-Csv "module_catalog.csv"

# Header
$content = "// AUTO-GENERATED PHYSICS TERMS`n"
$content += "// Generated: $(Get-Date)`n"
$content += "// Total: $($catalog.Count) modules`n`n"

$count = 0

foreach ($module in $catalog) {
    $count++
    $sourceNum = $module.Source
    $className = $module.ClassName
    $termName = "${className}Term"
    
    Write-Host "[$count/$($catalog.Count)] Generating $termName..."
    
    # Build class string line by line
    $classCode = "`n// Source$sourceNum : $className`n"
    $classCode += "class $termName : public PhysicsTerm {`n"
    $classCode += "private:`n"
    $classCode += "    $className instance;`n"
    $classCode += "public:`n"
    $classCode += "    $termName() {}`n"
    $classCode += "    double compute(double t, const std::map<std::string, double>`& params) const override {`n"
    $classCode += "        return 0.0; // TODO: Implement`n"
    $classCode += "    }`n"
    $classCode += "    std::string getName() const override { return `"$className`"; }`n"
    $classCode += "    std::string getDescription() const override { return `"Source$sourceNum`"; }`n"
    $classCode += "};`n"
    
    $content += $classCode
}

# Write all at once
$content | Out-File $outputFile -Encoding UTF8

Write-Host "`nGeneration complete!" -ForegroundColor Green
Write-Host "Output: $outputFile" -ForegroundColor Cyan
Write-Host "Total terms generated: $count" -ForegroundColor Yellow
