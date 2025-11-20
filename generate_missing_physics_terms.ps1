# PowerShell script to generate PhysicsTerm wrappers for all 169 standalone modules
# Based on physics_classes_extraction.csv inventory

param(
    [string]$CsvPath = ".\physics_classes_extraction.csv",
    [string]$OutputPath = ".\generated_physics_terms.cpp"
)

Write-Host "📊 Reading physics class inventory from: $CsvPath" -ForegroundColor Cyan

# Read CSV (assuming it exists from subagent task)
if (-Not (Test-Path $CsvPath)) {
    Write-Host "❌ ERROR: $CsvPath not found!" -ForegroundColor Red
    Write-Host "   Creating sample data..." -ForegroundColor Yellow
    
    # Since CSV wasn't created, let's scan the source files directly
    $classes = @()
    
    foreach ($file in Get-ChildItem -Filter "source*.cpp" | Sort-Object {[int]($_.Name -replace '\D',"")}) {
        Write-Host "  Scanning: $($file.Name)" -ForegroundColor Gray
        $content = Get-Content $file.FullName -Raw
        
        # Find class declarations
        $matches = [regex]::Matches($content, '(?ms)^\s*class\s+(\w+(?:Module|Term|Calculator|Solver|Field|Engine|Core|System))\s*(?::|{)')
        
        foreach ($m in $matches) {
            $className = $m.Groups[1].Value
            
            # Skip base classes
            if ($className -match '^(PhysicsTerm|UQFFModule|ModuleInterface)$') { continue }
            
            # Determine base class
            $baseClass = if ($content -match "class\s+$className\s*:\s*public\s+(\w+)") {
                $Matches[1]
            } else {
                "standalone"
            }
            
            # Determine physics type
            $physicsType = switch -Regex ($className) {
                'Vacuum|Quantum|Aether'  { "vacuum" }
                'UQFF|Unified'            { "unified" }
                'Gravity|SMBH|BlackHole' { "gravity" }
                'Magnetic|Magnetar'       { "magnetic" }
                'Galaxy|NGC|M[0-9]+'      { "galactic" }
                'Nebula|Star|Pulsar'      { "stellar" }
                'Hydrogen|Nuclear|LENR'   { "nuclear" }
                'Resonance|Frequency'     { "resonance" }
                default                   { "other" }
            }
            
            $classes += [PSCustomObject]@{
                SourceFile = $file.Name
                ClassName = $className
                LineNumber = 0
                BaseClass = $baseClass
                PhysicsType = $physicsType
            }
        }
    }
    
    Write-Host "`n✅ Found $($classes.Count) unique physics classes" -ForegroundColor Green
    
} else {
    # Import existing CSV
    $classes = Import-Csv $CsvPath
    Write-Host "✅ Loaded $($classes.Count) classes from CSV" -ForegroundColor Green
}

# Generate PhysicsTerm wrappers for standalone modules
$output = "/* AUTO-GENERATED PHYSICS TERMS */`n`n"

# Group classes by physics type
$groupedClasses = $classes | Group-Object PhysicsType

foreach ($group in $groupedClasses | Sort-Object Name) {
    $physicsType = $group.Name
    $typeClasses = $group.Group
    
    $output += "// $($physicsType.ToUpper()) PHYSICS TERMS`n"
    
    foreach ($class in $typeClasses) {
        $className = $class.ClassName
        $sourceFile = $class.SourceFile
        
        # Generate wrapper class name
        $wrapperName = if ($className -notmatch 'Term$') {
            $className -replace 'Module$|Calculator$|Solver$|System$|Core$', 'Term'
        } else {
            $className
        }
        
        # Generate PhysicsTerm wrapper using verbatim escaping
        $cppCode = @'
class {0} : public PhysicsTerm {{
public:
    double compute(double t, const std::map<std::string, double>& params) const override {{
        return 0.0; // TODO: Import from {1} in {2}
    }}
    std::string getName() const override {{ return "{0}"; }}
    std::string getDescription() const override {{ return "{3} from {2}"; }}
}};

'@
        
        $output += ($cppCode -f $wrapperName, $className, $sourceFile, $physicsType)
    }
}

$output += "// END AUTO-GENERATED`n"

