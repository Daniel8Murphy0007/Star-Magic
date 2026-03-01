# ============================================================================
# analyze_disassembly.ps1 - UQFF Assembly Analysis Script
# 
# Purpose: Analyze compiled MAIN_1_CoAnQi.exe disassembly for:
#   1. SIMD instruction usage (SSE2/AVX)
#   2. Vectorized physics computation patterns
#   3. Optimization opportunities
#   4. CPU cycle estimation
#
# Source: Grok opcode analysis (MIPS/x86-64)
# Created: 2026-03-01
# ============================================================================

param(
    [string]$ExePath = ".\build_msvc\Release\MAIN_1_CoAnQi.exe",
    [string]$OutputFile = "disassembly_analysis.txt",
    [switch]$Detailed,
    [switch]$SIMDOnly
)

Write-Host "============================================" -ForegroundColor Cyan
Write-Host "  UQFF Disassembly Analysis Tool" -ForegroundColor Cyan
Write-Host "============================================" -ForegroundColor Cyan

# Check if dumpbin is available
$dumpbin = Get-Command dumpbin -ErrorAction SilentlyContinue
if (-not $dumpbin) {
    Write-Host "`nSearching for Visual Studio dumpbin..." -ForegroundColor Yellow
    
    # Search common VS locations
    $vsLocations = @(
        "${env:ProgramFiles}\Microsoft Visual Studio\2022\*\VC\Tools\MSVC\*\bin\Hostx64\x64\dumpbin.exe",
        "${env:ProgramFiles(x86)}\Microsoft Visual Studio\2022\*\VC\Tools\MSVC\*\bin\Hostx64\x64\dumpbin.exe",
        "${env:ProgramFiles}\Microsoft Visual Studio\2019\*\VC\Tools\MSVC\*\bin\Hostx64\x64\dumpbin.exe"
    )
    
    foreach ($pattern in $vsLocations) {
        $found = Get-ChildItem -Path $pattern -ErrorAction SilentlyContinue | Select-Object -First 1
        if ($found) {
            $dumpbin = $found.FullName
            Write-Host "Found: $dumpbin" -ForegroundColor Green
            break
        }
    }
    
    if (-not $dumpbin) {
        Write-Host "ERROR: dumpbin.exe not found. Please run from VS Developer Command Prompt." -ForegroundColor Red
        exit 1
    }
}

# Verify executable exists
if (-not (Test-Path $ExePath)) {
    Write-Host "ERROR: Executable not found: $ExePath" -ForegroundColor Red
    Write-Host "Please build the project first with: cmake --build build_msvc --config Release" -ForegroundColor Yellow
    exit 1
}

Write-Host "`nAnalyzing: $ExePath" -ForegroundColor White

# ============================================================================
# SIMD Instruction Patterns
# ============================================================================

$SIMDPatterns = @{
    # SSE2 packed double operations
    "MOVAPD"  = @{ Cycles = 1;  Category = "SSE2-Load";  Description = "Move aligned packed double (128-bit)" }
    "MOVUPD"  = @{ Cycles = 1;  Category = "SSE2-Load";  Description = "Move unaligned packed double" }
    "ADDPD"   = @{ Cycles = 4;  Category = "SSE2-Math";  Description = "Add packed double" }
    "SUBPD"   = @{ Cycles = 4;  Category = "SSE2-Math";  Description = "Subtract packed double" }
    "MULPD"   = @{ Cycles = 4;  Category = "SSE2-Math";  Description = "Multiply packed double" }
    "DIVPD"   = @{ Cycles = 13; Category = "SSE2-Math";  Description = "Divide packed double" }
    "SQRTPD"  = @{ Cycles = 18; Category = "SSE2-Math";  Description = "Square root packed double" }
    "MAXPD"   = @{ Cycles = 4;  Category = "SSE2-Comp";  Description = "Maximum packed double" }
    "MINPD"   = @{ Cycles = 4;  Category = "SSE2-Comp";  Description = "Minimum packed double" }
    "CMPPD"   = @{ Cycles = 4;  Category = "SSE2-Comp";  Description = "Compare packed double" }
    
    # SSE2 scalar double operations
    "MOVSD"   = @{ Cycles = 1;  Category = "SSE2-Scalar"; Description = "Move scalar double" }
    "ADDSD"   = @{ Cycles = 4;  Category = "SSE2-Scalar"; Description = "Add scalar double" }
    "MULSD"   = @{ Cycles = 4;  Category = "SSE2-Scalar"; Description = "Multiply scalar double" }
    "DIVSD"   = @{ Cycles = 13; Category = "SSE2-Scalar"; Description = "Divide scalar double" }
    "SQRTSD"  = @{ Cycles = 18; Category = "SSE2-Scalar"; Description = "Square root scalar double" }
    
    # AVX packed double operations (256-bit)
    "VMOVAPD" = @{ Cycles = 1;  Category = "AVX-Load";  Description = "Move aligned packed double (256-bit)" }
    "VMOVUPD" = @{ Cycles = 1;  Category = "AVX-Load";  Description = "Move unaligned packed double" }
    "VADDPD"  = @{ Cycles = 4;  Category = "AVX-Math";  Description = "Add packed double (4-way)" }
    "VSUBPD"  = @{ Cycles = 4;  Category = "AVX-Math";  Description = "Subtract packed double" }
    "VMULPD"  = @{ Cycles = 4;  Category = "AVX-Math";  Description = "Multiply packed double (4-way)" }
    "VDIVPD"  = @{ Cycles = 13; Category = "AVX-Math";  Description = "Divide packed double" }
    "VSQRTPD" = @{ Cycles = 18; Category = "AVX-Math";  Description = "Square root packed double" }
    
    # AVX FMA (Fused Multiply-Add)
    "VFMADD132PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-add (a*b+c)" }
    "VFMADD213PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-add" }
    "VFMADD231PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-add" }
    "VFMSUB132PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-subtract" }
    "VFMSUB213PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-subtract" }
    "VFMSUB231PD" = @{ Cycles = 4; Category = "AVX-FMA"; Description = "Fused multiply-subtract" }
    
    # Special operations
    "XORPD"   = @{ Cycles = 1;  Category = "SSE2-Logic"; Description = "XOR packed double (zero registers)" }
    "VXORPD"  = @{ Cycles = 1;  Category = "AVX-Logic";  Description = "XOR packed double (256-bit)" }
    "ANDPD"   = @{ Cycles = 1;  Category = "SSE2-Logic"; Description = "AND packed double" }
    "VANDPD"  = @{ Cycles = 1;  Category = "AVX-Logic";  Description = "AND packed double (256-bit)" }
}

# Physics-related patterns to search for
$PhysicsPatterns = @(
    "vacuum",
    "density", 
    "buoyancy",
    "quantum",
    "gravity",
    "magnetic",
    "dipole",
    "cosmological",
    "hubble"
)

# ============================================================================
# Run Disassembly
# ============================================================================

Write-Host "`n[1/5] Running dumpbin disassembly..." -ForegroundColor Yellow

try {
    if ($dumpbin -is [string]) {
        $disasm = & $dumpbin /disasm:nobytes $ExePath 2>&1
    } else {
        $disasm = & dumpbin /disasm:nobytes $ExePath 2>&1
    }
}
catch {
    Write-Host "ERROR: Disassembly failed: $_" -ForegroundColor Red
    exit 1
}

$totalLines = ($disasm | Measure-Object).Count
Write-Host "  Total lines: $totalLines" -ForegroundColor Gray

# ============================================================================
# Analyze SIMD Instructions
# ============================================================================

Write-Host "`n[2/5] Analyzing SIMD instructions..." -ForegroundColor Yellow

$simdCounts = @{}
$categoryTotals = @{}
$estimatedCycles = 0

foreach ($line in $disasm) {
    foreach ($pattern in $SIMDPatterns.Keys) {
        if ($line -match "\b$pattern\b") {
            $simdCounts[$pattern] = ($simdCounts[$pattern] ?? 0) + 1
            
            $info = $SIMDPatterns[$pattern]
            $categoryTotals[$info.Category] = ($categoryTotals[$info.Category] ?? 0) + 1
            $estimatedCycles += $info.Cycles
        }
    }
}

# ============================================================================
# Output Results
# ============================================================================

Write-Host "`n[3/5] Generating report..." -ForegroundColor Yellow

$report = @"
============================================================================
  UQFF Disassembly Analysis Report
  Generated: $(Get-Date -Format "yyyy-MM-dd HH:mm:ss")
  Executable: $ExePath
============================================================================

SIMD INSTRUCTION SUMMARY
------------------------

"@

# Category summary
$report += "`n=== By Category ===`n"
foreach ($cat in ($categoryTotals.Keys | Sort-Object)) {
    $count = $categoryTotals[$cat]
    $report += "  $cat : $count`n"
}

# Detailed instruction counts
$report += "`n=== By Instruction ===`n"
foreach ($instr in ($simdCounts.Keys | Sort-Object)) {
    $count = $simdCounts[$instr]
    $info = $SIMDPatterns[$instr]
    $cycles = $count * $info.Cycles
    $report += "  $($instr.PadRight(15)) : $($count.ToString().PadLeft(6)) occurrences | ~$cycles cycles | $($info.Description)`n"
}

$totalSIMD = ($simdCounts.Values | Measure-Object -Sum).Sum
$report += "`n=== Totals ===`n"
$report += "  Total SIMD instructions: $totalSIMD`n"
$report += "  Estimated cycle cost:    $estimatedCycles`n"

# Vectorization ratio
$totalInstructions = ($disasm | Where-Object { $_ -match '^\s+[0-9A-F]+:' } | Measure-Object).Count
if ($totalInstructions -gt 0) {
    $vectorRatio = [math]::Round(($totalSIMD / $totalInstructions) * 100, 2)
    $report += "  Vectorization ratio:     $vectorRatio%`n"
}

# ============================================================================
# Physics Function Analysis
# ============================================================================

Write-Host "`n[4/5] Identifying physics functions..." -ForegroundColor Yellow

$physicsFunctions = @()
$currentFunction = ""

foreach ($line in $disasm) {
    # Function labels
    if ($line -match '^\s*([A-Za-z_][A-Za-z0-9_:]+):') {
        $currentFunction = $Matches[1]
        
        foreach ($pattern in $PhysicsPatterns) {
            if ($currentFunction -match $pattern) {
                $physicsFunctions += $currentFunction
                break
            }
        }
    }
}

$report += "`n============================================================================`n"
$report += "  PHYSICS-RELATED FUNCTIONS`n"
$report += "============================================================================`n"

$uniquePhysics = $physicsFunctions | Select-Object -Unique
foreach ($func in $uniquePhysics | Sort-Object) {
    $report += "  $func`n"
}
$report += "  Total: $($uniquePhysics.Count) physics functions`n"

# ============================================================================
# Optimization Recommendations
# ============================================================================

Write-Host "`n[5/5] Generating recommendations..." -ForegroundColor Yellow

$report += "`n============================================================================`n"
$report += "  OPTIMIZATION RECOMMENDATIONS`n"
$report += "============================================================================`n"

# Check for scalar operations that could be vectorized
$scalarOps = ($categoryTotals["SSE2-Scalar"] ?? 0)
$packedOps = ($categoryTotals["SSE2-Math"] ?? 0) + ($categoryTotals["AVX-Math"] ?? 0)

if ($scalarOps -gt $packedOps * 2) {
    $report += "`n[!] HIGH scalar-to-vector ratio detected.`n"
    $report += "    Consider using packed operations for multi-system calculations.`n"
    $report += "    Scalar ops: $scalarOps, Packed ops: $packedOps`n"
}

# Check for FMA usage
$fmaOps = ($categoryTotals["AVX-FMA"] ?? 0)
if ($fmaOps -eq 0) {
    $report += "`n[!] No FMA instructions detected.`n"
    $report += "    Enable /arch:AVX2 compiler flag for fused multiply-add operations.`n"
    $report += "    FMA reduces U_bi = -β × Ug × Ω × cos(πtn) to single instruction.`n"
}

# Check for division hotspots
$divOps = ($simdCounts["DIVPD"] ?? 0) + ($simdCounts["VDIVPD"] ?? 0) + ($simdCounts["DIVSD"] ?? 0)
if ($divOps -gt 100) {
    $report += "`n[!] High division count: $divOps`n"
    $report += "    Consider precomputing 1/ρ_UA for density ratio calculations.`n"
    $report += "    Division costs ~13 cycles vs multiplication ~4 cycles.`n"
}

# Check for sqrt hotspots
$sqrtOps = ($simdCounts["SQRTPD"] ?? 0) + ($simdCounts["VSQRTPD"] ?? 0) + ($simdCounts["SQRTSD"] ?? 0)
if ($sqrtOps -gt 50) {
    $report += "`n[!] High sqrt count: $sqrtOps`n"
    $report += "    Consider RSQRTPS (reciprocal sqrt) for less precise applications.`n"
}

if ($fmaOps -gt 0 -and $scalarOps -lt $packedOps) {
    $report += "`n[✓] Good vectorization detected with FMA support.`n"
}

$report += "`n============================================================================`n"
$report += "  END OF REPORT`n"
$report += "============================================================================`n"

# Save report
$report | Out-File -FilePath $OutputFile -Encoding UTF8
Write-Host "`nReport saved to: $OutputFile" -ForegroundColor Green

# Console summary
Write-Host "`n============================================" -ForegroundColor Cyan
Write-Host "  SUMMARY" -ForegroundColor Cyan
Write-Host "============================================" -ForegroundColor Cyan
Write-Host "  Total SIMD: $totalSIMD" -ForegroundColor White
Write-Host "  SSE2 ops:   $(($categoryTotals.Keys | Where-Object { $_ -like 'SSE2-*' } | ForEach-Object { $categoryTotals[$_] } | Measure-Object -Sum).Sum)" -ForegroundColor White
Write-Host "  AVX ops:    $(($categoryTotals.Keys | Where-Object { $_ -like 'AVX-*' } | ForEach-Object { $categoryTotals[$_] } | Measure-Object -Sum).Sum)" -ForegroundColor White
Write-Host "  Est cycles: $estimatedCycles" -ForegroundColor White
Write-Host "============================================`n" -ForegroundColor Cyan
