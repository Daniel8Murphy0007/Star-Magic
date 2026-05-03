$pandoc = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
# arXiv approved engine: pdflatex (NOT xelatex)
$pdflatex = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
$srcDir = "whitepapers"
$outDir = "pdf"

if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

$papers = @(
  "PAPER_447_OrionNebula_UQFF_MUGE_HAlpha_SFR_Evolution.md",
  "PAPER_448_MultiSystem_UQFF_Core_Compression_Framework_Fenv_Architecture.md",
  "PAPER_449_YoungStars_Outflows_Pressure_Vout100kms_MUGE_UQFF.md",
  "PAPER_450_Eagle_Nebula_Wind_Radiation_Pressure_NGC6611_UQFF.md",
  "PAPER_451_BigBang_Gravity_Evolution_MUGE_QG_DM_GW_Composite_Fcosmo.md",
  "PAPER_452_Compressed_UQFF_Env_Modular_7System_Cycle2_Registry.md",
  "PAPER_453_Magnetar_DualMode_UQFF_Compressed_Frequency_SGR1745_AetherResonance.md",
  "PAPER_454_MultiSystem_Compression_Cycle2_19System_Registry_Universal.md",
  "PAPER_455_UQFF_29System_Expanded_Registry_Saturn_Hydrogen_HRes_Session115Hub.md",
  "PAPER_456_MUGE_29System_Compressed_Unified_Gravity_Duniverse_4Factor_13Fenv.md",
  "PAPER_457_MUGE_38System_Extended_Env_Ftorque_Fshock_Fcosmo_AutoCascade.md",
  "PAPER_458_MUGE_Final_7System_Resonance_10Term_Acceleration_Suite_getSolution.md",
  "PAPER_459_UFE_Orb_Plasmoid_Dynamics_RedDwarf_tMinus_Transform_26QuantumLevels.md",
  "PAPER_460_Nebular_UQFF_Drawing32_LENR_NonLocal_Higgs_DNA_Energy.md",
  "PAPER_461_RedDwarf_LENR_Basel_Pi_Series_S2_Wmag_Cyclotron_Buoyancy.md",
  "PAPER_462_Inertia_UQFF_Wave_Energy_Inertial_Operator_ThreeLeg_Proofset.md",
  "PAPER_463_Hydrogen_Compressed_Space_Espace_7Factor_HiggsFreq_MayanPrecession.md",
  "PAPER_464_M51_WhirlpoolGalaxy_MUGE_UQFF_Tidal_NGC5195_DensityWaves.md",
  "PAPER_465_NGC1316_CosmicDustBunnies_MUGE_UQFF_Merger_AGNJets_ClusterDisruption.md",
  "PAPER_466_V838Mon_LightEcho_UQFF_Ug1_DustModulation_TRZ_VacuumCorrection.md",
  "PAPER_467_NGC1300_BarredSpiral_MUGE_UQFF_BarGasFunneling_DensityWaves_SFR.md",
  "PAPER_468_SMBHBinary_MUGE_UQFF_FrequencyDerived_DPM_THz_Coalescence.md",
  "PAPER_469_NGC346_Nebula_MUGE_UQFF_Ug3Collapse_ClusterEntanglement_Blueshifted.md",
  "PAPER_470_SMBH_Msigma_UQFF_Resonance_FeedbackCalibration_f0063_26StateModel.md",
  "PAPER_471_LENR_Keta_NeutronProduction_Calibration_NonLocal_SSq_26D_UmMediated.md",
  "PAPER_472_Abell2256_UQFF_FUBii_Enhanced.md",
  "PAPER_473_MUGEModule_7System_Compressed_Resonance.md",
  "PAPER_474_MUGEResonanceModule_12System_Superconductive.md",
  "PAPER_475_UQFF_SubTerm_Modules_Catalogue.md",
  "PAPER_476_DPM_PreBigBang_26Sphere_Birth_Model.md",
  "PAPER_477_Buoyancy_Coupling_Constants_Beta_i.md",
  "PAPER_478_Aether_Coupling_Background_Field.md",
  "PAPER_479_UQFFBuoyancyAstroModule_ComplexArithmetic_5System.md",
  "PAPER_480_UQFFBuoyancyCNBModule_CosmicNeutrinoBackground_6System.md",
  "PAPER_481_18System_UQFF_Module_Suite_Oct2025_Batch_FUBii_Implementations.md",
  "PAPER_482_HydrogenResonanceUQFFModule_PTOE_Nuclear_Resonance_UnifiedField.md",
  "PAPER_483_MultiSystem_UQFF_Compiler_Modules_4_5_7_8_AstroSystems.md"
)

$total = $papers.Count
$ok = 0; $fail = 0; $failList = @()

foreach ($md in $papers) {
    $n     = $ok + $fail + 1
    $src   = Join-Path $srcDir $md
    $stem  = $md -replace '\.md$', ''
    $pdf   = Join-Path $outDir "$stem.pdf"

    if (!(Test-Path $src)) {
        Write-Host "[$n/$total] SKIP (not found): $md"
        $fail++
        $failList += $md
        continue
    }

    Write-Host -NoNewline "[$n/$total] $md ... "

    $err = & $pandoc $src -o $pdf `
        "--pdf-engine=$pdflatex" `
        -V "geometry:margin=1in" `
        -V "fontsize=11pt" `
        -V "colorlinks=true" `
        --highlight-style=tango `
        2>&1

    if ($LASTEXITCODE -eq 0) {
        $sz = (Get-Item $pdf).Length
        Write-Host "OK ($([math]::Round($sz/1024,0)) KB)"
        $ok++
    } else {
        Write-Host "FAIL"
        Write-Host "  $err"
        $fail++
        $failList += $md
    }
}

Write-Host ""
Write-Host "=== RESULT: $ok OK, $fail FAILED ==="
if ($failList) { $failList | ForEach-Object { Write-Host "  FAILED: $_" } }
