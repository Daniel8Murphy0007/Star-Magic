$pandoc = "$env:LOCALAPPDATA\Pandoc\pandoc.exe"
$xe    = "C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\xelatex.exe"
$outDir = "pdf"
if (!(Test-Path $outDir)) { New-Item -ItemType Directory -Path $outDir | Out-Null }

$papers = @(
  "PAPER_430_SGR0501_4516_Magnetar_PerSystem_MUGE.md",
  "PAPER_431_SGR1745_2900_Complete_PerSystem_MUGE_BHProximity.md",
  "PAPER_432_SgrA_SMBH_PerSystem_MUGE_Accretion_DM_Precession.md",
  "PAPER_433_Tapestry_Starbirth_PerSystem_MUGE_WindFeedback.md",
  "PAPER_434_Westerlund2_PerSystem_MUGE_TauSF2Myr.md",
  "PAPER_435_PillarsOfCreation_PerSystem_MUGE_ErosionCoupling.md",
  "PAPER_436_RingsOfRelativity_PerSystem_MUGE_LensingAmplification.md",
  "PAPER_437_UQFFLearningAssessment_EvB_AdvancementMetric.md",
  "PAPER_438_NGC2525_PerSystem_MUGE_SNMassLoss_BHProximity.md",
  "PAPER_439_NGC3603_PerSystem_MUGE_CavityPressure_DualWind.md",
  "PAPER_440_BubbleNebula_NGC7635_PerSystem_MUGE_GrowingExpansion.md",
  "PAPER_441_AntennaeGalaxies_PerSystem_MUGE_MergerInteractionBoost.md",
  "PAPER_442_HorseheadNebula_PerSystem_MUGE_GrowingErosion5Myr.md",
  "PAPER_443_NGC1275_PerseusA_PerSystem_MUGE_BDecay_FilamentCoupling_CoolingFlow.md",
  "PAPER_444_HUDF_GalaxiesGalore_PerSystem_MUGE_CosmicScale_HighRedshift.md",
  "PAPER_445_NGC1792_StellarForge_PerSystem_MUGE_StarburstWindDominance.md",
  "PAPER_446_UQFFSource10_5ForceFramework_TriadicGravity_FirstPrimaryTextModule.md"
)

$ok = 0; $fail = 0; $failList = @()

foreach ($md in $papers) {
    $stem = $md -replace '\.md$', ''
    $pdf  = Join-Path $outDir "$stem.pdf"
    Write-Host -NoNewline "[$($ok+$fail+1)/17] $md ... "
    $err = & $pandoc $md -o $pdf `
        "--pdf-engine=$xe" `
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
