$names=@(
'PAPER_1165_UQFF_beta_i_Triangular_Closure',
'PAPER_1162_UQFF_KK_Tower_Mode_By_Mode_Closure',
'PAPER_1168_UQFF_Falsifiable_Predictions_Closed_Lagrangian',
'PAPER_1160_UQFF_F_TRZ_SO5_Closure',
'PAPER_1169_UQFF_Numerical_Confrontation_P1_P5_With_Archival_Data',
'PAPER_1164_UQFF_T22_Moduli_Stabilization_Closure',
'PAPER_1166_UQFF_V_UA_Polynomial_Closure',
'PAPER_1167_UQFF_All_8_Lagrangian_Gaps_Closed_Master_Synthesis',
'PAPER_1161_UQFF_26_Factorial_Pochhammer_Closure',
'PAPER_1163_UQFF_DPM_SO2_LightCone_Closure',
'PAPER_1159_UQFF_Phi_Res_Codimension_Closure',
'SCm_Mizuno_LENR_Transmutation',
'SCm_Rossi_ECat_Variants_Unified',
'PAPER_1178_UQFF_P13_DESI_Y5_w_Second_Derivative',
'PAPER_1179_UQFF_2027_2028_Quadruple_Falsifier',
'PAPER_1176_UQFF_P12_Euclid_Sigma8_R26_Saturation',
'PAPER_1177_UQFF_2027_Joint_Falsifier_Triple',
'SCm_Holmlid_Parkhomov_PonsFleischmann_Upgrade',
'SCm_PonsFleischmann_Derivation',
'PAPER_1175_UQFF_P11_LIGO_O5_Ringdown_Spectral_Offset',
'PAPER_1174_UQFF_Closed_Ledger_Falsifiability_P6_P10',
'PAPER_1170_UQFF_Vacuum_Energy_Ledger_R26_KK_BSFG_Saturation',
'PAPER_1171_UQFF_KK_Regulator_First_Principles_Derivation',
'uqff_production_arxiv',
'PAPER_1180_UQFF_P14_CMB_S4_mu_Distortion',
'PAPER_1172_UQFF_R26_Independent_Re_Derivation_Gauss_Bonnet',
'PAPER_1173_UQFF_KK_Tower_Hbar_Tracked_Derivation',
'PAPER_888_Et_Full_Lagrangian_Unified_Derivation',
'PAPER_882_Expansion_Lagrangian_Euler_Lagrange',
'PAPER_898_Phonon_Lagrangian_Phi_S26_Derivation',
'PAPER_957_Cooper_Pair_Lagrangian',
'PAPER_886_Erosion_Lagrangian_Euler_Lagrange',
'PAPER_1066_UQFF_Lagrangian_Derivation',
'PAPER_907_Stellar_Wind_Buoyancy_Lagrangian_Variation',
'PAPER_894_SCm_Et_Lagrangian_Variation',
'PAPER_503_UQFF_Lagrangian_Wolfram_Export',
'PAPER_1089_Inflation_Buoyancy_Lagrangian',
'PAPER_1095_Horizon_Buoyancy_Lagrangian',
'PAPER_1094_CMB_Buoyancy_Lagrangian',
'PAPER_1090_DarkEnergy_Buoyancy_Lagrangian',
'PAPER_1103_LQG_Buoyancy_Sector_Lagrangian_SpinFoam',
'PAPER_1065_Buoyancy_Lagrangian_EOM',
'PAPER_946_Merger_Lagrangian_Variation'
)
$miss=@(); $stale=@(); $ok=0
foreach($n in $names){
  $pdf = Join-Path 'pdf' ($n + '.pdf')
  $md  = Join-Path 'whitepapers' ($n + '.md')
  if(-not (Test-Path $pdf)){ $miss += $n; continue }
  if(Test-Path $md){
    $mt = (Get-Item $md).LastWriteTime
    $pt = (Get-Item $pdf).LastWriteTime
    if($mt -gt $pt){ $stale += ("{0}  md={1}  pdf={2}" -f $n, $mt.ToString('MM-dd HH:mm'), $pt.ToString('MM-dd HH:mm')); continue }
  }
  $ok++
}
Write-Host ("OK: {0} / {1}" -f $ok, $names.Count)
Write-Host ""
Write-Host ("MISSING PDFs ({0}):" -f $miss.Count)
$miss | ForEach-Object { Write-Host ("  " + $_) }
Write-Host ""
Write-Host ("STALE PDFs ({0}):" -f $stale.Count)
$stale | ForEach-Object { Write-Host ("  " + $_) }
