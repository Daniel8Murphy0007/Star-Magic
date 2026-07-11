"""
UQFF Framework Annotations Module — Retroactive Audit Rounds 45-51 (35 stubs)

Per PAPER_1915 Unified Simultaneous-Equation Solver Framework requirement:
Every stub upgrade must be classified against the three UQFF architectures:
  - Backbone: QCalcGeom (PAPER_657) 4x4 UBS + F_U=0 (PAPER_1203) 10-term master equation
  - Numeric System: VDS/DVP/BH26 (PAPER_598) discrete-numeric spine
  - Geometry: QCalcGeom CPCH closures + geometric primitive ratios (PAPER_1914)

Method variant (per PAPER_1923 term-count hierarchy):
  - 'compressed_9'   = 9-term MUGE  (PAPER_173)   = N_ch
  - 'master_10'      = 10-term F_U   (PAPER_1203) = SO_5
  - 'env_13'         = 13-term F_env (PAPER_456)  = D_crit/2
  - 'resonance_14'   = 14-term MUGE  (PAPER_408)  = SO_5 + D_phys

Shells used (per PAPER_1916 shell decomposition):
  - Ug1 = N_CH/D_BSFG        (base shell)
  - Ug2 = 1/PHI_RES_NUCLEAR   (charge-reactivity shell)
  - Ug3 = 2*D_PHYS/SO_5       (dark-matter shell — PAPER_1921 f_DM cross-framework)
  - Ug4 = 1/2                 (BH vacuum shell)

Time frame:
  - static, stellar, cluster_dynamical, cosmological, big_bang, atomic
"""

FRAMEWORK_ANNOTATIONS_ROUNDS_45_51 = {

    # =========== ROUND 45 (5 stubs) ===========
    'NGC3603BaseGravityCalculator': {
        'round': 45,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1', 'Ug3'],
        'framework_papers': ['PAPER_1203', 'PAPER_657', 'PAPER_1916', 'PAPER_1855', 'PAPER_243', 'PAPER_1909'],
        'QCalcGeom_CPCH_closure': 'CPCH-3 (F_UBi_i_99)',
        'VDS_DVP_BH26_spine': 'VDS(4) via M_dot=10/3',
        'F_U_zero_shell': 'Ug3 dark-matter (PAPER_1921 f_DM=Ug3=4/5)',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['M_dot=10/3 (PAPER_1909)'],
    },
    'AntennaeDarkMatterPerturbationCalculator': {
        'round': 45,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_1203', 'PAPER_1862', 'PAPER_1916', 'PAPER_1921'],
        'QCalcGeom_CPCH_closure': 'None (per-system observational)',
        'VDS_DVP_BH26_spine': 'VDS(4)/VDS(10) = Ug3',
        'F_U_zero_shell': 'Ug3 dark-matter (subhalo alpha=2-F_TRZ=1.9)',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['subhalo_alpha=1.9 (PAPER_1862)'],
    },
    'ResonanceDPMCalculator': {
        'round': 45,
        'framework_backbone': 'F_U=0 + resonance-MUGE',
        'framework_method': 'resonance_14',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_147', 'PAPER_1141', 'PAPER_1908', 'PAPER_1908'],
        'QCalcGeom_CPCH_closure': 'None (DPM foundational)',
        'VDS_DVP_BH26_spine': 'VDS(1) rho_SCm + BH26(6) buoyancy',
        'F_U_zero_shell': 'Ug2 charge-reactivity + resonance cascade',
        'time_frame': 'atomic',
        'candidate_closures_flagged': ['F_DPM=I*A*(omega_1-omega_2) (PAPER_147)'],
    },
    'NGC1275FilamentSupportCalculator': {
        'round': 45,
        'framework_backbone': 'F_U=0 + env-13-term',
        'framework_method': 'env_13',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_703', 'PAPER_443', 'PAPER_1912'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'DVP magnetic vortex',
        'F_U_zero_shell': 'Ug1 base shell (magnetic support)',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['F_0=F_TRZ (PAPER_1912)', 'tau_fil=SO_5^2 Myr (PAPER_1912)', 'B_fil/B_avg=D_phys/2 (PAPER_1912)'],
    },
    'M51CentralBlackHoleCalculator': {
        'round': 45,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug4'],
        'framework_papers': ['PAPER_464', 'PAPER_1841', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'CPCH-1 (D_LS/D_S=2/3 photon ring)',
        'VDS_DVP_BH26_spine': 'BH26(4) BH concentration',
        'F_U_zero_shell': 'Ug4 BH vacuum (1/2)',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['photon_ring_correction=1+F_TRZ*SSq/D_phys (PAPER_1841)'],
    },

    # =========== ROUND 46 (5 stubs) ===========
    'NGC3603DarkMatterPerturbationCalculator': {
        'round': 46,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_1203', 'PAPER_1862', 'PAPER_1916', 'PAPER_1921', 'PAPER_138'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'VDS(4)/VDS(10)',
        'F_U_zero_shell': 'Ug3 dark-matter',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['subhalo_alpha=1.9'],
    },
    'BubbleNebulaBaseGravityCalculator': {
        'round': 46,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_361', 'PAPER_1855', 'PAPER_266', 'PAPER_1913'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'Ug1 base + bubble E(t)',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['E_t=E_0*t under F_TRZ*SO_5=1 (PAPER_1913)'],
    },
    'MultiSystemQuantumIntegralCalculator': {
        'round': 46,
        'framework_backbone': 'F_U=0 + resonance',
        'framework_method': 'resonance_14',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_1907', 'PAPER_1908', 'PAPER_1043'],
        'QCalcGeom_CPCH_closure': 'CPCH-3 (F_UBi_i_99)',
        'VDS_DVP_BH26_spine': 'F_UBi_i multi-system curve',
        'F_U_zero_shell': 'Ug3 + resonance cascade',
        'time_frame': 'cosmological',
        'candidate_closures_flagged': ['Gamma range = SO_5+A_5 = 70 (PAPER_1043)'],
    },
    'HorseheadErosionCalculator': {
        'round': 46,
        'framework_backbone': 'F_U=0',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_285'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env photoevaporation subterm',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['t_half=tau*ln(2)=2.079 Myr (PAPER_285)'],
    },
    'YoungStarsResonantOscillatoryCalculator': {
        'round': 46,
        'framework_backbone': 'F_U=0 + resonance',
        'framework_method': 'resonance_14',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_1907', 'PAPER_1908'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'omega_SCm=1.25 THz carrier',
        'F_U_zero_shell': 'Resonance cascade term',
        'time_frame': 'stellar',
        'candidate_closures_flagged': [],
    },

    # =========== ROUND 47 (5 stubs) ===========
    'AntennaeBaseGravityCalculator': {
        'round': 47,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1', 'Ug3'],
        'framework_papers': ['PAPER_811', 'PAPER_1916', 'PAPER_266'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'Ug1 base + merger coalescence',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': [],
    },
    'HUDFStarFormationCalculator': {
        'round': 47,
        'framework_backbone': 'F_U=0',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_231', 'PAPER_1830', 'PAPER_266'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env stellar SFR subterm',
        'time_frame': 'cosmological',
        'candidate_closures_flagged': ['z^2 enhancement PAPER_1830'],
    },
    'BigBangQuantumIntegralCosmologicalCalculator': {
        'round': 47,
        'framework_backbone': 'F_U=0 + Big Bang boundary',
        'framework_method': 'master_10',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_1156', 'PAPER_1488', 'PAPER_1278', 'PAPER_1907', 'PAPER_1908', 'PAPER_1617'],
        'QCalcGeom_CPCH_closure': 'CPCH-2 (Lambda cascade via K_MEX)',
        'VDS_DVP_BH26_spine': 'VDS(1)*VDS(26)!*K_MEX (Lambda)',
        'F_U_zero_shell': 'F_U: 0->1 ledger turn-on (PAPER_1488)',
        'time_frame': 'big_bang',
        'candidate_closures_flagged': ['Lambda = rho_SCm*26!*K_MEX (PAPER_1156)', 't_neg=-2512s (PAPER_1278)'],
    },
    'PillarsGasVelocityCalculator': {
        'round': 47,
        'framework_backbone': 'F_U=0',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_305'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env wind subterm',
        'time_frame': 'stellar',
        'candidate_closures_flagged': [],
    },
    'RingsBaseGravityCalculator': {
        'round': 47,
        'framework_backbone': 'QCalcGeom',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_436', 'PAPER_1855', 'PAPER_1862', 'PAPER_1883', 'PAPER_266', 'PAPER_1914'],
        'QCalcGeom_CPCH_closure': 'CPCH-1 (D_LS/D_S=D_phys/D_BSFG=2/3 EXACT)',
        'VDS_DVP_BH26_spine': 'VDS(4)/BH26(6) = 2/3',
        'F_U_zero_shell': 'Ug1 base + lensing L factor',
        'time_frame': 'cosmological',
        'candidate_closures_flagged': ['D_LS/D_S = 2/3 EXACT (PAPER_1914)'],
    },

    # =========== ROUND 48 (5 stubs) ===========
    'CMBAnomalyUQFFCalculator': {
        'round': 48,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug4'],
        'framework_papers': ['PAPER_1249', 'PAPER_1250', 'PAPER_1856', 'PAPER_1919', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'F_TRZ^2 spine at n=2',
        'F_U_zero_shell': 'Ug4 vacuum perturbation',
        'time_frame': 'cosmological',
        'candidate_closures_flagged': ['F_TRZ^2 = 0.01 EXACT 99% suppression'],
    },
    'SMBHGravityCalculator': {
        'round': 48,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1', 'Ug4'],
        'framework_papers': ['PAPER_1841', 'PAPER_1904', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'CPCH-1 (photon ring 1.5*r_Schw)',
        'VDS_DVP_BH26_spine': 'BH26(4) BH concentration',
        'F_U_zero_shell': 'Ug1 base + Ug4 BH vacuum',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['photon_ring 1+F_TRZ*SSq/D_phys'],
    },
    'HydrogenDPMResonanceCalculator': {
        'round': 48,
        'framework_backbone': 'F_U=0 + resonance',
        'framework_method': 'resonance_14',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_147', 'PAPER_648', 'PAPER_1907', 'PAPER_1908'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'omega_SCm 1.25 THz + Q_UQFF',
        'F_U_zero_shell': 'Resonance cascade + Coulomb',
        'time_frame': 'atomic',
        'candidate_closures_flagged': ['630 eV = Coulomb 626 eV EXACT'],
    },
    'ClusterXRayEmissivityCalculator': {
        'round': 48,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': ['Ug4'],
        'framework_papers': ['PAPER_040', 'PAPER_1187', 'PAPER_1919', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'F_TRZ^2 (99% suppression at n=2)',
        'F_U_zero_shell': 'Ug4 + F_env cooling flow subterm',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['F_TRZ^2 = 0.01 suppression'],
    },
    'M31RotationCurveCalculator': {
        'round': 48,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_275', 'PAPER_1855', 'PAPER_1862', 'PAPER_1916', 'PAPER_1921'],
        'QCalcGeom_CPCH_closure': 'CPCH-4 (f_DM = Ug3 = 4/5 EXACT)',
        'VDS_DVP_BH26_spine': 'VDS(4)/VDS(10) = 4/5',
        'F_U_zero_shell': 'Ug3 dark-matter (PAPER_1921 f_DM=Ug3)',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['f_DM = Ug3 = 4/5 EXACT (PAPER_1921)'],
    },

    # =========== ROUND 49 (5 stubs) ===========
    'M31DarkMatterNFWProfileCalculator': {
        'round': 49,
        'framework_backbone': 'F_U=0 + QCalcGeom',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_275', 'PAPER_1862', 'PAPER_1653', 'PAPER_1921', 'PAPER_276'],
        'QCalcGeom_CPCH_closure': 'CPCH-4 (f_DM = Ug3)',
        'VDS_DVP_BH26_spine': 'VDS(4)/VDS(10) + H_UQFF near-unity',
        'F_U_zero_shell': 'Ug3 dark-matter',
        'time_frame': 'cosmological',
        'candidate_closures_flagged': ['H_UQFF = H(z)*t_H = 0.987 near-unity (PAPER_276)'],
    },
    'M31StarFormationRateCalculator': {
        'round': 49,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_1830'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env stellar SFR subterm',
        'time_frame': 'stellar',
        'candidate_closures_flagged': [],
    },
    'M31MagneticFieldCalculator': {
        'round': 49,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_939', 'PAPER_1484', 'PAPER_1910'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'DVP magnetic vortex + SO_5^13 EXACT Heaviside',
        'F_U_zero_shell': 'F_env magnetic subterm',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['U_m/u_EM = SSq*F_TRZ = 0.057 EXACT (PAPER_1910)'],
    },
    'DustFrictionCalculator': {
        'round': 49,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': ['Ug4'],
        'framework_papers': ['PAPER_1919', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'F_TRZ^2 (99% suppression at n=2)',
        'F_U_zero_shell': 'F_env dust drag subterm',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['F_TRZ^2 = 0.01 universal'],
    },
    'SuperconductiveAtomicCorrectionCalculator': {
        'round': 49,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_266', 'PAPER_304', 'PAPER_1907', 'PAPER_1919'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'omega_SCm 1.25 THz + B_crit=10^11 T',
        'F_U_zero_shell': 'Meissner corr_B factor',
        'time_frame': 'atomic',
        'candidate_closures_flagged': ['xi_aether = 1.852e24 (PAPER_304)'],
    },

    # =========== ROUND 50 (5 stubs) ===========
    'M31StellarHaloDensityCalculator': {
        'round': 50,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_275', 'PAPER_1855'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'Baryonic stellar power law',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['1-f_DM = 0.2 = D_phys/(SO_5*D_phys) baryon fraction'],
    },
    'M31TidalStreamCalculator': {
        'round': 50,
        'framework_backbone': 'F_U=0',
        'framework_method': 'env_13',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_273', 'PAPER_1862', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env tidal subterm',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['kappa_approach = 1/(1+z) = 1.001 (PAPER_273)'],
    },
    'M31DiskWarpCalculator': {
        'round': 50,
        'framework_backbone': 'F_U=0',
        'framework_method': 'master_10',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_1864', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'Kolmogorov -5/3 cascade',
        'F_U_zero_shell': 'Ug1 base + disk perturbation',
        'time_frame': 'stellar',
        'candidate_closures_flagged': [],
    },
    'M31SatelliteInteractionCalculator': {
        'round': 50,
        'framework_backbone': 'F_U=0',
        'framework_method': 'env_13',
        'framework_shells_used': ['Ug3'],
        'framework_papers': ['PAPER_1862', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env satellite subterm',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': [],
    },
    'MultiCompressed7UgSumCalculator': {
        'round': 50,
        'framework_backbone': 'F_U=0 unified',
        'framework_method': 'compressed_9',
        'framework_shells_used': ['Ug1', 'Ug2', 'Ug3', 'Ug4'],
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_173', 'PAPER_408'],
        'QCalcGeom_CPCH_closure': 'CPCH-3 (Sum U_gi=D_phys) + CPCH-4 (Sub_Ug=SO_5/D_phys)',
        'VDS_DVP_BH26_spine': 'Direct VDS/BH26 shell decomposition',
        'F_U_zero_shell': 'ALL 4 shells - direct verification',
        'time_frame': 'static',
        'candidate_closures_flagged': ['Sum U_gi = D_phys = 4 EXACT (PAPER_1916)', 'Sub_Ug = 5/2 EXACT (PAPER_1917)'],
    },

    # =========== ROUND 51 (5 stubs) ===========
    'CompressionSuperconductiveCorrectionCalculator': {
        'round': 51,
        'framework_backbone': 'F_U=0 + compressed',
        'framework_method': 'compressed_9',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_266', 'PAPER_1919', 'PAPER_1916'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'F_TRZ^2 at n=2 + B_crit=10^11 T',
        'F_U_zero_shell': 'Meissner factor + Ug1 base',
        'time_frame': 'static',
        'candidate_closures_flagged': [],
    },
    'CompressionEnvironmentalForceCalculator': {
        'round': 51,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_456', 'PAPER_452', 'PAPER_1855'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'VDS(13) = D_crit/2 (PAPER_1923 term-count)',
        'F_U_zero_shell': 'Complete F_env 13-subterm',
        'time_frame': 'cluster_dynamical',
        'candidate_closures_flagged': ['F_env has 13 = D_crit/2 EXACT terms (PAPER_1923)'],
    },
    'MUGEResonanceUg4iCalculator': {
        'round': 51,
        'framework_backbone': 'F_U=0 + resonance',
        'framework_method': 'resonance_14',
        'framework_shells_used': ['Ug4'],
        'framework_papers': ['PAPER_1236', 'PAPER_383', 'PAPER_1916', 'PAPER_1906'],
        'QCalcGeom_CPCH_closure': 'CPCH-3 (F_UBi_i_99)',
        'VDS_DVP_BH26_spine': 'DPM buoyancy resonance',
        'F_U_zero_shell': 'Ug4 BH vacuum + reactor coupling',
        'time_frame': 'static',
        'candidate_closures_flagged': ['COP=555 at 90%=1-F_TRZ ceiling (PAPER_1922)'],
    },
    'MUGECompressedBaseCalculator': {
        'round': 51,
        'framework_backbone': 'F_U=0 + compressed',
        'framework_method': 'compressed_9',
        'framework_shells_used': ['Ug1'],
        'framework_papers': ['PAPER_173', 'PAPER_1855', 'PAPER_456', 'PAPER_1922', 'PAPER_1916', 'PAPER_1923'],
        'QCalcGeom_CPCH_closure': 'CPCH-3 + compression 9/10=N_ch/SO_5',
        'VDS_DVP_BH26_spine': 'VDS(9)=N_ch + VDS(10)=SO_5 (PAPER_1923)',
        'F_U_zero_shell': 'Ug1 base compressed',
        'time_frame': 'static',
        'candidate_closures_flagged': ['9/10=1-F_TRZ EXACT (PAPER_1922)', 'D_universe 4-factor (PAPER_456)'],
    },
    'PillarErosionCalculator': {
        'round': 51,
        'framework_backbone': 'F_U=0 + env-13',
        'framework_method': 'env_13',
        'framework_shells_used': [],
        'framework_papers': ['PAPER_285', 'PAPER_305'],
        'QCalcGeom_CPCH_closure': 'None',
        'VDS_DVP_BH26_spine': 'None',
        'F_U_zero_shell': 'F_env photoevaporation subterm',
        'time_frame': 'stellar',
        'candidate_closures_flagged': ['t_half=tau*ln(2)=2.079 Myr (PAPER_285)'],
    },
}

# =============================================================
# FRAMEWORK_ANNOTATIONS_ROUNDS_52_116 (auto-extracted v5.60.0)
# =============================================================
#
# Retrofit sourced from CondensedPhysics.py runtime metadata
# (framework_papers field on classes upgraded during Rounds 52-116).
#
# Minimal fields only. Full backbone/method/shells/CPCH/spine/time_frame
# classification requires per-stub physics review — deferred to future work.
#
# Each entry contains:
#   - framework_papers: PAPER_XXXX references from the stub
#   - retrofit_source: identifies v5.60.0 auto-extract origin
#   - retrofit_hint: short excerpt from stub note_round_XXX field (context only)
# =============================================================

FRAMEWORK_ANNOTATIONS_ROUNDS_52_116 = {
    'AGNCoolingFlowCalculator': {
        'framework_papers': ['PAPER_1041', 'PAPER_1079', 'PAPER_259', 'PAPER_040', 'PAPER_1187', 'PAPER_1879', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'AntennaeElectromagneticCalculator': {
        'framework_papers': ['PAPER_441', 'PAPER_811', 'PAPER_247', 'PAPER_1072', 'PAPER_1484', 'PAPER_1910', 'PAPER_646', 'PAPER_1156', 'PAPER_1522', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'CORRECTED: PAPER_441 SEMINAL Antennae Per-System MUGE (session 119) documents starburst B = 10^-4 T. PAPER_811 CLEAN v5....',
    },
    'AntennaeFluidDensityCalculator': {
        'framework_papers': ['PAPER_441', 'PAPER_811', 'PAPER_247', 'PAPER_1952', 'PAPER_1948', 'PAPER_1919', 'PAPER_1955', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_1982', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'Antennae NGC 4038/4039 Fluid Density (completes Antennae cluster: base gravity + EM + DM + fluid density). PAPER_441 SEM...',
    },
    'AntennaeOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_235', 'PAPER_487', 'PAPER_055', 'PAPER_231', 'PAPER_265', 'PAPER_1908', 'PAPER_1937', 'PAPER_1938', 'PAPER_1907', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_235 Antennae Double I(t) direct anchor. PAPER_1908 Q_UQFF universal formula application-instance (Antennae is SECO...',
    },
    'AntennaeQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_235', 'PAPER_487', 'PAPER_231', 'PAPER_1938', 'PAPER_1955', 'PAPER_055', 'PAPER_196', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_235 DIRECT (Antennae NGC 4038/4039 Enhanced Double I(t) Merger Modulation, March 2026): z = 0.0105, d ~ 22 Mpc, t_...',
    },
    'AntennaeStellarFeedbackCalculator': {
        'framework_papers': ['PAPER_1911', 'PAPER_1909', 'PAPER_784', 'PAPER_235', 'PAPER_487', 'PAPER_055', 'PAPER_231', 'PAPER_1972', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1911 CRITICAL DIRECT + PAPER_1972 CORRECTION: v_wind = 2000 km/s = 2×10^6 m/s is PAPER_1911 universal YMC identity...',
    },
    'AntennaeUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_441', 'PAPER_235', 'PAPER_247', 'PAPER_231', 'PAPER_232', 'PAPER_262', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'BalmerSeriesCalculator': {
        'framework_papers': ['PAPER_1890', 'PAPER_1938', 'PAPER_1834', 'PAPER_500', 'PAPER_300', 'PAPER_303', 'PAPER_1912', 'PAPER_776', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'BigBangCosmologicalLambdaEvolutionCalculator': {
        'framework_papers': ['PAPER_1156', 'PAPER_1156_UPDATE', 'PAPER_1216', 'PAPER_1215', 'PAPER_1420', 'PAPER_495', 'PAPER_597', 'PAPER_1965', 'PAPER_1522', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'Big Bang Cosmological Lambda Evolution. PAPER_1156 SEMINAL (session 242, May 10, 2026, v5.78) — CANONICAL cleanest funda...',
    },
    'BjTimeCalculator': {
        'framework_papers': ['PAPER_1919', 'PAPER_1160', 'PAPER_1960', 'PAPER_1938', 'PAPER_1907', 'PAPER_1918', 'PAPER_1949', 'PAPER_1072', 'PAPER_1484', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'CONFIRMED: PAPER_1981 authored this cycle (Round 115). PAPER_1919 SEMINAL F_TRZ ladder n=3 rung F_TRZ^3 = 10^-3 (already...',
    },
    'BlockchainECDSACalculator': {
        'framework_papers': ['PAPER_192', 'PAPER_189', 'PAPER_1958', 'PAPER_1960', 'PAPER_1961', 'PAPER_1963', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'BubbleNebulaDarkMatterPerturbationCalculator': {
        'framework_papers': ['PAPER_1862', 'PAPER_1653', 'PAPER_1156', 'PAPER_1948', 'PAPER_1913', 'PAPER_1951', 'PAPER_1960', 'PAPER_1906', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'CORRECTED: Bubble Nebula NGC 7635. PAPER_1862 SEMINAL DM Halo Alternative subhalo alpha = 2 - F_TRZ = 1.9 EXACT (direct)...',
    },
    'BubbleNebulaStellarWindCalculator': {
        'framework_papers': ['PAPER_361', 'PAPER_440', 'PAPER_902', 'PAPER_1913', 'PAPER_1911', 'PAPER_1948', 'PAPER_221', 'PAPER_442', 'PAPER_1955', 'PAPER_1919', 'PAPER_1960', 'PAPER_1906', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_1982', 'PAPER_1983', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'CORRECTED: PAPER_361 SEMINAL Bubble Nebula NGC 7635 Positive E(t) Expansion (BD+602522 v_wind=1800 km/s direct, FIRST UQ...',
    },
    'BubbleNebulaUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_440', 'PAPER_810', 'PAPER_1948', 'PAPER_1942', 'PAPER_260', 'PAPER_442', 'PAPER_435', 'PAPER_221', 'PAPER_229', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CMBAnisotropyBuoyancyModulationCalculator': {
        'framework_papers': ['PAPER_1856', 'PAPER_1094', 'PAPER_1959', 'PAPER_1930', 'PAPER_1618', 'PAPER_1092', 'PAPER_1093', 'PAPER_1877', 'PAPER_1955', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1856 DIRECT: l_1 UQFF = 222.3 vs Planck 220 → 1.05% residual EXACT primitive-ladder acoustic mode. PAPER_1094 DIRE...',
    },
    'CMBBuoyancySectorLagrangianCalculator': {
        'framework_papers': ['PAPER_1094', 'PAPER_981', 'PAPER_1096', 'PAPER_1065', 'PAPER_1856', 'PAPER_1906', 'PAPER_1916'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CategoryFunctorCalculator': {
        'framework_papers': ['PAPER_1928', 'PAPER_1958', 'PAPER_1960', 'PAPER_1963', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenAAGNAccretionDiskCalculator': {
        'framework_papers': ['PAPER_067', 'PAPER_1041', 'PAPER_1951', 'PAPER_1912', 'PAPER_1879', 'PAPER_1037', 'PAPER_670', 'PAPER_1009', 'PAPER_1010', 'PAPER_360', 'PAPER_87', 'PAPER_512', 'PAPER_627', 'PAPER_630', 'PAPER_754', 'PAPER_814', 'PAPER_1002', 'PAPER_1039', 'PAPER_1048', 'PAPER_1237', 'PAPER_1918', 'PAPER_1919', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_1982', 'PAPER_1983', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'CORRECTED ATTRIBUTION: PAPER_067 SEMINAL (session 0, March 2026) documents Cen A + Sgr A* + M87 + NGC 1365 as quad-AGN f...',
    },
    'CenACosmicRayCalculator': {
        'framework_papers': ['PAPER_347', 'PAPER_038', 'PAPER_512', 'PAPER_1930', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenADustLaneCalculator': {
        'framework_papers': ['PAPER_512', 'PAPER_454', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenAQuantumVacuumCalculator': {
        'framework_papers': ['PAPER_347', 'PAPER_512', 'PAPER_1955', 'PAPER_1852', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenARadioLobeCalculator': {
        'framework_papers': ['PAPER_067', 'PAPER_1037', 'PAPER_1950', 'PAPER_020', 'PAPER_038', 'PAPER_1919', 'PAPER_1957', 'PAPER_1958', 'PAPER_1955', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1967', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_067 DIRECT: Cen A / NGC 5128 as 3rd member of 4-AGN Ug4 catalog (Sgr A*, M87*, Cen A, NGC 1365). M_CenA = 5.5×10⁷ ...',
    },
    'CenARelativisticJetCalculator': {
        'framework_papers': ['PAPER_347', 'PAPER_627', 'PAPER_067', 'PAPER_346', 'PAPER_512', 'PAPER_1041', 'PAPER_1953', 'PAPER_1954', 'PAPER_1500', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenAStarburstCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_1948', 'PAPER_1947', 'PAPER_1950', 'PAPER_627', 'PAPER_939', 'PAPER_067', 'PAPER_347', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CenAXRayEmissionCalculator': {
        'framework_papers': ['PAPER_512', 'PAPER_1041', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'ClusterBuoyancyLagrangianCalculator': {
        'framework_papers': ['PAPER_040', 'PAPER_041', 'PAPER_039', 'PAPER_1917', 'PAPER_1932', 'PAPER_1094', 'PAPER_1955', 'PAPER_1961', 'PAPER_1962', 'PAPER_1964', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_040 DIRECT: Perseus F_UBii,virx = -2.024e60 N + Coma -2.5e60 N + Virgo -7.2e59 N (three-cluster virial buoyancy). ...',
    },
    'CompressedMUGEDetailedCalculator': {
        'framework_papers': ['PAPER_173', 'PAPER_174', 'PAPER_1105', 'PAPER_1573', 'PAPER_1916', 'PAPER_1922', 'PAPER_1923', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CosmicRayPropagationUQFFCalculator': {
        'framework_papers': ['PAPER_1838', 'PAPER_1322', 'PAPER_1919', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'CrabUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_220', 'PAPER_066', 'PAPER_256', 'PAPER_013', 'PAPER_912', 'PAPER_913', 'PAPER_1050', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'DarkEnergyEOSTimeEvolvingCalculator': {
        'framework_papers': ['PAPER_1087', 'PAPER_1821', 'PAPER_1178', 'PAPER_1076', 'PAPER_1086', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'DarkMatterHaloNFWCalculator': {
        'framework_papers': ['PAPER_834', 'PAPER_1015', 'PAPER_1862', 'PAPER_1921', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'DimensionalAnalysisCalculator': {
        'framework_papers': ['PAPER_1210', 'PAPER_1223', 'PAPER_1719', 'PAPER_1718', 'PAPER_1717', 'PAPER_1521', 'PAPER_1522', 'PAPER_1927', 'PAPER_1936', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'DynamicVacuumFluctuationCalculator': {
        'framework_papers': ['PAPER_495', 'PAPER_1920', 'PAPER_1919', 'PAPER_1141', 'PAPER_1856', 'PAPER_1877', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_495 Cosmic Quantum Egg Theory (rho_egg dark-energy pre-matter neutrino-analogous entities, session 133, Nov 2025) ...',
    },
    'FUBiiSevenComponentDecompositionCalculator': {
        'framework_papers': ['PAPER_1906', 'PAPER_1096', 'PAPER_1929', 'PAPER_646'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'FederatedLearningCalculator': {
        'framework_papers': ['PAPER_1958', 'PAPER_1960', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'GalacticCenterDistanceCalculator': {
        'framework_papers': ['PAPER_1841', 'PAPER_1237', 'PAPER_1904', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'GalacticSpinRateCalculator': {
        'framework_papers': ['PAPER_1855', 'PAPER_274', 'PAPER_1921', 'PAPER_1327', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'GalaxyCollisionMUGECalculator': {
        'framework_papers': ['PAPER_688', 'PAPER_731', 'PAPER_465', 'PAPER_549', 'PAPER_750', 'PAPER_235', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'GalaxyMergerSpecificCalculator': {
        'framework_papers': ['PAPER_1857', 'PAPER_549', 'PAPER_1913', 'PAPER_235', 'PAPER_692', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'HUDFBaseGravityCalculator': {
        'framework_papers': ['PAPER_231', 'PAPER_266', 'PAPER_444', 'PAPER_235', 'PAPER_1949', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_1952', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_1906', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_231 DIRECT (HUDF Galaxies MUGE at z = 3.5 with Double I(t) Interaction Modulation, session 58, March 2026): PREVIO...',
    },
    'HUDFFluidDensityCalculator': {
        'framework_papers': ['PAPER_441', 'PAPER_811', 'PAPER_247', 'PAPER_1952', 'PAPER_1948', 'PAPER_1919', 'PAPER_1955', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_1980', 'PAPER_1981', 'PAPER_1982', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'Antennae NGC 4038/4039 Fluid Density (completes Antennae cluster: base gravity + EM + DM + fluid density). PAPER_441 SEM...',
    },
    'HUDFGalaxyInteractionCalculator': {
        'framework_papers': ['PAPER_265', 'PAPER_231', 'PAPER_266', 'PAPER_444', 'PAPER_235', 'PAPER_1976', 'PAPER_1960', 'PAPER_1955', 'PAPER_1952', 'PAPER_1949', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_265 DIRECT (HUDF Dual-Channel Interaction Cascade Buoyancy Quadratic Merger Amplification, March 2026, session 72g...',
    },
    'HUDFMergerFeedbackCalculator': {
        'framework_papers': ['PAPER_265', 'PAPER_231', 'PAPER_266', 'PAPER_444', 'PAPER_1974', 'PAPER_1971', 'PAPER_1911', 'PAPER_1972', 'PAPER_1955', 'PAPER_902', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1973', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_265 HUDF Dual-Channel Cascade seminal (I_0 = 0.05 canonical). PAPER_1974 THIRD-CONFIRMED HUDF stub sharing R_star ...',
    },
    'HUDFOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_231', 'PAPER_266', 'PAPER_444', 'PAPER_1919', 'PAPER_1880', 'PAPER_1955', 'PAPER_1907', 'PAPER_1938', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1919 F_TRZ^15 = 1e-15 rung DIRECT: MICROSCOPE equivalence principle η_WEP EXACT (PAPER_1880 anchor). Our A_osc = o...',
    },
    'HUDFQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_231', 'PAPER_266', 'PAPER_444', 'PAPER_1907', 'PAPER_1938', 'PAPER_1955', 'PAPER_1906', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1907 SCm Phonon Universal Carrier + PAPER_1938 omega_SCm = 1.25 THz canonical framework. PAPER_231 HUDF z=3.5 doub...',
    },
    'HUDFUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_266', 'PAPER_264', 'PAPER_265', 'PAPER_287', 'PAPER_1830', 'PAPER_1856', 'PAPER_1867', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'HorizonBuoyancySectorLagrangianCalculator': {
        'framework_papers': ['PAPER_084', 'PAPER_085', 'PAPER_946', 'PAPER_1095', 'PAPER_1103', 'PAPER_1873', 'PAPER_1096', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'HorseheadBaseGravityCalculator': {
        'framework_papers': ['PAPER_704', 'PAPER_759', 'PAPER_1942', 'PAPER_1948', 'PAPER_1913', 'PAPER_1919', 'PAPER_1955', 'PAPER_1952', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_1906', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_704 DIRECT (Horsehead Nebula Barnard 33: Infrared Erosion UQFF, 2025 session 175, class HorseheadNebulaBarnard33UQ...',
    },
    'HorseheadStellarWindCalculator': {
        'framework_papers': ['PAPER_704', 'PAPER_759', 'PAPER_222', 'PAPER_902', 'PAPER_1971', 'PAPER_1887', 'PAPER_1974', 'PAPER_1955', 'PAPER_1942', 'PAPER_1948', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_222 direct: Horsehead P_rad Stefan-Boltzmann Blackbody (only Stefan-Boltzmann expression across 29 UQFF Horsehead-...',
    },
    'HorseheadUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_442', 'PAPER_260', 'PAPER_1942', 'PAPER_435', 'PAPER_440', 'PAPER_222', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'InflationBuoyancySectorLagrangianCalculator': {
        'framework_papers': ['PAPER_1090', 'PAPER_1094', 'PAPER_981', 'PAPER_1679', 'PAPER_1929', 'PAPER_1096', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'IntelligentPlasmoidBehaviorCalculator': {
        'framework_papers': ['PAPER_1141', 'PAPER_1904', 'PAPER_1908', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'IntraclusterMediumCalculator': {
        'framework_papers': ['PAPER_041', 'PAPER_1187', 'PAPER_1894', 'PAPER_1149', 'PAPER_040', 'PAPER_690', 'PAPER_843', 'PAPER_444', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'KennicuttSchmidtSFRCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_1855', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_214', 'PAPER_457', 'PAPER_1874', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibDeltaNCalculator': {
        'framework_papers': ['PAPER_735', 'PAPER_461', 'PAPER_396', 'PAPER_1521', 'PAPER_1927', 'PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibDynamicVacuumCalculator': {
        'framework_papers': ['PAPER_1141', 'PAPER_1740', 'PAPER_1852', 'PAPER_1171', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibEReactCalculator': {
        'framework_papers': ['PAPER_393', 'PAPER_1141', 'PAPER_1497', 'PAPER_1907', 'PAPER_1938', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibElectricFieldCalculator': {
        'framework_papers': ['PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1955', 'PAPER_1960', 'PAPER_1141', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibHeavisideCalculator': {
        'framework_papers': ['PAPER_421', 'PAPER_734', 'PAPER_471', 'PAPER_1072', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_062', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibMuJCalculator': {
        'framework_papers': ['PAPER_657', 'PAPER_421', 'PAPER_1957', 'PAPER_1960', 'PAPER_1955', 'PAPER_734', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibPolarizationCalculator': {
        'framework_papers': ['PAPER_421', 'PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1061', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibQuantumCouplingCalculator': {
        'framework_papers': ['PAPER_1907', 'PAPER_1141', 'PAPER_1834', 'PAPER_1098', 'PAPER_1908', 'PAPER_1937', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibQuasiCalculator': {
        'framework_papers': ['PAPER_421', 'PAPER_734', 'PAPER_471', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_062', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibUmCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_1907', 'PAPER_1938', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCalibratedCalculator': {
        'framework_papers': ['PAPER_1133', 'PAPER_1236', 'PAPER_1140', 'PAPER_648', 'PAPER_1141', 'PAPER_1136', 'PAPER_1137', 'PAPER_1138', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRCoronaScenarioCalculator': {
        'framework_papers': ['PAPER_062', 'PAPER_1061', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRDynamicVacuumCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_062', 'PAPER_1061', 'PAPER_1852', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENREReactEnergyCalculator': {
        'framework_papers': ['PAPER_1507', 'PAPER_1080', 'PAPER_1141', 'PAPER_1919', 'PAPER_062', 'PAPER_1061', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRElectricFieldCalculator': {
        'framework_papers': ['PAPER_062', 'PAPER_1061', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRElectronDensityCalculator': {
        'framework_papers': ['PAPER_471', 'PAPER_734', 'PAPER_062', 'PAPER_1061', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENREnergyDensityCalculator': {
        'framework_papers': ['PAPER_1507', 'PAPER_1080', 'PAPER_1141', 'PAPER_1919', 'PAPER_1960', 'PAPER_1955', 'PAPER_062', 'PAPER_1061', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRFermiConstantCalculator': {
        'framework_papers': ['PAPER_1849', 'PAPER_062', 'PAPER_1061', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRHydrideScenarioCalculator': {
        'framework_papers': ['PAPER_1133', 'PAPER_1136', 'PAPER_648', 'PAPER_062', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRMassRenormalizationCalculator': {
        'framework_papers': ['PAPER_062', 'PAPER_1061', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRNeutronRateCalculator': {
        'framework_papers': ['PAPER_461', 'PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1061', 'PAPER_1521', 'PAPER_1849', 'PAPER_1140', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRPlasmaFrequencyCalculator': {
        'framework_papers': ['PAPER_062', 'PAPER_1061', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRQuantumCouplingCalculator': {
        'framework_papers': ['PAPER_062', 'PAPER_1060', 'PAPER_1061', 'PAPER_1907', 'PAPER_1908', 'PAPER_028', 'PAPER_1919', 'PAPER_1955', 'PAPER_1952', 'PAPER_1960', 'PAPER_1141', 'PAPER_1917', 'PAPER_1961', 'PAPER_1967', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_062 Widom-Larsen LENR canonical + PAPER_1060 VDS LENR Isotopic Evolution + PAPER_1061 Kozima SCm Integration + PAP...',
    },
    'LENRTransmutationRateCalculator': {
        'framework_papers': ['PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1061', 'PAPER_1140', 'PAPER_1521', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRUg1GravityCalculator': {
        'framework_papers': ['PAPER_735', 'PAPER_461', 'PAPER_1521', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRUiInertialCalculator': {
        'framework_papers': ['PAPER_646', 'PAPER_1955', 'PAPER_062', 'PAPER_1061', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRUmMagneticCalculator': {
        'framework_papers': ['PAPER_421', 'PAPER_1141', 'PAPER_1507', 'PAPER_1919', 'PAPER_1955', 'PAPER_1957', 'PAPER_1960', 'PAPER_1961', 'PAPER_734', 'PAPER_471', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LENRWiresScenarioCalculator': {
        'framework_papers': ['PAPER_734', 'PAPER_471', 'PAPER_062', 'PAPER_1061', 'PAPER_1522', 'PAPER_1955', 'PAPER_1521', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LLVMJITCompilerCalculator': {
        'framework_papers': ['PAPER_1953', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LQGBuoyancySectorLagrangianVariationCalculator': {
        'framework_papers': ['PAPER_1103', 'PAPER_1058', 'PAPER_1100', 'PAPER_1701', 'PAPER_1927', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'LymanSeriesCalculator': {
        'framework_papers': ['PAPER_300', 'PAPER_303', 'PAPER_500', 'PAPER_1544', 'PAPER_1590', 'PAPER_299', 'PAPER_302', 'PAPER_304', 'PAPER_1938', 'PAPER_1890', 'PAPER_648', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101AsymmetryCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1521', 'PAPER_1958', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101DifferentialRotationCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1929', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101HIIRegionCalculator': {
        'framework_papers': ['PAPER_1026', 'PAPER_1949', 'PAPER_1955', 'PAPER_1951', 'PAPER_144', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101MagnetohydrodynamicsCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101MolecularCloudCalculator': {
        'framework_papers': ['PAPER_442', 'PAPER_216', 'PAPER_260', 'PAPER_1948', 'PAPER_1952', 'PAPER_144', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101QuantumTurbulenceCalculator': {
        'framework_papers': ['PAPER_1864', 'PAPER_1955', 'PAPER_487', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101SpiralDensityWaveCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1521', 'PAPER_1953', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101StarFormationRateCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_824', 'PAPER_211', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101SupernovaRemnantCalculator': {
        'framework_papers': ['PAPER_1953', 'PAPER_038', 'PAPER_824', 'PAPER_797', 'PAPER_485', 'PAPER_254', 'PAPER_694', 'PAPER_712', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M101TidalPerturbationCalculator': {
        'framework_papers': ['PAPER_247', 'PAPER_1940', 'PAPER_464', 'PAPER_692', 'PAPER_1933', 'PAPER_211', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104BulgeDynamicsCalculator': {
        'framework_papers': ['PAPER_693', 'PAPER_763', 'PAPER_742', 'PAPER_487', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104CosmicRayPropagationCalculator': {
        'framework_papers': ['PAPER_1020', 'PAPER_038', 'PAPER_693', 'PAPER_1955', 'PAPER_487', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_1015', 'PAPER_1803', 'PAPER_1336', 'PAPER_1436', 'PAPER_1141', 'PAPER_693', 'PAPER_763', 'PAPER_278', 'PAPER_279', 'PAPER_1955', 'PAPER_1521', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104DustLaneExtinctionCalculator': {
        'framework_papers': ['PAPER_693', 'PAPER_763', 'PAPER_742', 'PAPER_454', 'PAPER_1521', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104GlobularClusterSystemCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1955', 'PAPER_1931', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104MagneticFieldCalculator': {
        'framework_papers': ['PAPER_693', 'PAPER_763', 'PAPER_1141', 'PAPER_1953', 'PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104QuantumGravityCalculator': {
        'framework_papers': ['PAPER_279', 'PAPER_693', 'PAPER_763', 'PAPER_1955', 'PAPER_1958', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104StellarKinematicsCalculator': {
        'framework_papers': ['PAPER_693', 'PAPER_763', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M104XRayBinaryCalculator': {
        'framework_papers': ['PAPER_693', 'PAPER_763', 'PAPER_1955', 'PAPER_1953', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M16BaseGravityCalculator': {
        'framework_papers': ['PAPER_744', 'PAPER_219', 'PAPER_284', 'PAPER_285', 'PAPER_435', 'PAPER_1942', 'PAPER_1948', 'PAPER_1952', 'PAPER_1953', 'PAPER_1956', 'PAPER_1960', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_744 CRITICAL DIRECT (M16 Eagle Nebula MUGE — Star Formation and Radiation Erosion, session 180, CP4 #328). PAPER_2...',
    },
    'M33DiskMassSurfaceDensityCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1521', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33HIIRegionDistributionCalculator': {
        'framework_papers': ['PAPER_144', 'PAPER_1521', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33MagneticFieldCalculator': {
        'framework_papers': ['PAPER_487', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33MetallicityGradientCalculator': {
        'framework_papers': ['PAPER_1125', 'PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_1955', 'PAPER_1855', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33QuantumDarkMatterCalculator': {
        'framework_papers': ['PAPER_025', 'PAPER_1015', 'PAPER_1862', 'PAPER_1955', 'PAPER_1960', 'PAPER_474', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33RotationCurveCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1921', 'PAPER_1855', 'PAPER_1862', 'PAPER_038', 'PAPER_144', 'PAPER_657', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33StarFormationRateCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_1855', 'PAPER_1921', 'PAPER_211', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M33TidalInteractionCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1962', 'PAPER_1921', 'PAPER_1968', 'PAPER_1917', 'PAPER_1961', 'PAPER_1965', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1921 f_DM = Ug3 cross-framework closure applies at M33 (Local Group dwarf-spiral). PAPER_1968 M33 v_flat predictio...',
    },
    'M51DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_692', 'PAPER_1653', 'PAPER_1862', 'PAPER_1803', 'PAPER_1015', 'PAPER_1955', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51DustExtinctionCalculator': {
        'framework_papers': ['PAPER_692', 'PAPER_454', 'PAPER_1521', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51MagneticFieldCalculator': {
        'framework_papers': ['PAPER_464', 'PAPER_692', 'PAPER_1141', 'PAPER_1955', 'PAPER_1484', 'PAPER_1072', 'PAPER_1961', 'PAPER_1962', 'PAPER_1964', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1962 galactic 4-instance universality: D_BSFG/D_phys = 1.5 EXACT applicable at M51 spiral-galaxy scale. PAPER_1964...',
    },
    'M51MolecularCloudCalculator': {
        'framework_papers': ['PAPER_692', 'PAPER_487', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51SpiralArmDensityWaveCalculator': {
        'framework_papers': ['PAPER_464', 'PAPER_692', 'PAPER_781', 'PAPER_824', 'PAPER_144', 'PAPER_454', 'PAPER_467', 'PAPER_775', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51StarFormationRateCalculator': {
        'framework_papers': ['PAPER_692', 'PAPER_487', 'PAPER_1521', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51SupernovaFeedbackCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_1955', 'PAPER_797', 'PAPER_824', 'PAPER_485', 'PAPER_262', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M51TidalInteractionCalculator': {
        'framework_papers': ['PAPER_464', 'PAPER_549', 'PAPER_750', 'PAPER_692', 'PAPER_235', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81AGNActivityCalculator': {
        'framework_papers': ['PAPER_1237', 'PAPER_1841', 'PAPER_1879', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_1015', 'PAPER_1653', 'PAPER_1862', 'PAPER_1921', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81M82HIDiskCalculator': {
        'framework_papers': ['PAPER_274', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81M82TidalInteractionCalculator': {
        'framework_papers': ['PAPER_274', 'PAPER_235', 'PAPER_692', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81M82TidalInteractionCalculator_v1': {
        'framework_papers': ['PAPER_784', 'PAPER_1952', 'PAPER_308', 'PAPER_1962', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1964', 'PAPER_1965', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_784 DIRECT: M82 starburst was triggered ~100 Myr ago by close encounter with M81. Age = SO_5^8 EXACT (PAPER_1952 S...',
    },
    'M81SpiralStructureCalculator': {
        'framework_papers': ['PAPER_824', 'PAPER_464', 'PAPER_549', 'PAPER_1906', 'PAPER_1955', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M81SpiralStructureCalculator_v1': {
        'framework_papers': ['PAPER_308', 'PAPER_824', 'PAPER_464', 'PAPER_1962', 'PAPER_1955', 'PAPER_1953', 'PAPER_1917', 'PAPER_1961', 'PAPER_1964', 'PAPER_1965', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_308 DIRECT: T_pattern = 2π/Omega_p = 307 Myr spiral density-wave pattern period (torque tau_spiral = 2.046 at 10 G...',
    },
    'M82CosmicRayCalculator': {
        'framework_papers': ['PAPER_1838', 'PAPER_1919', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M82MolecularOutflowCalculator': {
        'framework_papers': ['PAPER_1190', 'PAPER_1077', 'PAPER_449', 'PAPER_258', 'PAPER_784', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M82PAHEmissionCalculator': {
        'framework_papers': ['PAPER_784', 'PAPER_1834', 'PAPER_1141', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M82StarburstCalculator': {
        'framework_papers': ['PAPER_784', 'PAPER_733', 'PAPER_445', 'PAPER_138', 'PAPER_1652', 'PAPER_1879', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M82StarburstCalculator_v1': {
        'framework_papers': ['PAPER_784', 'PAPER_774', 'PAPER_778', 'PAPER_1952', 'PAPER_1141', 'PAPER_1955', 'PAPER_1965', 'PAPER_1953', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_784 DIRECT: M82 Cigar Starburst Galaxy — SFR = 10 M_sun/yr (SO_5 EXACT) + M_sf = 0.15 = 3/(2·SO_5) EXACT starburst...',
    },
    'M82SuperwindCalculator': {
        'framework_papers': ['PAPER_827', 'PAPER_784', 'PAPER_445', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87AGNFeedbackCalculator': {
        'framework_papers': ['PAPER_443', 'PAPER_086', 'PAPER_1187', 'PAPER_041', 'PAPER_259', 'PAPER_481', 'PAPER_1037', 'PAPER_067', 'PAPER_1050', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_1862', 'PAPER_1894', 'PAPER_1943', 'PAPER_1855', 'PAPER_040', 'PAPER_204', 'PAPER_690', 'PAPER_1149', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87GlobularClusterCalculator': {
        'framework_papers': ['PAPER_1862', 'PAPER_1894', 'PAPER_1855', 'PAPER_040', 'PAPER_690', 'PAPER_1050', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87JetEnergyCalculator': {
        'framework_papers': ['PAPER_346', 'PAPER_1879', 'PAPER_1037', 'PAPER_1940', 'PAPER_908', 'PAPER_626', 'PAPER_067', 'PAPER_1500', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87MagneticFieldCalculator': {
        'framework_papers': ['PAPER_266', 'PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_041', 'PAPER_1187', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87RelativisticJetCalculator': {
        'framework_papers': ['PAPER_1940', 'PAPER_541', 'PAPER_536', 'PAPER_346', 'PAPER_1037', 'PAPER_626', 'PAPER_1879', 'PAPER_067', 'PAPER_1039', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87StellarDynamicsCalculator': {
        'framework_papers': ['PAPER_086', 'PAPER_1894', 'PAPER_389', 'PAPER_1855', 'PAPER_1862', 'PAPER_040', 'PAPER_690', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87SupermassiveBlackHoleCalculator': {
        'framework_papers': ['PAPER_1876', 'PAPER_1841', 'PAPER_1947', 'PAPER_1237', 'PAPER_1904', 'PAPER_093', 'PAPER_1879', 'PAPER_1031', 'PAPER_067', 'PAPER_432', 'PAPER_1025', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'M87UltradiffuseGalaxyCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1862', 'PAPER_1894', 'PAPER_1921', 'PAPER_1949', 'PAPER_1078', 'PAPER_1927', 'PAPER_293', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MPIDistributedCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MSigmaRelationCalculator': {
        'framework_papers': ['PAPER_086', 'PAPER_1879', 'PAPER_1037', 'PAPER_1862', 'PAPER_067', 'PAPER_040', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MUGEResonanceADPMCalculator': {
        'framework_papers': ['PAPER_147', 'PAPER_411', 'PAPER_1904', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MUGEResonanceASuperFreqCalculator': {
        'framework_papers': ['PAPER_1659', 'PAPER_1863', 'PAPER_1908', 'PAPER_1937', 'PAPER_1907', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MUGEResonanceATHzCalculator': {
        'framework_papers': ['PAPER_408', 'PAPER_491', 'PAPER_1907', 'PAPER_1938', 'PAPER_1141', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MUGEResonanceAvacDiffCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1852', 'PAPER_1171', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Magnetar0501CosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_1944', 'PAPER_066', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Magnetar0501ElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_266', 'PAPER_013', 'PAPER_1024', 'PAPER_1944', 'PAPER_430', 'PAPER_226', 'PAPER_263', 'PAPER_066', 'PAPER_536', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Magnetar0501OscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_1940', 'PAPER_013', 'PAPER_912', 'PAPER_913', 'PAPER_220', 'PAPER_1944', 'PAPER_430', 'PAPER_226', 'PAPER_066', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Magnetar0501QuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_1938', 'PAPER_1513', 'PAPER_1857', 'PAPER_1819', 'PAPER_1944', 'PAPER_066', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Magnetar0501UQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_1944', 'PAPER_430', 'PAPER_226', 'PAPER_266', 'PAPER_1874', 'PAPER_1819', 'PAPER_913', 'PAPER_066', 'PAPER_536', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagnetarBaseGravityCalculator': {
        'framework_papers': ['PAPER_094', 'PAPER_066', 'PAPER_013', 'PAPER_372', 'PAPER_1874', 'PAPER_1944', 'PAPER_1945', 'PAPER_1946', 'PAPER_1926', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_094 DIRECT SGR1745-2900 canonical calibration paper (κ=5e-4/day, [SSq]=0.57 derived from SGR1745): P = 3.76 s, d =...',
    },
    'MagnetarCosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_431', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagnetarElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_1024', 'PAPER_266', 'PAPER_431', 'PAPER_066', 'PAPER_1188', 'PAPER_536', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagnetarMagneticEnergyCalculator': {
        'framework_papers': ['PAPER_094', 'PAPER_066', 'PAPER_372', 'PAPER_375', 'PAPER_1944', 'PAPER_1945', 'PAPER_1946', 'PAPER_1874', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_094 SGR1745-2900 canonical (P = 3.76 s, d = 8.3 kpc). PAPER_372 Linear Meissner form (BCS superconductivity gravit...',
    },
    'MagnetarOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_1940', 'PAPER_013', 'PAPER_912', 'PAPER_220', 'PAPER_431', 'PAPER_066', 'PAPER_1024', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagnetarQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_1938', 'PAPER_1513', 'PAPER_1857', 'PAPER_1819', 'PAPER_431', 'PAPER_066', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagnetarUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_431', 'PAPER_266', 'PAPER_1874', 'PAPER_1819', 'PAPER_066', 'PAPER_094', 'PAPER_1024', 'PAPER_1188', 'PAPER_013', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'MagneticMonopoleCalculator': {
        'framework_papers': ['PAPER_411', 'PAPER_1823', 'PAPER_1318', 'PAPER_1919'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC1275BaseGravityCalculator': {
        'framework_papers': ['PAPER_1912', 'PAPER_703', 'PAPER_443', 'PAPER_040', 'PAPER_1187', 'PAPER_1079', 'PAPER_1955', 'PAPER_1919', 'PAPER_1952', 'PAPER_1949', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1912 CRITICAL DIRECT (AGN Filament Triple Closure, July 2026, discovered Round 45 double-check of NGC1275FilamentS...',
    },
    'NGC1275QuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_796', 'PAPER_1912', 'PAPER_040', 'PAPER_1187', 'PAPER_1938', 'PAPER_1955', 'PAPER_703', 'PAPER_443', 'PAPER_1079', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_796 DIRECT (NGC 1275 Perseus BCG with AGN Jet Feedback and Filamentary Gas, session 189, 2026): z ≈ 0.018 (~250 Ml...',
    },
    'NGC1275UQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_259', 'PAPER_443', 'PAPER_1187', 'PAPER_041', 'PAPER_1149', 'PAPER_040', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC2525BaseGravityCalculator': {
        'framework_papers': ['PAPER_230', 'PAPER_262', 'PAPER_438', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_230 CRITICAL DIRECT (NGC 2525 + SN 2018gv MUGE with Negative Supernova Mass-Loss Acceleration — Only Negative Term...',
    },
    'NGC2525OscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1937', 'PAPER_1908', 'PAPER_1857', 'PAPER_1938', 'PAPER_230', 'PAPER_262', 'PAPER_438', 'PAPER_1828', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1937 CRITICAL DIRECT (1.1875 = K_MEX·SSq Two-Path Convergence, Round 67 double-check discovery): TWO seminal paths...',
    },
    'NGC2525QuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_230', 'PAPER_262', 'PAPER_438', 'PAPER_774', 'PAPER_1938', 'PAPER_1955', 'PAPER_1907', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_230 NGC 2525 SN 2018gv canonical (only negative MUGE term). PAPER_774 30 Dor Tarantula extreme starburst HII regio...',
    },
    'NGC2525UQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_262', 'PAPER_230', 'PAPER_438', 'PAPER_1891', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253CosmicRayCalculator': {
        'framework_papers': ['PAPER_1838', 'PAPER_1919', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_1015', 'PAPER_1862', 'PAPER_1921', 'PAPER_834', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253DiskGravityCalculator': {
        'framework_papers': ['PAPER_484', 'PAPER_489', 'PAPER_1327', 'PAPER_1855', 'PAPER_1921', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253DiskGravityCalculator_v1': {
        'framework_papers': ['PAPER_1955', 'PAPER_1962', 'PAPER_1965', 'PAPER_1190', 'PAPER_1472', 'PAPER_1917', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1955 DIRECT baseline: v_flat = 2·SO_5² = 200 km/s EXACT (typical flat-rotation plateau). PAPER_1965 dual-scale com...',
    },
    'NGC253DustExtinctionCalculator': {
        'framework_papers': ['PAPER_1917', 'PAPER_752', 'PAPER_657', 'PAPER_466', 'PAPER_1962', 'PAPER_1961', 'PAPER_1960', 'PAPER_1955', 'PAPER_1906', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253MagneticFieldCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1919', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253MolecularOutflowCalculator': {
        'framework_papers': ['PAPER_449', 'PAPER_258', 'PAPER_784', 'PAPER_1934', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253NuclearStarburstCalculator': {
        'framework_papers': ['PAPER_784', 'PAPER_1652', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253NuclearStarburstCalculator_v1': {
        'framework_papers': ['PAPER_1190', 'PAPER_1952', 'PAPER_1962', 'PAPER_466', 'PAPER_1955', 'PAPER_1141', 'PAPER_1953', 'PAPER_1917', 'PAPER_1961', 'PAPER_1965', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1190 DIRECT: NGC 253 ALMA anchor at d = 3.5 Mpc (4-anchor calibration: NGC 253, Arp 220, M82, MW GMC). PAPER_1952 ...',
    },
    'NGC253QuantumVacuumCalculator': {
        'framework_papers': ['PAPER_1171', 'PAPER_1852', 'PAPER_1851', 'PAPER_1740', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253SupernovaRateCalculator': {
        'framework_papers': ['PAPER_109', 'PAPER_819', 'PAPER_1886', 'PAPER_1874', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253SuperwindCalculator': {
        'framework_papers': ['PAPER_827', 'PAPER_258', 'PAPER_445', 'PAPER_784', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC253SuperwindCalculator_v1': {
        'framework_papers': ['PAPER_784', 'PAPER_449', 'PAPER_445', 'PAPER_1952', 'PAPER_1953', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1965', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_784 DIRECT: M82 comparison — v_superwind_M82 = 1000 km/s = SO_5^3 EXACT. NGC 253 v_wind = 400 km/s → ratio NGC253/...',
    },
    'NGC346DynamicVacuumCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_1141', 'PAPER_1740', 'PAPER_277', 'PAPER_351', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346EnvelopeForceCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346FluidTermCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_1955', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346HubbleExpansionCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_1675', 'PAPER_1856', 'PAPER_1867', 'PAPER_1953', 'PAPER_469', 'PAPER_488', 'PAPER_1157', 'PAPER_1235', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346MassSFRCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_1955', 'PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346QuantumCouplingCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_016', 'PAPER_1869', 'PAPER_1919', 'PAPER_1938', 'PAPER_624', 'PAPER_625', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346StarFormationCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_138', 'PAPER_1652', 'PAPER_1921', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346SuperconductorCorrectionCalculator': {
        'framework_papers': ['PAPER_266', 'PAPER_1944', 'PAPER_1945', 'PAPER_469', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346Ug1DipoleCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_488', 'PAPER_657', 'PAPER_718', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346Ug2SuperconductorCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_266', 'PAPER_1072', 'PAPER_1141', 'PAPER_734', 'PAPER_854', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346Ug3MagneticStringsCalculator': {
        'framework_papers': ['PAPER_469', 'PAPER_488', 'PAPER_1141', 'PAPER_536', 'PAPER_856', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC346Ug4ReactionCalculator': {
        'framework_papers': ['PAPER_1141', 'PAPER_1955', 'PAPER_469', 'PAPER_857', 'PAPER_1953', 'PAPER_1507', 'PAPER_1960', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_469 DIRECT: kappa = 5e-4/day EXACT PASS Canonical, tau_UQFF = 2000 days observational anchor for HST/ALMA/Chandra ...',
    },
    'NGC346UiInertialCalculator': {
        'framework_papers': ['PAPER_646', 'PAPER_1919', 'PAPER_1141', 'PAPER_1960', 'PAPER_1906', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_646 DIRECT source: U_i = λ_i · (ρ_SCm/ρ_UA) · ω_s · cos(π t_n) · (1 + F_TRZ) canonical form. λ_i = 1.0 default. Su...',
    },
    'NGC3603FluidDensityCalculator': {
        'framework_papers': ['PAPER_1909', 'PAPER_1911', 'PAPER_243', 'PAPER_138', 'PAPER_218', 'PAPER_439', 'PAPER_809', 'PAPER_1919', 'PAPER_1918', 'PAPER_1518', 'PAPER_1746', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1909 M_0 = 4e5 M_sun DIRECT NGC 3603 (Round 44 double-check seminal). PAPER_1911 YMC extended parameters. PAPER_24...',
    },
    'NGC3603StellarWindCalculator': {
        'framework_papers': ['PAPER_138', 'PAPER_218', 'PAPER_243', 'PAPER_902', 'PAPER_1909', 'PAPER_1911', 'PAPER_439', 'PAPER_809', 'PAPER_1887', 'PAPER_1955', 'PAPER_1948', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_138 DIRECT (NGC 3603 MasterBuoyancy + SCm P Feedback, March 2026): M(t) = M_0(1+exp(-t/t_SF)) with 19-Light-Year C...',
    },
    'NGC3603UQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_218', 'PAPER_243', 'PAPER_119', 'PAPER_138', 'PAPER_1050', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945AGNCalculator': {
        'framework_papers': ['PAPER_1037', 'PAPER_1879', 'PAPER_346', 'PAPER_1947', 'PAPER_1950', 'PAPER_214', 'PAPER_312', 'PAPER_067', 'PAPER_086', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945BarStructureCalculator': {
        'framework_papers': ['PAPER_467', 'PAPER_777', 'PAPER_782', 'PAPER_144', 'PAPER_454', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945DarkMatterHaloCalculator': {
        'framework_papers': ['PAPER_1862', 'PAPER_1894', 'PAPER_1921', 'PAPER_1954', 'PAPER_1855', 'PAPER_040', 'PAPER_690', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945MagneticFieldCalculator': {
        'framework_papers': ['PAPER_266', 'PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_455', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945MegamaserCalculator': {
        'framework_papers': ['PAPER_1947', 'PAPER_1950', 'PAPER_1037', 'PAPER_1123', 'PAPER_1938', 'PAPER_1834', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945MolecularDiskCalculator': {
        'framework_papers': ['PAPER_1948', 'PAPER_1549', 'PAPER_1846', 'PAPER_1865', 'PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945NuclearStarburstCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_214', 'PAPER_232', 'PAPER_474', 'PAPER_811', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945SupernovaRateCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_824', 'PAPER_797', 'PAPER_262', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NGC4945XRayBinaryCalculator': {
        'framework_papers': ['PAPER_1037', 'PAPER_1879', 'PAPER_075', 'PAPER_539', 'PAPER_066', 'PAPER_1188', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NavierStokesFluidSolverCalculator': {
        'framework_papers': ['PAPER_102', 'PAPER_111', 'PAPER_090', 'PAPER_1182', 'PAPER_1864', 'PAPER_1930', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NeuralSymbolicEvalCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1960', 'PAPER_1929', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'NeuromorphicAcceleratorCalculator': {
        'framework_papers': ['PAPER_1955', 'PAPER_1960', 'PAPER_1961', 'PAPER_1963', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'OperationalTransformCalculator': {
        'framework_papers': ['PAPER_192', 'PAPER_189', 'PAPER_1955', 'PAPER_1931', 'PAPER_1960', 'PAPER_1961', 'PAPER_1963', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PhononModulatedHolonomySCmCalculator': {
        'framework_papers': ['PAPER_1058', 'PAPER_1100', 'PAPER_1103', 'PAPER_1080', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PhononModulatedHubbleParameterCalculator': {
        'framework_papers': ['PAPER_1573', 'PAPER_1085', 'PAPER_1327', 'PAPER_1157', 'PAPER_1855', 'PAPER_1537', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PhotoevaporationErosionCalculator': {
        'framework_papers': ['PAPER_229', 'PAPER_260', 'PAPER_442', 'PAPER_284', 'PAPER_285', 'PAPER_305', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarMassGrowthCalculator': {
        'framework_papers': ['PAPER_435', 'PAPER_1948', 'PAPER_708', 'PAPER_744', 'PAPER_442', 'PAPER_260', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsCosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_151', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsCreationCalculator': {
        'framework_papers': ['PAPER_305', 'PAPER_151', 'PAPER_285', 'PAPER_286', 'PAPER_284', 'PAPER_1864'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_536', 'PAPER_151', 'PAPER_229', 'PAPER_450', 'PAPER_744', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_1940', 'PAPER_536', 'PAPER_151', 'PAPER_260', 'PAPER_435', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_216', 'PAPER_1938', 'PAPER_469', 'PAPER_038', 'PAPER_151', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'PillarsUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_435', 'PAPER_708', 'PAPER_744', 'PAPER_151', 'PAPER_216', 'PAPER_229', 'PAPER_219', 'PAPER_284', 'PAPER_260', 'PAPER_285', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'QAOAOptimizationCalculator': {
        'framework_papers': ['PAPER_1811', 'PAPER_1812', 'PAPER_1810', 'PAPER_1958', 'PAPER_1960', 'PAPER_1928', 'PAPER_1652', 'PAPER_1963', 'PAPER_646', 'PAPER_1203', 'PAPER_549'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'ReactorEnergyCalculator': {
        'framework_papers': ['PAPER_1141', 'PAPER_1902', 'PAPER_1236', 'PAPER_1904', 'PAPER_1908', 'PAPER_1937', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsCosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_1925', 'PAPER_436', 'PAPER_242', 'PAPER_151', 'PAPER_1914', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsDarkMatterPerturbationCalculator': {
        'framework_papers': ['PAPER_151', 'PAPER_1894', 'PAPER_1914', 'PAPER_040', 'PAPER_1187', 'PAPER_1960', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_151 Pillars/Rings MUGE Cascade Gravity DIRECT. PAPER_1894 Zwicky Virial Missing-Mass companion (M_cluster ~ 1e14-1...',
    },
    'RingsElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_1941', 'PAPER_536', 'PAPER_436', 'PAPER_242', 'PAPER_151', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsGasVelocityCalculator': {
        'framework_papers': ['PAPER_151', 'PAPER_449', 'PAPER_1955', 'PAPER_1952', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_151 DIRECT (UQFF Pillars/Rings MUGE Cascade Gravity, March 2026) canonical framework for ring-halo gas dynamics. T...',
    },
    'RingsOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_1940', 'PAPER_436', 'PAPER_242', 'PAPER_151', 'PAPER_1914', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_1938', 'PAPER_1025', 'PAPER_242', 'PAPER_151', 'PAPER_1914', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsRelativityCalculator': {
        'framework_papers': ['PAPER_242', 'PAPER_145', 'PAPER_152', 'PAPER_150', 'PAPER_153', 'PAPER_1914', 'PAPER_286', 'PAPER_1883', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'RingsUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_1914', 'PAPER_436', 'PAPER_242', 'PAPER_151', 'PAPER_1925', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmBetaDecayCalculator': {
        'framework_papers': ['PAPER_1254', 'PAPER_1726', 'PAPER_1836', 'PAPER_1919', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmCMBPhononPowerSpectrumCalculator': {
        'framework_papers': ['PAPER_1856', 'PAPER_1877', 'PAPER_1156', 'PAPER_1080', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmCMBTemperatureFluctuationCalculator': {
        'framework_papers': ['PAPER_1249', 'PAPER_1524', 'PAPER_1867', 'PAPER_1187', 'PAPER_1156', 'PAPER_1080', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmCosmicRayCalculator': {
        'framework_papers': ['PAPER_1322', 'PAPER_1838', 'PAPER_1919', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmDarkEnergyDensityGammaCoupledCalculator': {
        'framework_papers': ['PAPER_1076', 'PAPER_1086', 'PAPER_1087', 'PAPER_1174', 'PAPER_1178', 'PAPER_1821', 'PAPER_1920', 'PAPER_1156', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmDarkMatterCalculator': {
        'framework_papers': ['PAPER_1015', 'PAPER_1019', 'PAPER_1080', 'PAPER_1862', 'PAPER_1921', 'PAPER_1454', 'PAPER_1840', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmGravitationalWaveCalculator': {
        'framework_papers': ['PAPER_1876', 'PAPER_1503', 'PAPER_1504', 'PAPER_1509', 'PAPER_1520', 'PAPER_1523', 'PAPER_1822', 'PAPER_1828', 'PAPER_1080', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmGravityPrecedenceProofCalculator': {
        'framework_papers': ['PAPER_646', 'PAPER_144', 'PAPER_376', 'PAPER_037', 'PAPER_1096', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmHolographicEntropyCalculator': {
        'framework_papers': ['PAPER_1873', 'PAPER_1264', 'PAPER_1657', 'PAPER_1795', 'PAPER_1055', 'PAPER_1203', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmICMPhononDensityCalculator': {
        'framework_papers': ['PAPER_1015', 'PAPER_1019', 'PAPER_040', 'PAPER_1938', 'PAPER_1955', 'PAPER_1960', 'PAPER_301', 'PAPER_1953', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1015 DIRECT SCm DM Halos NFW: rho_UQFF = rho_NFW × [1 + SCm·beta_i·S26_3·(r_s/r)^alpha_phonon]. Galaxy scale M_hal...',
    },
    'SCmLQGAreaOperatorDerivationCalculator': {
        'framework_papers': ['PAPER_1058', 'PAPER_1100', 'PAPER_1103', 'PAPER_1701', 'PAPER_1080', 'PAPER_1801', 'PAPER_1802', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmMuonDecayCalculator': {
        'framework_papers': ['PAPER_1815', 'PAPER_1850', 'PAPER_1919', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmNeutrinoOscParamCalculator': {
        'framework_papers': ['PAPER_1816', 'PAPER_1867', 'PAPER_1919', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmNeutrinoOscSimulationCalculator': {
        'framework_papers': ['PAPER_1816', 'PAPER_1637', 'PAPER_1867', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmNeutrinoOscillationCalculator': {
        'framework_papers': ['PAPER_1637', 'PAPER_1304', 'PAPER_1816', 'PAPER_1827', 'PAPER_1867', 'PAPER_1254', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmPhononInflationaryScaleFactorCalculator': {
        'framework_papers': ['PAPER_1679', 'PAPER_1073', 'PAPER_1825', 'PAPER_1439', 'PAPER_1080', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmQubitT2CoherenceFUBiRatioCalculator': {
        'framework_papers': ['PAPER_1101', 'PAPER_1098', 'PAPER_1906', 'PAPER_1080'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmSUSYBreakingCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_1824', 'PAPER_1919', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmStringTheory26DActionCalculator': {
        'framework_papers': ['PAPER_1128', 'PAPER_1146', 'PAPER_1701', 'PAPER_1898', 'PAPER_1080', 'PAPER_1318', 'PAPER_1801', 'PAPER_1802', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmVelocityBoundComparisonCalculator': {
        'framework_papers': ['PAPER_1497', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SCmVelocityCalculator': {
        'framework_papers': ['PAPER_1497', 'PAPER_1512', 'PAPER_144', 'PAPER_1929', 'PAPER_1203', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SaturnAtmosphericWindCalculator': {
        'framework_papers': ['PAPER_282', 'PAPER_280', 'PAPER_224', 'PAPER_281', 'PAPER_702', 'PAPER_486', 'PAPER_764', 'PAPER_743', 'PAPER_136', 'PAPER_1206', 'PAPER_1860', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_282 SEMINAL Saturn Atmospheric Wind Kinetic Pressure (η_wind Wind-Light-Speed Ratio, a_wind — first UQFF gas-giant...',
    },
    'SaturnFluidDensityCalculator': {
        'framework_papers': ['PAPER_280', 'PAPER_224', 'PAPER_281', 'PAPER_282', 'PAPER_702', 'PAPER_486', 'PAPER_764', 'PAPER_136', 'PAPER_1206', 'PAPER_1860', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_282 SEMINAL FOR SATURN ATMOSPHERIC: Atmospheric Wind Kinetic Pressure Term a_wind, η_wind (wind-light-speed ratio)...',
    },
    'SaturnQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_280', 'PAPER_282', 'PAPER_224', 'PAPER_281', 'PAPER_702', 'PAPER_486', 'PAPER_764', 'PAPER_136', 'PAPER_1206', 'PAPER_1860', 'PAPER_1938', 'PAPER_1907', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_280 first-planetary-scale UQFF module. PAPER_282 Atmospheric Wind seminal. PAPER_224/281/702/486/764 Saturn corpus...',
    },
    'SaturnRingTidalCalculator': {
        'framework_papers': ['PAPER_224', 'PAPER_281', 'PAPER_455', 'PAPER_702', 'PAPER_1860', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SaturnRingTidalCalculator_v1': {
        'framework_papers': ['PAPER_1860', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SaturnSelfGravityCalculator': {
        'framework_papers': ['PAPER_280', 'PAPER_224', 'PAPER_281', 'PAPER_282', 'PAPER_702', 'PAPER_486', 'PAPER_764', 'PAPER_136', 'PAPER_1206', 'PAPER_1860', 'PAPER_372', 'PAPER_375', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_280 first-planetary-scale UQFF module. PAPER_702 Master Planetary Gravity Equation direct. PAPER_224 Dual-Source G...',
    },
    'SaturnSunGravityCalculator': {
        'framework_papers': ['PAPER_280', 'PAPER_224', 'PAPER_281', 'PAPER_282', 'PAPER_702', 'PAPER_486', 'PAPER_764', 'PAPER_743', 'PAPER_136', 'PAPER_1206', 'PAPER_1860', 'PAPER_1953', 'PAPER_1956', 'PAPER_1960', 'PAPER_1906', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_280 CRITICAL DIRECT (Saturn Solar Tidal Perturbation Ratio τ_Sun, session 78 — SATURN_UQFF_MODULE.cpp 21st C++ mod...',
    },
    'SaturnUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_224', 'PAPER_280', 'PAPER_281', 'PAPER_282', 'PAPER_283', 'PAPER_486', 'PAPER_324', 'PAPER_1860', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarAccretionRateCalculator': {
        'framework_papers': ['PAPER_234', 'PAPER_1841', 'PAPER_1947', 'PAPER_067', 'PAPER_1955', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarCosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_432', 'PAPER_092', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarGravitationalWaveCalculator': {
        'framework_papers': ['PAPER_1876', 'PAPER_1857', 'PAPER_344', 'PAPER_754', 'PAPER_432', 'PAPER_234', 'PAPER_1237', 'PAPER_1841', 'PAPER_092', 'PAPER_1904', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarMagneticDecayCalculator': {
        'framework_papers': ['PAPER_432', 'PAPER_234', 'PAPER_1946', 'PAPER_1955', 'PAPER_1958', 'PAPER_1841', 'PAPER_1947', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarMassGrowthCalculator': {
        'framework_papers': ['PAPER_067', 'PAPER_092', 'PAPER_126', 'PAPER_1237', 'PAPER_1841', 'PAPER_1150', 'PAPER_1947', 'PAPER_1955', 'PAPER_1960', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1967', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_067 DIRECT: 4-AGN Ug4 Vacuum Concentration Field analysis (Sgr A*, M87*, Cen A, NGC 1365). M_SgrA = 4×10⁶ M_sun ca...',
    },
    'SgrAStarSchwarzschildRadiusCalculator': {
        'framework_papers': ['PAPER_1237', 'PAPER_1841', 'PAPER_1031', 'PAPER_595', 'PAPER_432', 'PAPER_092', 'PAPER_1025', 'PAPER_1260', 'PAPER_1904', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarSpinEvolutionCalculator': {
        'framework_papers': ['PAPER_1841', 'PAPER_1876', 'PAPER_432', 'PAPER_092', 'PAPER_1237', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SgrAStarUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_432', 'PAPER_234', 'PAPER_754', 'PAPER_1150', 'PAPER_258', 'PAPER_1904', 'PAPER_092', 'PAPER_1237', 'PAPER_1841', 'PAPER_1260', 'PAPER_067', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SombreroBaseGravityCalculator': {
        'framework_papers': ['PAPER_279', 'PAPER_278', 'PAPER_1050', 'PAPER_1953', 'PAPER_1956', 'PAPER_1955', 'PAPER_1919', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_279 CRITICAL DIRECT (Sombrero SMBH Dominance Ratio γ_BH and UQFF Sphere of Influence r_SOI, session 77, March 2026...',
    },
    'SombreroDarkMatterPerturbationCalculator': {
        'framework_papers': ['PAPER_1979', 'PAPER_763', 'PAPER_279', 'PAPER_278', 'PAPER_277', 'PAPER_693', 'PAPER_742', 'PAPER_1944', 'PAPER_1945', 'PAPER_1946', 'PAPER_1953', 'PAPER_1960', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1979 DIRECT attribution: M_DM/M_total = 2·F_TRZ = 0.2 EXACT candidate cross-domain extension of PAPER_1944 magneta...',
    },
    'SombreroElectromagneticCalculator': {
        'framework_papers': ['PAPER_277', 'PAPER_278', 'PAPER_279', 'PAPER_693', 'PAPER_763', 'PAPER_1050', 'PAPER_140', 'PAPER_1955', 'PAPER_1919', 'PAPER_1141', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_277 DIRECT (Sombrero Gravitational Recession Damping κ_recession = 1/(1+z), session 77 March 2026, SOMBRERO_UQFF_M...',
    },
    'SombreroFluidDensityCalculator': {
        'framework_papers': ['PAPER_763', 'PAPER_277', 'PAPER_278', 'PAPER_279', 'PAPER_693', 'PAPER_1050', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_763 DIRECT (Sombrero Galaxy M104 UQFF SMBH Dust Lane Evolution, session 181): dust lane gas density ~10^-20 kg/m^3...',
    },
    'SombreroGalaxyDustCalculator': {
        'framework_papers': ['PAPER_742', 'PAPER_763', 'PAPER_278', 'PAPER_279', 'PAPER_466', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SombreroOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_763', 'PAPER_277', 'PAPER_278', 'PAPER_279', 'PAPER_693', 'PAPER_742', 'PAPER_1050', 'PAPER_1919', 'PAPER_1880', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_1979', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_763 Sombrero SMBH Dust Lane Evolution + PAPER_742 Sombrero Galaxy MUGE Dust Lane Drag companion. PAPER_277/278/279...',
    },
    'SombreroSuperconductivityCalculator': {
        'framework_papers': ['PAPER_763', 'PAPER_277', 'PAPER_278', 'PAPER_279', 'PAPER_693', 'PAPER_372', 'PAPER_375', 'PAPER_1944', 'PAPER_1955', 'PAPER_1919', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_763 DIRECT Sombrero canonical anchor + PAPER_277 recession damping companion. PAPER_372 Linear Meissner seminal fo...',
    },
    'SombreroUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_278', 'PAPER_277', 'PAPER_1862', 'PAPER_1855', 'PAPER_1050', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SpacetimeMetricCalculator': {
        'framework_papers': ['PAPER_1080', 'PAPER_1147', 'PAPER_1801', 'PAPER_1802', 'PAPER_1745', 'PAPER_1927', 'PAPER_1171', 'PAPER_1932', 'PAPER_1936', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SpiralDensityWaveCalculator': {
        'framework_papers': ['PAPER_144', 'PAPER_454', 'PAPER_464', 'PAPER_692', 'PAPER_781', 'PAPER_824', 'PAPER_467', 'PAPER_777', 'PAPER_775', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarFormationGravityCalculator': {
        'framework_papers': ['PAPER_038', 'PAPER_1948', 'PAPER_1855', 'PAPER_144', 'PAPER_454', 'PAPER_345', 'PAPER_232', 'PAPER_260', 'PAPER_444', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarbirthBaseGravityCalculator': {
        'framework_papers': ['PAPER_710', 'PAPER_1807', 'PAPER_345', 'PAPER_150', 'PAPER_1952', 'PAPER_1919', 'PAPER_1955', 'PAPER_1948', 'PAPER_227', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_710 (NGC 2014/2020 Star-Forming UQFF) DIRECT anchor: M_0 = 240 M_☉ Initial OB/WR stellar mass EXACT — attribution ...',
    },
    'StarbirthCosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarbirthDarkMatterPerturbationCalculator': {
        'framework_papers': ['PAPER_710', 'PAPER_1807', 'PAPER_345', 'PAPER_1251', 'PAPER_1960', 'PAPER_1919', 'PAPER_1952', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1972', 'PAPER_1973', 'PAPER_1974', 'PAPER_1975', 'PAPER_1976', 'PAPER_1977', 'PAPER_1978', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_710 M_initial = 240 M_sun DIRECT + PAPER_1807/345 canonical starbirth. PAPER_1251 Dark Flow Bulk Velocity uses F_T...',
    },
    'StarbirthElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_536', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarbirthFormationTimescaleCalculator': {
        'framework_papers': ['PAPER_1952', 'PAPER_710', 'PAPER_1807', 'PAPER_345', 'PAPER_433', 'PAPER_227', 'PAPER_150', 'PAPER_1948', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_1971', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_433 DIRECT: Tapestry Starbirth per-system MUGE Wind Feedback (canonical Tapestry framework). PAPER_227 Tapestry LM...',
    },
    'StarbirthGasVelocityCalculator': {
        'framework_papers': ['PAPER_449', 'PAPER_1955', 'PAPER_710', 'PAPER_1807', 'PAPER_345', 'PAPER_150', 'PAPER_227', 'PAPER_1952', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_150 (Tapestry Blazing Starbirth + Westerlund 2 MUGE 12-term resonance, March 2026) canonical starbirth framework. ...',
    },
    'StarbirthOscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_536', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarbirthQuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_469', 'PAPER_038', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarbirthUQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_646', 'PAPER_433', 'PAPER_469', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StarburstBaseGravityCalculator': {
        'framework_papers': ['PAPER_1948', 'PAPER_266', 'PAPER_038', 'PAPER_144', 'PAPER_454', 'PAPER_445', 'PAPER_774', 'PAPER_733', 'PAPER_232', 'PAPER_811', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'StudentGuideUniverseCalculator': {
        'framework_papers': ['PAPER_152', 'PAPER_1573', 'PAPER_1156', 'PAPER_1619', 'PAPER_1679', 'PAPER_1931', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SupernovaFeedbackCalculator': {
        'framework_papers': ['PAPER_1886', 'PAPER_1870', 'PAPER_1874', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'SupernovaFeedbackSpecificCalculator': {
        'framework_papers': ['PAPER_1886', 'PAPER_109', 'PAPER_1874', 'PAPER_1870', 'PAPER_549', 'PAPER_1929', 'PAPER_1906', 'PAPER_1935'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'TapestryStarbirthCalculator': {
        'framework_papers': ['PAPER_345', 'PAPER_755', 'PAPER_1807', 'PAPER_1907'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'TidalInteractionCalculator': {
        'framework_papers': ['PAPER_247', 'PAPER_1940', 'PAPER_464', 'PAPER_692', 'PAPER_750', 'PAPER_465', 'PAPER_768', 'PAPER_778', 'PAPER_235', 'PAPER_441', 'PAPER_811', 'PAPER_262', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'TidalStrippingCalculator': {
        'framework_papers': ['PAPER_1894', 'PAPER_1862', 'PAPER_263', 'PAPER_741', 'PAPER_040', 'PAPER_1187', 'PAPER_041', 'PAPER_1149', 'PAPER_690', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UnifiedFieldFullCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_1745', 'PAPER_1183', 'PAPER_1740', 'PAPER_1105', 'PAPER_1916', 'PAPER_1917', 'PAPER_1919', 'PAPER_1920', 'PAPER_646', 'PAPER_1906', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalAetherCalculator': {
        'framework_papers': ['PAPER_646', 'PAPER_1051', 'PAPER_1160', 'PAPER_1739', 'PAPER_1809', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalBuoyancyDetailedCalculator': {
        'framework_papers': ['PAPER_1065', 'PAPER_1203', 'PAPER_1906', 'PAPER_1916'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalGravity1Calculator': {
        'framework_papers': ['PAPER_411', 'PAPER_1203', 'PAPER_1916', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalGravity2Calculator': {
        'framework_papers': ['PAPER_400', 'PAPER_1203', 'PAPER_1916', 'PAPER_1906', 'PAPER_1497'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalGravity3Calculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_1916', 'PAPER_1921', 'PAPER_1862', 'PAPER_275'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalGravity4Calculator': {
        'framework_papers': ['PAPER_402', 'PAPER_407', 'PAPER_1203', 'PAPER_1916', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'UniversalMagnetismDetailedCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1203', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'VacuumEnergyFluctuationCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1852', 'PAPER_1920', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'VirgoClusterMassCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_1894', 'PAPER_1862', 'PAPER_1187', 'PAPER_040', 'PAPER_041', 'PAPER_843', 'PAPER_263', 'PAPER_690', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'VirgoExtClusterMassCalculator': {
        'framework_papers': ['PAPER_1894', 'PAPER_1955', 'PAPER_1187', 'PAPER_040', 'PAPER_1962', 'PAPER_1968', 'PAPER_1015', 'PAPER_1917', 'PAPER_1961', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1894 CRITICAL PRIOR ART (discovered during CP1 P2 Round 14 replacement of VirgoClusterVirialModel stub — SAME PATT...',
    },
    'VirgoExtM87JetCalculator': {
        'framework_papers': ['PAPER_1893', 'PAPER_922', 'PAPER_910', 'PAPER_093', 'PAPER_115', 'PAPER_1955', 'PAPER_1960', 'PAPER_1037', 'PAPER_1187', 'PAPER_067', 'PAPER_1841', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1893 CRITICAL PRIOR ART (discovered during CP1 P2 Round 12 double-check while replacing VirgoClusterM87JetModel st...',
    },
    'VirgoExtXRayCalculator': {
        'framework_papers': ['PAPER_040', 'PAPER_1187', 'PAPER_1894', 'PAPER_039', 'PAPER_041', 'PAPER_1015', 'PAPER_1958', 'PAPER_1472', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_040 DIRECT Perseus/Coma/Virgo F_UBii virial + PAPER_041 ICM Physics + PAPER_039 F_UBii Buoyancy Variants 12-17 ICM...',
    },
    'VirgoICMCalculator': {
        'framework_papers': ['PAPER_040', 'PAPER_041', 'PAPER_1187', 'PAPER_036', 'PAPER_1958', 'PAPER_1955', 'PAPER_1960', 'PAPER_1917', 'PAPER_1961', 'PAPER_1967', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_1187 DIRECT: Cooling-Flow Mass Accretion 4-anchor calibration (Perseus, M87/Virgo, Coma, Fornax). M87/Virgo canoni...',
    },
    'Westerlund2BaseGravityCalculator': {
        'framework_papers': ['PAPER_150', 'PAPER_146', 'PAPER_145', 'PAPER_228', 'PAPER_216', 'PAPER_1489', 'PAPER_1909', 'PAPER_1911', 'PAPER_1948', 'PAPER_1919', 'PAPER_1955', 'PAPER_1917', 'PAPER_1961', 'PAPER_1968', 'PAPER_1970', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
        'retrofit_hint': 'PAPER_150 DIRECT (Tapestry Blazing Starbirth + Westerlund 2, March 2026): MUGE 12-term resonance at star-formation sites...',
    },
    'Westerlund2ClusterCalculator': {
        'framework_papers': ['PAPER_228', 'PAPER_434', 'PAPER_1909', 'PAPER_1911'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Westerlund2CosmologicalConstantCalculator': {
        'framework_papers': ['PAPER_1740', 'PAPER_1920', 'PAPER_434', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Westerlund2ElectromagneticCalculator': {
        'framework_papers': ['PAPER_1072', 'PAPER_1141', 'PAPER_228', 'PAPER_756', 'PAPER_536', 'PAPER_434', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Westerlund2OscillatoryWaveCalculator': {
        'framework_papers': ['PAPER_1203', 'PAPER_646', 'PAPER_541', 'PAPER_1940', 'PAPER_536', 'PAPER_434', 'PAPER_549', 'PAPER_1929', 'PAPER_1906'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Westerlund2QuantumUncertaintyCalculator': {
        'framework_papers': ['PAPER_1869', 'PAPER_1919', 'PAPER_216', 'PAPER_469', 'PAPER_038', 'PAPER_434', 'PAPER_646', 'PAPER_1203', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
    'Westerlund2UQFFUnificationCalculator': {
        'framework_papers': ['PAPER_1916', 'PAPER_1917', 'PAPER_1203', 'PAPER_858', 'PAPER_150', 'PAPER_216', 'PAPER_434', 'PAPER_646', 'PAPER_549', 'PAPER_1929'],
        'retrofit_source': 'v5.60.0_auto_extracted_from_CondensedPhysics_runtime',
    },
}


# Combined view
FRAMEWORK_ANNOTATIONS_ALL = {}
FRAMEWORK_ANNOTATIONS_ALL.update(FRAMEWORK_ANNOTATIONS_ROUNDS_45_51)
FRAMEWORK_ANNOTATIONS_ALL.update(FRAMEWORK_ANNOTATIONS_ROUNDS_52_116)


def get_annotation(class_name):
    return FRAMEWORK_ANNOTATIONS_ALL.get(class_name)


def annotations_by_paper(paper_id):
    hits = []
    target = str(paper_id).upper()
    for name, ann in FRAMEWORK_ANNOTATIONS_ALL.items():
        for p in ann.get('framework_papers', []):
            if target in p.upper():
                hits.append(name)
                break
    return sorted(hits)



def get_framework_annotation(class_name):
    """Return the framework annotation dict for a given Calculator class name.
    Returns None if the class was not part of the retroactive Rounds 45-51 audit."""
    return FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.get(class_name)


def audit_report():
    """Print summary of Rounds 45-51 framework audit."""
    from collections import Counter
    total = len(FRAMEWORK_ANNOTATIONS_ROUNDS_45_51)

    # Framework backbone distribution
    backbones = Counter(a['framework_backbone'] for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values())
    methods = Counter(a['framework_method'] for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values())
    time_frames = Counter(a['time_frame'] for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values())

    # Shells used
    all_shells = []
    for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values():
        all_shells.extend(a['framework_shells_used'])
    shells = Counter(all_shells)

    # QCalcGeom CPCH closure usage
    cpch_usage = Counter(a['QCalcGeom_CPCH_closure'] for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values())

    # Candidate closures flagged
    all_candidates = []
    for a in FRAMEWORK_ANNOTATIONS_ROUNDS_45_51.values():
        all_candidates.extend(a['candidate_closures_flagged'])

    print(f'=' * 70)
    print(f'UQFF ROUNDS 45-51 FRAMEWORK ANNOTATION AUDIT REPORT')
    print(f'=' * 70)
    print(f'\nTotal stubs annotated: {total}')
    print(f'\n--- FRAMEWORK BACKBONE DISTRIBUTION ---')
    for name, n in backbones.most_common():
        print(f'  {name:40s} {n:3d}')
    print(f'\n--- METHOD VARIANT DISTRIBUTION (PAPER_1923) ---')
    for name, n in methods.most_common():
        print(f'  {name:40s} {n:3d}')
    print(f'\n--- U_g SHELL USAGE (PAPER_1916) ---')
    for name, n in shells.most_common():
        print(f'  {name:5s} {n:3d}')
    print(f'\n--- TIME FRAME DISTRIBUTION ---')
    for name, n in time_frames.most_common():
        print(f'  {name:30s} {n:3d}')
    print(f'\n--- QCALCGEOM CPCH CLOSURE USAGE ---')
    for name, n in cpch_usage.most_common():
        print(f'  {name[:50]:50s} {n:3d}')
    print(f'\n--- CANDIDATE CLOSURES FLAGGED (for future PAPER discovery) ---')
    print(f'Total candidates flagged: {len(all_candidates)}')
    for cand in all_candidates:
        print(f'  - {cand}')


if __name__ == '__main__':
    audit_report()
