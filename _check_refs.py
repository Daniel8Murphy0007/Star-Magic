import re
files = [
    'whitepapers/PAPER_003_GW150914_UQFF_vs_LIGO_Strain.md',
    'whitepapers/PAPER_004_GW170817_BNS_Chirp_Phase_Evolution.md',
    'whitepapers/PAPER_005_BH_Merger_Energy_Retention_UQFF.md',
    'whitepapers/PAPER_069_Radio_Transient_ASKAP_J1832_UQFF.md',
    'whitepapers/PAPER_092_SgrA_MUGE_Comparison.md',
    'whitepapers/PAPER_094_Magnetar_SGR1745_UQFF_Calibration.md',
    'whitepapers/PAPER_228_Westerlund2_OB_StellarWind_MUGE.md',
    'whitepapers/PAPER_660_WhiteHoleRadiationUQFF.md',
    'whitepapers/PAPER_694_CrabNebulaPWNUQFF.md',
]
for f in files:
    c = open(f, encoding='utf-8', errors='replace').read()
    m = re.search(r'## References(.*?)(?=\n## |\Z)', c, re.DOTALL)
    if m:
        ref = m.group(0)[:600]
        name = f.split('/')[-1]
        print(f'=== {name} ===')
        print(ref)
        print()
