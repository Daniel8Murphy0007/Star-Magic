from pathlib import Path
import re

samples = [
    'whitepapers/PAPER_040_XRay_Cluster_Buoyancy_Perseus_Coma_Virgo.md',
    'whitepapers/PAPER_047_Nuclear_Binding_Energy_26Level_Polynomial.md',
    'whitepapers/PAPER_070_GRB_Afterglow_UQFF_Vacuum_Drag.md',
    'whitepapers/PAPER_095_Hawking_Radiation_UQFF_Corrections.md',
    'whitepapers/PAPER_100_BEC_UQFF_Condensate_Phase_Transition.md',
    'whitepapers/PAPER_200_Yang_Mills_Mass_Gap_UQFF.md',
    'whitepapers/PAPER_300_Riemann_Hypothesis_UQFF_Zeta.md',
    'whitepapers/PAPER_500_LENR_Lattice_Confinement_UQFF.md',
]
for fpath in samples:
    try:
        c = Path(fpath).read_text(encoding='utf-8', errors='replace')
        m_tags = re.search(r'^tags:\s*\[([^\]]+)\]', c, re.MULTILINE)
        tags = m_tags.group(1) if m_tags else 'none'
        print(Path(fpath).name + ' | ' + tags[:100])
    except Exception as e:
        print('MISSING: ' + fpath + ' | ' + str(e))
