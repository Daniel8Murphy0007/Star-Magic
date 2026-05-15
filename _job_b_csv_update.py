"""Update _job_b_categorization.csv with completion status for the 10 updated papers,
and reclassify them from bucket H to bucket H* (Damping/GW-application papers requiring
G1-G8 ledger closure block per Session 262 audit)."""
import csv, os

CSV = '_job_b_categorization.csv'
DONE = {
    'PAPER_001_GW170817_UQFF_Damping_Analysis.md',
    'PAPER_002_GW190425_Mass_Gap_Interpretation.md',
    'PAPER_003_GW150914_UQFF_vs_LIGO_Strain.md',
    'PAPER_004_GW170817_BNS_Chirp_Phase_Evolution.md',
    'PAPER_005_BH_Merger_Energy_Retention_UQFF.md',
    'PAPER_006_GW170817_Multi_Messenger_Full_Inspiral.md',
    'PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF.md',
    'PAPER_008_UQFF_Waveform_Phase_Evolution_Template_Mismatch.md',
    'PAPER_008b_Full_Inspiral_Waveform_UQFF.md',
    'PAPER_009_Damping_Mechanism_Decomposition_UQFF.md',
}

rows = []
with open(CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames + ['status', 'session', 'commit']
    for r in reader:
        if r['filename'] in DONE:
            r['bucket'] = 'H*'  # reclassified: damping/GW-application papers DO need G1-G8 + ledger block
            r['suggested_action'] = 'DONE - G1-G8 + 27-decade ledger block inserted; tailored P-suite hook'
            r['status'] = 'DONE v5.78'
            r['session'] = '262'
            r['commit'] = '7fadd27d'
        else:
            r.setdefault('status', '')
            r.setdefault('session', '')
            r.setdefault('commit', '')
        rows.append(r)

with open(CSV, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f'Updated {len(DONE)} rows -> bucket H*, status=DONE v5.78')
print(f'Total rows: {len(rows)}')
