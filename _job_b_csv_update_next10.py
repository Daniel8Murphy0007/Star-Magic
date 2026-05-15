"""Mark next 10 papers DONE in _job_b_categorization.csv (commit forthcoming)."""
import csv

CSV = '_job_b_categorization.csv'
DONE = {
    'PAPER_010_Post_Merger_Oscillations_Remnant_Mass_UQFF.md',
    'PAPER_010b_Time_Domain_Chirp_23Hz_UQFF.md',
    'PAPER_011_Stochastic_GW_Background_UQFF_Implications.md',
    'PAPER_011b_Amplitude_Reduction_Factor_UQFF.md',
    'PAPER_012_Eccentric_Binary_Circularization_UQFF.md',
    'PAPER_012b_GW150914_Waveform_Validation.md',
    'PAPER_013_Magnetar_Spin_Down_UQFF_Framework.md',
    'PAPER_013b_LISA_SMBH_Merger_Rate_UQFF.md',
    'PAPER_014_Primordial_Black_Holes_UQFF_Formation.md',
    'PAPER_014b_EMRI_Aether_Damping_UQFF.md',
}

rows = []
with open(CSV, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fields = reader.fieldnames
    for r in reader:
        if r['filename'] in DONE:
            r['bucket'] = 'H*'
            r['suggested_action'] = 'DONE - G1-G8 + 27-decade ledger block inserted; tailored P-suite hook'
            r['status'] = 'DONE v5.78'
            r['session'] = '262'
            r['commit'] = 'PENDING'  # will set after commit
        rows.append(r)

with open(CSV, 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=fields)
    w.writeheader()
    w.writerows(rows)

print(f'Marked {len(DONE)} rows.')
