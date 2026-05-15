import csv
DONE = {
    'PAPER_009b_Aether_String_TRZ_Damping_GW.md',
    'PAPER_015_Cosmological_Implications_UQFF_Modified_GW_Propagation.md',
    'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF.md',
    'PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations.md',
    'PAPER_016b_White_Dwarf_Foreground_UQFF.md',
    'PAPER_017_Redshift_Corrections_z1_in_UQFF_GW_Propagation.md',
    'PAPER_018_Aether_Noise_Spectrum_Characterization_for_LISA.md',
    'PAPER_019_Pulsar_Timing_Array_Anomalies_UQFF.md',
    'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md',
    'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density.md',
}
CSV='_job_b_categorization.csv'
with open(CSV,'r',encoding='utf-8') as f:
    rows=list(csv.DictReader(f))
    fields=list(rows[0].keys())
for r in rows:
    if r['filename'] in DONE:
        r['bucket']='H*' if r['bucket']=='H' else r['bucket']  # keep A/D bucket info
        if r['bucket'] in ('A','D'):
            r['bucket']=r['bucket']+'*'
        r['suggested_action']='DONE - G1-G8 + 27-decade ledger block inserted; tailored P-suite/ledger/xi hook'
        r['status']='DONE v5.78'
        r['session']='262'
        r['commit']='PENDING'
with open(CSV,'w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
print(f'Marked {len(DONE)} rows.')
