import os
files=['PAPER_010_Post_Merger_Oscillations_Remnant_Mass_UQFF','PAPER_010b_Time_Domain_Chirp_23Hz_UQFF','PAPER_011_Stochastic_GW_Background_UQFF_Implications','PAPER_011b_Amplitude_Reduction_Factor_UQFF','PAPER_012_Eccentric_Binary_Circularization_UQFF','PAPER_012b_GW150914_Waveform_Validation','PAPER_013_Magnetar_Spin_Down_UQFF_Framework','PAPER_013b_LISA_SMBH_Merger_Rate_UQFF','PAPER_014_Primordial_Black_Holes_UQFF_Formation','PAPER_014b_EMRI_Aether_Damping_UQFF']
for f in files:
    md='whitepapers/'+f+'.md'
    pdf='pdf/'+f+'.pdf'
    txt=open(md,encoding='utf-8').read()
    has='v5.78 Closure' in txt
    sz=os.path.getsize(pdf) if os.path.exists(pdf) else 0
    tag='YES' if has else 'NO'
    print(f'{f:60s} block={tag} pdf={sz//1024}KB')
