"""Regenerate PDFs for the 10 updated GW papers using pdflatex via pandoc."""
import os, subprocess, sys
sys.stdout.reconfigure(encoding='utf-8')

PANDOC = os.path.expandvars(r'%LOCALAPPDATA%\Pandoc\pandoc.exe')
PDFLATEX = r'C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe'
HEADER = '_pdf_unicode_header.tex'

papers = [
    'PAPER_001_GW170817_UQFF_Damping_Analysis',
    'PAPER_002_GW190425_Mass_Gap_Interpretation',
    'PAPER_003_GW150914_UQFF_vs_LIGO_Strain',
    'PAPER_004_GW170817_BNS_Chirp_Phase_Evolution',
    'PAPER_005_BH_Merger_Energy_Retention_UQFF',
    'PAPER_006_GW170817_Multi_Messenger_Full_Inspiral',
    'PAPER_007_Tidal_Deformability_Constraints_BNS_UQFF',
    'PAPER_008_UQFF_Waveform_Phase_Evolution_Template_Mismatch',
    'PAPER_008b_Full_Inspiral_Waveform_UQFF',
    'PAPER_009_Damping_Mechanism_Decomposition_UQFF',
]

ok = fail = 0
for stem in papers:
    src = os.path.join('whitepapers', stem + '.md')
    dst = os.path.join('pdf', stem + '.pdf')
    cmd = [PANDOC, src, '-o', dst,
           '--pdf-engine=' + PDFLATEX,
           '-H', HEADER,
           '-V', 'geometry:margin=1in',
           '-V', 'fontsize=11pt',
           '-V', 'colorlinks=true',
           '--highlight-style=tango']
    print(f'>>> {stem}')
    r = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
    if r.returncode == 0 and os.path.exists(dst):
        size = os.path.getsize(dst)
        print(f'    OK ({size} bytes)')
        ok += 1
    else:
        print(f'    FAIL rc={r.returncode}')
        # Print last 30 lines of stderr
        err = (r.stderr or '') + (r.stdout or '')
        tail = '\n'.join(err.splitlines()[-30:])
        print(tail)
        fail += 1

print(f'\nDone. OK={ok}  FAIL={fail}')
