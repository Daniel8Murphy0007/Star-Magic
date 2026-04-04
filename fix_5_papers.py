import subprocess, tempfile, os, glob, re

HEADER = 'pdf_header.tex'
BASE_CMD = [
    'pandoc', '--pdf-engine=xelatex',
    '-V', 'geometry:a4paper,top=0.75in,bottom=0.75in,left=0.75in,right=0.75in',
    '-V', 'fontsize=11pt', '-V', 'documentclass=article',
    '-H', HEADER,
    '--pdf-engine-opt=-interaction=nonstopmode',
    '--from=markdown-yaml_metadata_block-raw_tex+smart',
    '--standalone', '--wrap=none',
]
wp_dir = 'whitepapers'
pdf_dir = 'pdf'

def aggressive_clean(text):
    cleaned = []
    for ch in text:
        cp = ord(ch)
        if cp < 128: cleaned.append(ch); continue
        if 0x00B0 <= cp <= 0x00FF: cleaned.append(ch); continue
        if 0x0370 <= cp <= 0x03FF: cleaned.append(ch); continue
        if 0x2200 <= cp <= 0x22FF: cleaned.append(ch); continue
        if ch in '\u2013\u2014': cleaned.append('--'); continue
        if ch in '\u2018\u2019\u201c\u201d': cleaned.append('"'); continue
        if ch == '\u00d7': cleaned.append('x'); continue
        if ch == '\u2248': cleaned.append('~='); continue
        if ch == '\u226a': cleaned.append('<<'); continue
        if ch == '\u2609': cleaned.append('(Sun)'); continue
        if ch == '\u2713': cleaned.append('[OK]'); continue
        if ch == '\u2026': cleaned.append('...'); continue
        if 0x2070 <= cp <= 0x2079: cleaned.append(str(cp - 0x2070)); continue
        if cp == 0x207A: cleaned.append('+'); continue
        if cp == 0x207B: cleaned.append('-'); continue
        if 0x2080 <= cp <= 0x2089: cleaned.append(str(cp - 0x2080)); continue
        if cp == 0x1D62: cleaned.append('i'); continue
        if cp == 0x2099: cleaned.append('n'); continue
        if cp == 0x2081: cleaned.append('1'); continue
        if cp == 0x2082: cleaned.append('2'); continue
        if cp == 0x2083: cleaned.append('3'); continue
        if cp == 0x2084: cleaned.append('4'); continue
        if cp == 0x2212: cleaned.append('-'); continue
        if cp == 0x2264: cleaned.append('<='); continue
        if cp == 0x2265: cleaned.append('>='); continue
        if cp == 0x2260: cleaned.append('!='); continue
    return ''.join(cleaned)

for p in ['072','114','129','137','239']:
    files = glob.glob(f'{wp_dir}/PAPER_{p}_*.md')
    if not files: print(f'PAPER_{p}: not found'); continue
    md = files[0]
    with open(md, encoding='utf-8', errors='replace') as f: text = f.read()
    text2 = aggressive_clean(text)
    tmp = tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', suffix='.md', delete=False)
    tmp.write(text2); tmp.close()
    pdf_name = os.path.basename(md).replace('.md','.pdf')
    pdf_path = os.path.join(pdf_dir, pdf_name)
    r = subprocess.run(BASE_CMD + [tmp.name, '-o', pdf_path], capture_output=True, timeout=120)
    os.unlink(tmp.name)
    if r.returncode == 0 and os.path.exists(pdf_path):
        print(f'PAPER_{p}: OK ({os.path.getsize(pdf_path)//1024}KB)')
    else:
        err = r.stderr.decode('utf-8', errors='replace')
        ekeys = [l for l in err.split('\n') if l.startswith('!')][:3]
        print(f'PAPER_{p}: FAIL {ekeys}')
