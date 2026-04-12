import os, re
pdf_nums = set()
for f in os.listdir('pdf'):
    if f.startswith('PAPER_') and f.endswith('.pdf'):
        m = re.match(r'PAPER_(\d+[a-z]?)', f)
        if m: pdf_nums.add(m.group(1))
wp_nums = set()
for f in os.listdir('whitepapers'):
    if f.startswith('PAPER_') and f.endswith('.md'):
        m = re.match(r'PAPER_(\d+[a-z]?)', f)
        if m: wp_nums.add(m.group(1))
pdf_ct = len([f for f in os.listdir('pdf') if f.startswith('PAPER_') and f.endswith('.pdf')])
wp_ct = len([f for f in os.listdir('whitepapers') if f.startswith('PAPER_') and f.endswith('.md')])
print(f"PDFs: {pdf_ct} files, {len(pdf_nums)} unique numbers")
print(f"WPs:  {wp_ct} files, {len(wp_nums)} unique numbers")
print(f"Orphan PDFs: {sorted(pdf_nums - wp_nums)}")
print(f"Orphan WPs:  {sorted(wp_nums - pdf_nums)}")
status = "PERFECT" if not (pdf_nums - wp_nums) and not (wp_nums - pdf_nums) else "MISMATCH"
print(f"Parity: {status}")
