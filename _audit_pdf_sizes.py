import glob, re, collections, os

pdf_files = sorted(glob.glob('pdf/PAPER_*.pdf'))
sizes = collections.Counter()
other_list = []

pat = re.compile(rb'/MediaBox\s*\[\s*([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s*\]')

for pdf in pdf_files:
    try:
        with open(pdf, 'rb') as f:
            data = f.read()
        m = pat.search(data)
        if m:
            w = round(float(m.group(3)))
            h = round(float(m.group(4)))
            key = '%dx%d' % (w, h)
            sizes[key] += 1
            if key not in ('595x842', '612x792'):
                other_list.append(os.path.basename(pdf))
        else:
            sizes['no_mediabox'] += 1
    except Exception as e:
        sizes['error'] += 1

print('Page size distribution (width x height pts):')
for k, v in sizes.most_common():
    label = ''
    if k == '595x842':
        label = '(A4)'
    elif k == '612x792':
        label = '(Letter)'
    print('  %s  %s -> %d PDFs' % (k, label, v))

if other_list:
    print('Non-standard (%d): %s' % (len(other_list), other_list[:10]))
print('Total checked: %d' % sum(sizes.values()))
