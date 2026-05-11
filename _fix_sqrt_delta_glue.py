"""Bulk: fix \sqrt{}(...) -> \sqrt{...}; \Deltax->\Delta x; \Deltap->\Delta p; \omegat->\omega t."""
import pathlib, re

PATS = [
    (re.compile(r'\\sqrt\{\}\(([^()]+)\)'), r'\\sqrt{\1}'),
    (re.compile(r'\\Deltax(?![A-Za-z])'), r'\\Delta x'),
    (re.compile(r'\\Deltap(?![A-Za-z])'), r'\\Delta p'),
    (re.compile(r'\\omegat(?![A-Za-z])'), r'\\omega t'),
    (re.compile(r'\\proptoomega(?![A-Za-z])'), r'\\propto \\omega'),
    (re.compile(r'\\pit(?![A-Za-z])'), r'\\pi t'),
    (re.compile(r'\\Deltat(?![A-Za-z])'), r'\\Delta t'),
    (re.compile(r'\\deltaa(?![A-Za-z])'), r'\\delta a'),
    (re.compile(r'\\DeltaBR(?![A-Za-z])'), r'\\Delta \\mathrm{BR}'),
    (re.compile(r'\\gammat(?![A-Za-z])'), r'\\gamma t'),
    (re.compile(r'\\cdotT(?![A-Za-z])'), r'\\cdot T'),
    (re.compile(r'\\cdotcos(?![A-Za-z])'), r'\\cdot \\cos'),
]

count = 0
for p in pathlib.Path('whitepapers').glob('*.md'):
    s = p.read_text(encoding='utf-8')
    o = s
    for r, rep in PATS:
        s = r.sub(rep, s)
    if s != o:
        p.write_text(s, encoding='utf-8')
        count += 1
        print('fixed', p.name)
print('total', count)
