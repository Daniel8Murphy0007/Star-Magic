"""Fix \\rightr<delim> -> \\right<delim> and \\leftl<delim> -> \\left<delim>"""
import re, glob, pathlib
PAT_R = re.compile(r'\\rightr(floor|ceil|angle|brace|bracket|paren|vert|Vert)\b')
PAT_L = re.compile(r'\\leftl(floor|ceil|angle|brace|bracket|paren|vert|Vert)\b')
count = 0
for fp in glob.glob('whitepapers/PAPER_*.md'):
    p = pathlib.Path(fp)
    s = p.read_text(encoding='utf-8')
    new = PAT_R.sub(r'\\right\\\1', s)
    new = PAT_L.sub(r'\\left\\\1', new)
    if new != s:
        p.write_text(new, encoding='utf-8')
        count += 1
        print(f'  fixed {p.name}')
print(f'Updated {count} files.')
