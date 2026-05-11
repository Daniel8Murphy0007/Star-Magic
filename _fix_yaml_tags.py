"""Fix YAML tag arrays containing unescaped { } in whitepaper frontmatter.
Replaces tokens like F_{U\_Bi\_i} with F_U_Bi_i (safe YAML scalar)."""
import re, glob, pathlib

YAML_TAG_RE = re.compile(r'^tags:\s*\[(.*)\]\s*$', re.MULTILINE)

def clean_token(tok: str) -> str:
    t = tok.strip().strip("'\"")
    # Remove braces and backslash-underscore escapes
    t = t.replace('\\_', '_').replace('{', '').replace('}', '')
    return t

count = 0
for fp in glob.glob('whitepapers/PAPER_*.md'):
    p = pathlib.Path(fp)
    text = p.read_text(encoding='utf-8')
    # Only first 30 lines (frontmatter)
    head_end = 0
    lines = text.split('\n')
    # find frontmatter end
    if lines and lines[0].strip() == '---':
        for i in range(1, min(40, len(lines))):
            if lines[i].strip() == '---':
                head_end = i
                break
    if head_end == 0:
        continue
    head = '\n'.join(lines[:head_end+1])
    rest = '\n'.join(lines[head_end+1:])
    m = YAML_TAG_RE.search(head)
    if not m:
        continue
    inside = m.group(1)
    if '{' not in inside and '\\_' not in inside:
        continue
    toks = [clean_token(t) for t in inside.split(',')]
    new_line = 'tags: [' + ', '.join(f'"{t}"' for t in toks) + ']'
    new_head = head[:m.start()] + new_line + head[m.end():]
    p.write_text(new_head + '\n' + rest, encoding='utf-8')
    count += 1
    print(f'  fixed {p.name}')
print(f'Updated {count} files.')
