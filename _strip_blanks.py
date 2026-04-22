"""Strip spurious blank lines: blank between def/class header and its docstring."""
import ast

with open('CondensedPhysics.py', encoding='utf-8-sig', errors='ignore') as f:
    lines = f.readlines()

before = len(lines)

result = []
i = 0
while i < len(lines):
    line = lines[i]
    if line.strip() == '':
        prev = result[-1].strip() if result else ''
        # Find next non-blank
        j = i + 1
        while j < len(lines) and lines[j].strip() == '':
            j += 1
        nxt = lines[j].strip() if j < len(lines) else ''

        # Remove blank between def/class header and docstring
        if prev.endswith(':') and (nxt.startswith('"""') or nxt.startswith("'''")):
            i += 1
            continue
    result.append(line)
    i += 1

content = ''.join(result)
after = len(result)
blank_after = sum(1 for l in result if l.strip() == '')
print(f'Before: {before:,}  After: {after:,}  Blank: {blank_after:,}  Removed: {before-after:,}')

try:
    ast.parse(content)
    print('Syntax: OK')
    with open('CondensedPhysics.py', 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)
    print('WRITTEN.')
except SyntaxError as e:
    print(f'SyntaxError at line {e.lineno}: {e.msg}')
    print('NOT written.')
