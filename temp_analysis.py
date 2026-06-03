from pathlib import Path
import re
from pprint import pprint
p = Path('Star-MagicProofEngine.py')
text = p.read_text(encoding='utf-8')
start = text.index('self.PROOF_DERIVATION_MODES')
# find end of dict by matching braces? simple heuristic: find the first '}' at column 8 after start where next line starts with 'def '
end = text.find('\ndef ', start)
block = text[start:end]
entries = re.findall(r"'([^']+)':\s*\{[^}]*?'callable':\s*([^,\n]+)", block, re.S)
print('entries count', len(entries))
for name, callable_expr in entries:
    if callable_expr.strip().startswith('lambda'):
        print('LAMBDA', name, callable_expr.strip())
    elif callable_expr.strip().startswith('self.'):
        print('METHOD', name, callable_expr.strip())
    else:
        print('OTHER', name, callable_expr.strip())
# find methods with simple return string bodies
meths = re.findall(r'def\s+([A-Za-z0-9_]+)\s*\([^\)]*\):\n((?:\s+[^\n]+\n)+)', text)
for name, body in meths:
    body_strip = body.strip()
    if body_strip.startswith('return') and not ('math.' in body or 'np.' in body or 'exp(' in body or 'cos(' in body or 'sin(' in body or 'pow(' in body or '*' in body or '+' in body or '-' in body or '/' in body):
        print('POSSIBLE-STUB', name)
