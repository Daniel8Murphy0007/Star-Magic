from pathlib import Path
import re
p = Path('Star-MagicProofEngine.py')
text = p.read_text(encoding='utf-8')
start = text.index('self.PROOF_DERIVATION_MODES')
end = text.find('\ndef ', start)
block = text[start:end]
entries = re.findall(r"'([^']+)':\s*\{[^}]*?'callable':\s*([^,\n]+)", block, re.S)
for name, expr in entries:
    if expr.strip().startswith('lambda'):
        print(name, '=>', expr.strip())
print('COUNT', sum(1 for _, expr in entries if expr.strip().startswith('lambda')))
