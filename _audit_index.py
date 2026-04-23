import re

with open('index.js', encoding='utf-8', errors='replace') as f:
    content = f.read()
lines = content.split('\n')

# Backtick-n injections (PowerShell newline corruption)
bad = [(i+1, line) for i, line in enumerate(lines) if ';`n ' in line or ';`n\t' in line]
print(f'Backtick-n injections: {len(bad)}')
for lineno, line in bad:
    print(f'  {lineno}: {line[:110]}')

# Mojibake
moji = re.findall(r'ï¿½', content)
print(f'\nMojibake (ï¿½) count: {len(moji)}')

# Duplicate exports block
early_exports = len(re.findall(r'module\.exports\s*=\s*module\.exports\s*\|\|', content))
print(f'\nEarly module.exports||{{}} assignments: {early_exports}')

# r==0 division risks (functions that divide by r without guard)
funcs_with_r_div = re.findall(r'function \w+\([^)]*\)[^{]*\{[^}]*/ r\b', content)
print(f'Top-level funcs dividing by r (sample): {len(funcs_with_r_div)}')

# console.log style (emoji chars as unicode escapes vs literal)
print(f'Emoji in console.log: {len(re.findall(chr(127), content))}')

# undefined variable references (SN1987A etc referenced in exports but defined?)
exported_names = re.findall(r'module\.exports\.(\w+)\s*=\s*(\w+)', content)
defined_classes = set(re.findall(r'^class (\w+)', content, re.MULTILINE))
defined_funcs = set(re.findall(r'^function (\w+)', content, re.MULTILINE))
defined_consts = set(re.findall(r'^const (\w+)\s*=', content, re.MULTILINE))
defined_all = defined_classes | defined_funcs | defined_consts
missing = [(k, v) for k, v in exported_names if v not in defined_all]
print(f'\nExported symbols NOT defined as class/function/const at top level: {len(missing)}')
for k, v in missing[:20]:
    print(f'  module.exports.{k} = {v}  <-- UNDEFINED')

print(f'\nTotal classes: {len(defined_classes)}')
print(f'Total functions: {len(defined_funcs)}')
print(f'Total top-level consts: {len(defined_consts)}')
