with open('index.js', encoding='utf-8', errors='replace') as f:
    content = f.read()
import re

remaining_moji = len(re.findall('ï¿½', content))
remaining_A = content.count('Â²') + content.count('Â³') + content.count('Â¹') + content.count('Ã—')
backtick_n = content.count(';`n ')
print(f'Remaining ï¿½ mojibake: {remaining_moji}')
print(f'Remaining Â-type mojibake: {remaining_A}')
print(f'Backtick-n injections: {backtick_n}')
print(f'File lines: {content.count(chr(10))}')
print(f'File size: {len(content)} bytes')
use_strict = "'use strict'" in content
print(f"'use strict' present: {use_strict}")
r_guard = 'if (!r || r <= 0) r = CONSTANTS.SOLAR_RADIUS;' in content
print(f'r=0 guard in calculateDipMomentumEnergy: {r_guard}')
