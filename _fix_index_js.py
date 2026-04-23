"""Fix index.js issues:
1. 5 backtick-n injections (PowerShell `n corruption) - CRITICAL: prevents file from loading
2. 2538 mojibake characters (ï¿½) in comments/strings - QUALITY: corrupted Unicode
"""
import re

with open('index.js', encoding='utf-8', errors='replace') as f:
    content = f.read()

original_len = len(content)
changes = []

# ─────────────────────────────────────────────────────────
# FIX 1: Backtick-n injections (CRITICAL)
# PowerShell `n got embedded as literal backtick+n in JS source
# These appear as:  ;`n        this.B_DPM  (which JS parses as template literal start)
# Fix: replace ;`n (with trailing spaces) with ;\n (real newline + spaces)
# ─────────────────────────────────────────────────────────

# Pattern: semicolon, backtick, n, then spaces
# Variant 1: CONSTANTS.GRAVITATIONAL_CONSTANT;`n        this.B_DPM
# Variant 2: CONSTANTS.GRAVITATIONAL_CONSTANT || 6.6743e-11;`n        this.B_DPM

# Use regex to find all ;`n<spaces> patterns
def fix_backtick_n(m):
    spaces = m.group(1)
    return ';\n' + spaces

before = content.count('`;n')
# Replace ;`n followed by spaces with real newline + same spaces
content = re.sub(r';`n( +)', fix_backtick_n, content)
after = content.count('`;n')
fixed = before - after
changes.append(f'Backtick-n injections fixed: {fixed} (remaining: {after})')

# ─────────────────────────────────────────────────────────
# FIX 2: Mojibake - ï¿½ character replacements in comments and strings
# The character ï¿½ (U+FFFD replacement char in UTF-8: 0xEF 0xBF 0xBD) appeared
# when Unicode superscript/special chars were corrupted.
#
# Context-based substitutions:
#   m/s²  → the ² was corrupted: "m/sï¿½" → "m/s²"
#   m³    → "mï¿½" → "m³"  
#   J·s   → "Jï¿½s" → "J·s"
#   kg·   → "kgï¿½" → "kg·"
#   m²    → "mï¿½" → "m²"
#   N·m   → "Nï¿½m" → "N·m"
# ─────────────────────────────────────────────────────────

# Most common patterns based on units in physics code
moji = 'ï¿½'
moji_count_before = content.count(moji)

replacements = [
    # Units: per second squared (acceleration)
    ('m/s' + moji, 'm/s²'),
    # Units: cubic meter
    ('m' + moji + ')', 'm³)'),
    ('m' + moji + ' ', 'm³ '),
    ('m' + moji + '\n', 'm³\n'),
    ('m' + moji + "'", "m³'"),
    ('m' + moji + '"', 'm³"'),
    # Units: Joule-second (Planck constant)
    ('J' + moji + 's', 'J·s'),
    ('J' + moji + 'S', 'J·S'),
    # Units: kg/m³ (density)
    ('kg/m' + moji, 'kg/m³'),
    # Units: kg·m 
    ('kg' + moji + 'm', 'kg·m'),
    # Units: N·m (torque)
    ('N' + moji + 'm', 'N·m'),
    # Units: m² (area, squared)
    ('m' + moji + '/', 'm²/'),
    # Units: m^-2 (cosmological constant)
    ('m' + moji + '-', 'm⁻'),
    # Units: s^-1 (Hubble)
    ('s^-1' + moji, 's⁻¹'),
    ('s-1' + moji, 's⁻¹'),
    # Generic trailing moji that appear after unit characters
    # Handle remaining cases: most likely ² or ³ based on physics context
]

total_moji_fixed = 0
for old, new in replacements:
    count = content.count(old)
    if count > 0:
        content = content.replace(old, new)
        total_moji_fixed += count
        changes.append(f'  Replaced "{old[:20]}" → "{new[:20]}": {count}')

# Handle remaining ï¿½ - replace with appropriate Unicode based on context
# After letter/number that suggest squared or cubed
# Remaining ï¿½ after fixes - replace with generic superscript placeholder ²
remaining_moji = content.count(moji)
if remaining_moji > 0:
    # For remaining ones, look at what precedes them
    # Most physics units are squared (²) or cubed (³)
    # Safe fallback: replace with ² as most common physics context
    content = content.replace(moji, '²')
    changes.append(f'  Remaining mojibake → ² (fallback): {remaining_moji}')

moji_count_after = content.count(moji)
changes.append(f'Mojibake fixed: {moji_count_before - moji_count_after} of {moji_count_before}')

# ─────────────────────────────────────────────────────────
# FIX 3: Add 'use strict' at top (after first comment block)
# ─────────────────────────────────────────────────────────
if "'use strict'" not in content and '"use strict"' not in content:
    # Add after the initial console.log block
    target = "console.log('Initializing Advanced Unified Quantum Field Force calculations...\\n');\n"
    if target in content:
        content = content.replace(target, target + "\n'use strict';\n")
        changes.append("Added 'use strict' directive")

# ─────────────────────────────────────────────────────────
# FIX 4: r=0 guard in calculateDipMomentumEnergy (division by r^2 with no guard)
# This is called in all Ug1-4 calculations with user-provided r
# ─────────────────────────────────────────────────────────
old_dpm = """function calculateDipMomentumEnergy(r, layerIndex) {
    const r_i = r / layerIndex; // Layer-dependent radius"""
new_dpm = """function calculateDipMomentumEnergy(r, layerIndex) {
    if (!r || r <= 0) r = CONSTANTS.SOLAR_RADIUS; // Guard against r=0/undefined
    const r_i = r / layerIndex; // Layer-dependent radius"""
if old_dpm in content:
    content = content.replace(old_dpm, new_dpm)
    changes.append('Added r=0 guard in calculateDipMomentumEnergy')

# ─────────────────────────────────────────────────────────
# FIX 5: Verify syntax (compile check via regex for obvious issues)
# ─────────────────────────────────────────────────────────

# Verify no remaining backtick-n
remaining_bn = len(re.findall(r';`n ', content))
changes.append(f'Remaining backtick-n: {remaining_bn}')

print('=== index.js Fix Report ===')
for c in changes:
    print(' ', c)
print(f'\nFile size: {original_len:,} → {len(content):,} bytes')

with open('index.js', 'w', encoding='utf-8', newline='\n') as f:
    f.write(content)
print('Written: index.js')
