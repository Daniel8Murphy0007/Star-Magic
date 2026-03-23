"""
Extract all 38 module pairs from grok_share_bdfb3a05b06.txt.
Writes .h and .cpp files for each unique class.
Skips files that already exist.
"""
import os
import re

FNAME = 'grok_share_bdfb3a05b06.txt'

with open(FNAME, encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# ─── Find all 'cpp' separator lines (1-based) ───────────────────────────────
sep_lines = []
for i, line in enumerate(lines):
    if line.strip() == 'cpp':
        sep_lines.append(i)   # 0-based index

print(f"Found {len(sep_lines)} 'cpp' separators")

# ─── Pair them up: even index = .h start, odd index = .cpp start ────────────
# Each block pair: h_start, cpp_start, next_h_start (or EOF)
pairs = []
for k in range(0, len(sep_lines) - 1, 2):
    h_sep = sep_lines[k]      # line with 'cpp' before header
    c_sep = sep_lines[k + 1]  # line with 'cpp' before implementation
    # End of .cpp block = start of next h_sep (or end of file)
    if k + 2 < len(sep_lines):
        end_idx = sep_lines[k + 2]
    else:
        end_idx = len(lines)
    pairs.append((h_sep, c_sep, end_idx))

print(f"Total module pairs: {len(pairs)}")

# ─── Helper: infer filename from comment line in block ───────────────────────
def infer_filename(block_lines):
    """Return filename from first // comment, or None."""
    for line in block_lines[:5]:
        m = re.match(r'//\s+([\w\-]+\.(h|cpp|hpp))', line.strip())
        if m:
            return m.group(1)
    return None

def infer_classname(block_lines):
    """Return class name from 'class FooBar {' line."""
    for line in block_lines:
        m = re.match(r'class\s+(\w+)\s*[:{]', line.strip())
        if m:
            return m.group(1)
    return None

# ─── Extract and write each pair ─────────────────────────────────────────────
SKIP_EXISTING = True
created = []
skipped = []
errors = []

# Track class names to handle duplicates (keep later/larger version)
class_registry = {}  # classname -> (pair_index, cpp_size)

# First pass: collect all class names and sizes
for idx, (h_sep, c_sep, end_idx) in enumerate(pairs):
    h_lines = lines[h_sep + 1 : c_sep]
    c_lines = lines[c_sep + 1 : end_idx]
    classname = infer_classname(h_lines) or infer_classname(c_lines)
    if classname:
        size = len(c_lines)
        if classname not in class_registry or size > class_registry[classname][1]:
            class_registry[classname] = (idx, size)

print(f"\nUnique classes found: {len(class_registry)}")
for cn, (idx, sz) in sorted(class_registry.items()):
    print(f"  [{idx:2d}] {cn} ({sz} cpp lines)")

# Second pass: write canonical (largest) version of each class
print("\n--- Writing files ---")
for classname, (idx, size) in sorted(class_registry.items(), key=lambda x: x[1][0]):
    h_sep, c_sep, end_idx = pairs[idx]
    h_lines = lines[h_sep + 1 : c_sep]
    c_lines = lines[c_sep + 1 : end_idx]
    
    h_fname = infer_filename(h_lines)
    c_fname = infer_filename(c_lines)
    
    # Fallback filename from classname
    if not h_fname:
        h_fname = f"{classname}.h"
    if not c_fname:
        c_fname = f"{classname}.cpp"
    
    # Normalize: strip leading whitespace from lines
    h_text = ''.join(h_lines)
    c_text = ''.join(c_lines)
    
    for fname, text in [(h_fname, h_text), (c_fname, c_text)]:
        if SKIP_EXISTING and os.path.exists(fname):
            skipped.append(fname)
            print(f"  SKIP (exists): {fname}")
        else:
            with open(fname, 'w', encoding='utf-8') as f_out:
                f_out.write(text)
            created.append(fname)
            print(f"  CREATED: {fname} ({len(text)} chars)")

print(f"\n=== Summary ===")
print(f"Created: {len(created)} files")
print(f"Skipped (exist): {len(skipped)} files")
print(f"\nCreated files:")
for fn in created:
    print(f"  {fn}")
print(f"\nSkipped files:")
for fn in skipped:
    print(f"  {fn}")
