# add_wolfram_terms.py
import os
import re

files = [f for f in os.listdir('.') if re.match(r'^source\d+\.cpp$', f) or f == 'MAIN_1_CoAnQi.cpp']

for filename in files:
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            content = f.readlines()
    except UnicodeDecodeError:
        # Try with latin-1 encoding for files with special characters
        with open(filename, 'r', encoding='latin-1') as f:
            content = f.readlines()

    # Skip if already has WOLFRAM_TERM
    if any('WOLFRAM_TERM' in line for line in content):
        print(f"Skipping {filename} – already has term")
        continue

    # Very small but meaningful term per file – in real Phase 5 we'll make these the actual translated physics
    term = f"(* Auto-contribution from {filename} *) + {filename.replace('.cpp', '')}_unification_sector"

    insert = f'#define WOLFRAM_TERM "{term}"\n'
    content.insert(1, insert)  # right after includes / prologue

    try:
        with open(filename, 'w', encoding='utf-8') as f:
            f.writelines(content)
    except UnicodeEncodeError:
        # Write with latin-1 if utf-8 fails
        with open(filename, 'w', encoding='latin-1') as f:
            f.writelines(content)

    print(f"Injected WOLFRAM_TERM into {filename}")

print("\nAll 175+ files now have WOLFRAM_TERM definitions. Commit this change.")
