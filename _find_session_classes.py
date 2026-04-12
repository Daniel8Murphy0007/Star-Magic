#!/usr/bin/env python3
"""Find all class definitions between _SESSION blocks to populate empty registrations."""
import re

t = open('CondensedPhysics4.py', encoding='utf-8').read()
lines = t.split('\n')

# Find all _SESSION block positions
session_blocks = []
for i, line in enumerate(lines):
    m = re.match(r'^(_SESSION_\w+_CLASSES)\s*=\s*\[', line)
    if m:
        session_blocks.append((i+1, m.group(1)))  # 1-indexed

# Find all class definitions with their line numbers
class_defs = []
for i, line in enumerate(lines):
    m = re.match(r'^class\s+(\w+)', line)
    if m:
        class_defs.append((i+1, m.group(1)))  # 1-indexed

# For each empty session block, find the classes that belong to it
# (classes between the previous session block and this one)
for idx, (block_line, block_name) in enumerate(session_blocks):
    # Check if block is empty
    block_match = re.search(re.escape(block_name) + r'\s*=\s*\[(.*?)\]', t, re.DOTALL)
    if block_match:
        entries = re.findall(r'"(\w+)"', block_match.group(1))
        if entries:
            continue  # already populated
    
    # Find the start of relevant classes (after previous block's last class or after previous block)
    if idx > 0:
        prev_block_line = session_blocks[idx-1][0]
    else:
        prev_block_line = 0
    
    # Find classes between prev_block_line and this block_line
    relevant_classes = [(ln, cn) for ln, cn in class_defs if prev_block_line < ln < block_line]
    
    # But we want classes AFTER the previous session's registration block, not before it
    # Actually the pattern is: classes appear BEFORE their registration block
    # Let's scan backwards from this block to find which classes go here
    
    # Better approach: the header comment says which session they belong to
    # Let's just find classes between the previous registration block end and this one
    
    if relevant_classes:
        print(f"\n{block_name} (line {block_line}) — {len(relevant_classes)} classes to add:")
        for ln, cn in relevant_classes:
            print(f'    "{cn}",')
    else:
        # No classes between prev block and this one — check if classes are AFTER this block
        if idx < len(session_blocks) - 1:
            next_block_line = session_blocks[idx+1][0]
        else:
            next_block_line = len(lines)
        relevant_after = [(ln, cn) for ln, cn in class_defs if block_line < ln < (next_block_line if idx < len(session_blocks)-1 else len(lines)+1)]
        if relevant_after:
            print(f"\n{block_name} (line {block_line}) — {len(relevant_after)} classes AFTER block:")
            for ln, cn in relevant_after:
                print(f'    "{cn}",')
        else:
            print(f"\n{block_name} (line {block_line}) — NO classes found")
