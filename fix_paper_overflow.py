#!/usr/bin/env python3
"""Fix text overflow in MUGE whitepaper markdown files for PDF generation.

Breaks long display equations into \\begin{aligned} blocks and wraps
long text paragraphs to fit within PDF margins (A4 with 0.75in margins).
"""
import re
import sys
import os

def wrap_text_line(line, max_width=80):
    """Wrap a long text line at word boundaries."""
    if len(line) <= max_width:
        return line
    # Don't wrap lines that are part of aligned blocks or tables
    if line.strip().startswith(('|', '>', '<!--', '```', '\\', '&')):
        return line
    # Don't wrap lines that are pure display math
    if line.strip().startswith('$$'):
        return line
    
    words = line.split(' ')
    lines = []
    current = ''
    for word in words:
        if current and len(current) + 1 + len(word) > max_width:
            lines.append(current)
            current = word
        else:
            current = current + ' ' + word if current else word
    if current:
        lines.append(current)
    return '\n'.join(lines)


def fix_display_equation(match):
    """Break a long display equation into aligned block."""
    content = match.group(1).strip()
    
    # Skip if already an aligned block
    if '\\begin{aligned}' in content:
        return match.group(0)
    
    # Skip short equations
    if len(content) < 100:
        return match.group(0)
    
    # Handle boxed equations
    is_boxed = content.startswith('\\boxed{') and content.endswith('}')
    if is_boxed:
        inner = content[7:-1]
        if '\\begin{aligned}' in inner:
            return match.group(0)
        # Try to break at = sign
        if '=' in inner:
            parts = inner.split('=', 1)
            lhs = parts[0].strip()
            rhs = parts[1].strip()
            # Break RHS at + signs if long
            rhs_parts = split_at_operators(rhs)
            if len(rhs_parts) > 1:
                aligned = f"\\boxed{{\\begin{{aligned}}\n{lhs} &= {rhs_parts[0]} \\\\\n"
                for p in rhs_parts[1:]:
                    aligned += f"&\\quad {p} \\\\\n"
                # Remove last \\
                aligned = aligned.rstrip('\\\\\n') + '\n'
                aligned += "\\end{aligned}}"
                return f"$${aligned}$$"
        return match.group(0)
    
    # Regular equation with = sign
    if '=' in content and not content.startswith('\\begin'):
        parts = content.split('=', 1)
        lhs = parts[0].strip()
        rhs = parts[1].strip()
        
        # Break RHS at + or major operators
        rhs_parts = split_at_operators(rhs)
        if len(rhs_parts) > 1:
            aligned = f"\\begin{{aligned}}\n{lhs} &= {rhs_parts[0]} \\\\\n"
            for p in rhs_parts[1:]:
                aligned += f"&\\quad {p} \\\\\n"
            aligned = aligned.rstrip('\\\\\n') + '\n'
            aligned += "\\end{aligned}"
            return f"$${aligned}$$"
    
    # Break at + signs for sums
    if ' + ' in content and len(content) > 110:
        parts = split_at_operators(content)
        if len(parts) > 1:
            aligned = f"\\begin{{aligned}}\n& {parts[0]} \\\\\n"
            for p in parts[1:]:
                aligned += f"&\\quad {p} \\\\\n"
            aligned = aligned.rstrip('\\\\\n') + '\n'
            aligned += "\\end{aligned}"
            return f"$${aligned}$$"
    
    return match.group(0)


def split_at_operators(expr):
    """Split expression at + or major grouping points."""
    # Split at top-level + signs (not inside braces)
    parts = []
    depth = 0
    current = ''
    i = 0
    while i < len(expr):
        ch = expr[i]
        if ch in ('{', '(', '['):
            depth += 1
            current += ch
        elif ch in ('}', ')', ']'):
            depth -= 1
            current += ch
        elif depth == 0 and i > 0 and expr[i:i+3] == ' + ':
            if len(current.strip()) > 20:  # Don't split tiny fragments
                parts.append(current.strip())
                current = '+ '
                i += 3
                continue
            else:
                current += ch
        elif depth == 0 and i > 0 and expr[i:i+5] == ' \\;+':
            parts.append(current.strip())
            current = '+ '
            i += 5
            # skip trailing whitespace
            while i < len(expr) and expr[i] == ' ':
                i += 1
            continue
        else:
            current += ch
        i += 1
    if current.strip():
        parts.append(current.strip())
    
    # If only 1 part, try splitting at \cdot
    if len(parts) <= 1 and ' \\cdot ' in expr and len(expr) > 110:
        parts = []
        tokens = expr.split(' \\cdot ')
        current = tokens[0]
        for t in tokens[1:]:
            if len(current) + len(t) + 7 > 80:
                parts.append(current.strip())
                current = '\\cdot ' + t
            else:
                current += ' \\cdot ' + t
        if current.strip():
            parts.append(current.strip())
    
    return parts if len(parts) > 1 else [expr]


def fix_file(filepath):
    """Fix all overflow issues in a markdown file."""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    original = content
    
    # 1. Fix display equations ($$...$$)
    # Match single-line display equations
    content = re.sub(
        r'\$\$((?:(?!\$\$).)+)\$\$',
        fix_display_equation,
        content,
        flags=re.DOTALL
    )
    
    # 2. Wrap long text paragraphs
    lines = content.split('\n')
    new_lines = []
    for line in lines:
        if len(line) > 110:
            # Skip table rows, comments, blockquotes, code, math
            if line.strip().startswith(('|', '>', '<!--', '```', '$$', '\\', '&')):
                new_lines.append(line)
            else:
                new_lines.append(wrap_text_line(line, 80))
        else:
            new_lines.append(line)
    content = '\n'.join(new_lines)
    
    if content != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False


if __name__ == '__main__':
    papers = sys.argv[1:] if len(sys.argv) > 1 else []
    if not papers:
        print("Usage: python fix_paper_overflow.py <paper_file> [paper_file2] ...")
        sys.exit(1)
    
    for paper in papers:
        path = os.path.join('whitepapers', paper) if not paper.startswith('whitepapers') else paper
        if not os.path.exists(path):
            print(f"  SKIP {paper} (not found)")
            continue
        changed = fix_file(path)
        print(f"  {'FIXED' if changed else 'OK   '} {paper}")
