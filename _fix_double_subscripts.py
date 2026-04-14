"""Fix double-subscript errors in LaTeX math environments.
Patterns like F_U_Bi_i inside $$ blocks cause 'Double subscript' LaTeX errors.
Fix: wrap multi-underscore identifiers in \text{} inside math environments."""
import os, re

def fix_math_double_subscripts(content):
    """Fix double/triple/quad subscript identifiers inside math environments."""
    
    # Pattern to match multi-underscore identifiers (2+ underscores)
    # like F_U_Bi_i, g_compressed, rho_SCm, etc.
    multi_underscore = re.compile(
        r'(?<![\\])([A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+){2,})'
    )
    
    def fix_in_math(match):
        """Replace multi-underscore identifiers with \\text{escaped} version."""
        ident = match.group(1)
        # Don't touch things that look like LaTeX commands
        if ident.startswith('\\'):
            return match.group(0)
        # Escape underscores for \text{} context
        escaped = ident.replace('_', r'\_')
        return r'\text{' + escaped + '}'
    
    # Process display math blocks ($$...$$)
    def fix_display_math(m):
        block = m.group(0)
        prefix = m.group(1)  # opening $$
        inner = m.group(2)
        suffix = m.group(3)  # closing $$
        fixed_inner = multi_underscore.sub(fix_in_math, inner)
        return prefix + fixed_inner + suffix
    
    # Process inline math ($...$) - but not $$
    def fix_inline_math(m):
        block = m.group(0)
        if block.startswith('$$'):
            return block  # skip display math, handled separately
        inner = m.group(1)
        fixed_inner = multi_underscore.sub(fix_in_math, inner)
        return '$' + fixed_inner + '$'
    
    # Fix display math first ($$...$$)
    content = re.sub(r'(\$\$)(.*?)(\$\$)', fix_display_math, content, flags=re.DOTALL)
    
    # Fix inline math ($...$) - careful not to match $$
    content = re.sub(r'(?<!\$)\$(?!\$)(.*?)(?<!\$)\$(?!\$)', fix_inline_math, content)
    
    return content


def fix_table_double_subscripts(content):
    """Fix double subscript identifiers in markdown table cells.
    Pandoc can convert table content to math mode in some cases."""
    lines = content.split('\n')
    fixed_lines = []
    
    for line in lines:
        if '|' in line and '_' in line:
            # In table cells, escape multiple underscores
            parts = line.split('|')
            new_parts = []
            for part in parts:
                # Count underscores in potential identifiers
                # Only fix if it looks like a code identifier (not valid LaTeX)
                part = re.sub(
                    r'(?<![`\\$])([A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+){2,})(?![`$])',
                    lambda m: '`' + m.group(1) + '`',
                    part
                )
                new_parts.append(part)
            fixed_lines.append('|'.join(new_parts))
        else:
            fixed_lines.append(line)
    
    return '\n'.join(fixed_lines)


fixed_count = 0
for fn in sorted(os.listdir('whitepapers')):
    if not fn.endswith('.md'):
        continue
    path = os.path.join('whitepapers', fn)
    with open(path, encoding='utf-8') as f:
        content = f.read()
    
    original = content
    content = fix_math_double_subscripts(content)
    content = fix_table_double_subscripts(content)
    
    if content != original:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        fixed_count += 1

print(f"Fixed double subscripts in {fixed_count} papers")
