#!/usr/bin/env python3
"""Convert PAPER_*.tex files to PAPER_*.md format"""
import os
import re
from pathlib import Path

def convert_tex_to_md(tex_content):
    """Convert LaTeX content to Markdown"""
    
    # Extract title, author, date from \title, \author, \date
    title_match = re.search(r'\\title\{([^}]+)\}', tex_content)
    author_match = re.search(r'\\author\{([^}]+)\}', tex_content)
    date_match = re.search(r'\\date\{([^}]+)\}', tex_content)
    
    title = title_match.group(1) if title_match else "Untitled"
    author = author_match.group(1) if author_match else "Star-Magic UQFF Program"
    date = date_match.group(1) if date_match else "May 2026"
    
    # Remove everything before \begin{document}
    content = re.sub(r'^.*?\\begin\{document\}', '', tex_content, flags=re.DOTALL)
    
    # Remove \end{document} and everything after
    content = re.sub(r'\\end\{document\}.*$', '', content, flags=re.DOTALL)
    
    # Remove \maketitle
    content = re.sub(r'\\maketitle\s*', '', content)
    
    # Convert \section*{...} and \section{...} to ##
    content = re.sub(r'\\section\*?\{([^}]+)\}', r'## \1', content)
    
    # Convert \subsection{...} to ###
    content = re.sub(r'\\subsection\*?\{([^}]+)\}', r'### \1', content)
    
    # Convert \subsubsection{...} to ####
    content = re.sub(r'\\subsubsection\*?\{([^}]+)\}', r'#### \1', content)
    
    # Convert \textbf{...} to **...**
    content = re.sub(r'\\textbf\{([^}]+)\}', r'**\1**', content)
    
    # Convert \textit{...} to *...*
    content = re.sub(r'\\textit\{([^}]+)\}', r'*\1*', content)
    
    # Convert \texttt{...} to `...`
    content = re.sub(r'\\texttt\{([^}]+)\}', r'`\1`', content)
    
    # Convert \[ ... \] to $$ ... $$
    content = re.sub(r'\\\[(.*?)\\\]', r'$$\1$$', content, flags=re.DOTALL)
    
    # Remove remaining LaTeX commands that are simple
    content = re.sub(r'\\(?:noindent|hspace\{[^}]*\}|vspace\{[^}]*\}|\\)\s*', '', content)
    
    # Convert \newcommand{...}{...} definitions to comments
    content = re.sub(r'\\newcommand.*?\}.*?\}', '', content, flags=re.DOTALL)
    
    # Handle \href{url}{text}
    content = re.sub(r'\\href\{([^}]+)\}\{([^}]+)\}', r'[\2](\1)', content)
    
    # Handle line breaks
    content = re.sub(r'\\\\\s*\n', '\n\n', content)
    
    # Remove remaining backslash escapes for special characters
    content = re.sub(r'\\([_&%$#{}])', r'\1', content)
    
    # Clean up multiple blank lines
    content = re.sub(r'\n\n\n+', '\n\n', content)
    
    # Create markdown header
    header = f"# {title}\n\n"
    header += f"**Author:** {author}\n\n"
    header += f"**Date:** {date}\n\n"
    
    # Add a horizontal rule
    header += "---\n\n"
    
    return header + content.strip()

def get_paper_info(filename):
    """Extract paper number and extract a clean title"""
    # Match PAPER_1234_Title_Words_Here.tex
    match = re.match(r'PAPER_(\d+)_(.+)\.tex', filename)
    if match:
        paper_num = match.group(1)
        title_text = match.group(2).replace('_', ' ')
        return paper_num, title_text
    return None, None

# Process all .tex files in whitepapers folder
whitepapers_dir = Path('whitepapers')
tex_files = sorted(whitepapers_dir.glob('PAPER_1*.tex'))

print(f"Found {len(tex_files)} .tex files to convert\n")

converted = []
skipped = []

for tex_file in tex_files:
    # Skip if .md already exists
    md_name = tex_file.stem + '.md'
    md_path = whitepapers_dir / md_name
    
    if md_path.exists():
        skipped.append(tex_file.name)
        continue
    
    try:
        # Read .tex file
        with open(tex_file, 'r', encoding='utf-8') as f:
            tex_content = f.read()
        
        # Convert to markdown
        md_content = convert_tex_to_md(tex_content)
        
        # Write .md file
        with open(md_path, 'w', encoding='utf-8') as f:
            f.write(md_content)
        
        converted.append((tex_file.name, md_name))
        print(f"✓ {tex_file.name} → {md_name}")
    
    except Exception as e:
        print(f"✗ {tex_file.name}: {e}")

print(f"\n{'='*70}")
print(f"CONVERSION SUMMARY")
print(f"{'='*70}")
print(f"Converted: {len(converted)} files")
print(f"Skipped (already exists): {len(skipped)} files")

if converted:
    print(f"\nFiles created:")
    for tex_name, md_name in converted:
        print(f"  {md_name}")

if skipped:
    print(f"\nAlready existing (skipped):")
    for name in skipped[:5]:
        print(f"  {name}")
    if len(skipped) > 5:
        print(f"  ... and {len(skipped)-5} more")

print(f"\nTotal new .md files: {len(converted)}")
