#!/usr/bin/env python3
"""
Convert SCm .tex files from /pdf to /whitepapers as .md format
Then copy original .tex files to whitepapers for reference
"""
import os
import re
import shutil
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
    header += "---\n\n"
    
    return header + content.strip()

# Paths
pdf_folder = Path("pdf")
whitepapers_folder = Path("whitepapers")

# Files to process (6 SCm files)
files_to_process = [
    "SCm_Holmlid_KER_Validation.tex",
    "SCm_Holmlid_Parkhomov_PonsFleischmann_Upgrade.tex",
    "SCm_Holmlid_Rossi_Parkhomov_Validation.tex",
    "SCm_Mizuno_LENR_Transmutation.tex",
    "SCm_PonsFleischmann_Derivation.tex",
    "SCm_Rossi_ECat_Variants_Unified.tex"
]

print(f"Processing {len(files_to_process)} SCm .tex files\n")

converted = []
failed = []
copied_tex = []

for tex_file in files_to_process:
    tex_path = pdf_folder / tex_file
    
    # Check if source .tex exists
    if not tex_path.exists():
        print(f"[SKIP] {tex_file} (not found in /pdf)")
        failed.append(tex_file)
        continue
    
    try:
        # Read .tex file
        with open(tex_path, 'r', encoding='utf-8') as f:
            tex_content = f.read()
        
        # Convert to markdown
        md_content = convert_tex_to_md(tex_content)
        
        # Create .md filename
        md_name = tex_file.replace('.tex', '.md')
        md_path = whitepapers_folder / md_name
        
        # Write .md file to whitepapers
        with open(md_path, 'w', encoding='utf-8') as f:
            f.write(md_content)
        
        converted.append(md_name)
        print(f"[CONVERT] {tex_file} -> {md_name}")
        
        # Also copy the original .tex file to whitepapers for reference
        tex_dest = whitepapers_folder / tex_file
        shutil.copy2(tex_path, tex_dest)
        copied_tex.append(tex_file)
        print(f"[COPY] {tex_file} -> whitepapers/")
        
    except Exception as e:
        print(f"[ERROR] {tex_file}: {e}")
        failed.append(tex_file)

print(f"\n{'='*70}")
print("CONVERSION SUMMARY")
print(f"{'='*70}")
print(f"[OK] Converted to .md: {len(converted)} files")
print(f"[OK] Copied .tex to whitepapers/: {len(copied_tex)} files")
print(f"[FAILED] {len(failed)} files")

if converted:
    print(f"\nNew .md files created:")
    for name in converted:
        print(f"  {name}")

if copied_tex:
    print(f"\n.tex files copied to whitepapers/:")
    for name in copied_tex:
        print(f"  {name}")

if failed:
    print(f"\nFailed files:")
    for name in failed:
        print(f"  {name}")

print(f"\nNext step: Generate PDFs from the new .md files")
print(f"Run: powershell -ExecutionPolicy Bypass -File generate_scm_pdfs.ps1")
