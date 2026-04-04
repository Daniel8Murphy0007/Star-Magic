import argparse
import shutil
import os
from pathlib import Path
import re
import subprocess
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare LaTeX paper source for arXiv upload (embodies the guide).")
    parser.add_argument("source_dir", help="Path to your paper source directory")
    parser.add_argument("--tmp_dir", default="arxiv_tmp", help="Name of temporary working directory (default: arxiv_tmp)")
    parser.add_argument("--main_tex", default="main.tex", help="Name of main .tex file (default: main.tex)")
    parser.add_argument("--compile", action="store_true", help="Attempt automatic compilation (requires LaTeX installed)")
    return parser.parse_args()

def copy_to_tmp(source_dir, tmp_dir):
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
    shutil.copytree(source_dir, tmp_dir)
    print(f"✅ Deep copy created: {tmp_dir}/")

def flatten_subdirs(tmp_dir):
    """Flatten common asset subdirectories to root and update paths in all .tex files."""
    common_asset_dirs = ['figures', 'figs', 'img', 'images', 'supplement', 'appendix', 'tables']
    path_replacements = {}
    for dname in common_asset_dirs:
        sub_path = os.path.join(tmp_dir, dname)
        if os.path.exists(sub_path) and os.path.isdir(sub_path):
            for file in os.listdir(sub_path):
                src = os.path.join(sub_path, file)
                dst = os.path.join(tmp_dir, file)
                if os.path.exists(dst):
                    print(f"⚠️  Name conflict: {file} already exists in root. Skipping move.")
                    continue
                shutil.move(src, dst)
                # Record both slashed and double-backslash versions for safety
                path_replacements[f"{dname}/{file}"] = file
                path_replacements[f"{dname}\\\\{file}"] = file
            # Remove empty directory
            if not os.listdir(sub_path):
                os.rmdir(sub_path)
            print(f"✅ Flattened {dname}/ → root")
    if path_replacements:
        update_tex_paths(tmp_dir, path_replacements)
    else:
        print("ℹ️  No common subdirectories found to flatten.")

def update_tex_paths(tmp_dir, replacements):
    """Simple but effective string replacement for \includegraphics, \input, etc."""
    for tex_file in Path(tmp_dir).glob("**/*.tex"):
        with open(tex_file, 'r', encoding='utf-8') as f:
            content = f.read()
        original = content
        for old, new in replacements.items():
            content = content.replace(old, new)
        if content != original:
            with open(tex_file, 'w', encoding='utf-8') as f:
                f.write(content)
            print(f"✅ Updated paths in {tex_file.name}")

def clean_unnecessary_files(tmp_dir):
    """Remove hidden folders, generated files, etc."""
    junk = ['.git', '.gitignore', '__pycache__', '.DS_Store', '.vscode', '.idea', '.svn']
    for item in junk:
        path = os.path.join(tmp_dir, item)
        if os.path.exists(path):
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
            print(f"✅ Removed junk: {item}")
    # Remove common generated LaTeX files (we will regenerate .bbl later)
    for ext in ['.pdf', '.aux', '.log', '.out', '.blg', '.bbl', '.synctex.gz', '.toc', '.lof', '.lot']:
        for f in Path(tmp_dir).glob(f"**/*{ext}"):
            if f.is_file():
                f.unlink()
                print(f"✅ Removed generated file: {f.name}")
    print("⚠️  Review and manually delete any unused .tex / .sty / .cls files now.")

def remove_comments_from_tex(tmp_dir):
    """Strip LaTeX comments (% ...) from every .tex file. Handles escaped \% safely."""
    for tex_file in Path(tmp_dir).glob("**/*.tex"):
        with open(tex_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        cleaned = []
        for line in lines:
            if '%' in line:
                # Split on first unescaped %
                parts = re.split(r'(?<!\\)%', line, 1)
                cleaned_line = parts[0].rstrip() + '\n'
                if cleaned_line.strip():
                    cleaned.append(cleaned_line)
            else:
                cleaned.append(line)
        with open(tex_file, 'w', encoding='utf-8') as f:
            f.writelines(cleaned)
        print(f"✅ Comments stripped from {tex_file.name}")

def add_arxiv_typeout(tmp_dir, main_tex="main.tex"):
    main_path = os.path.join(tmp_dir, main_tex)
    if not os.path.exists(main_path):
        print(f"⚠️  {main_tex} not found – skipping typeout line.")
        return
    with open(main_path, 'r', encoding='utf-8') as f:
        content = f.read()
    if "\\end{document}" in content and "get arXiv to do 4 passes" not in content:
        content = content.rstrip() + "\n\n\\typeout{get arXiv to do 4 passes: Label(s) may have changed. Rerun}\n"
        with open(main_path, 'w', encoding='utf-8') as f:
            f.write(content)
        print("✅ Added arXiv 4-pass typeout line to main.tex")
    else:
        print("ℹ️  Typeout line already present or \\end{document} not found.")

def try_compile(tmp_dir, main_tex="main.tex"):
    """Optional: run pdflatex + bibtex (up to 4 passes). Requires LaTeX in PATH."""
    cwd = os.getcwd()
    os.chdir(tmp_dir)
    try:
        base = main_tex.replace('.tex', '')
        for _ in range(3):
            subprocess.run(["pdflatex", "-interaction=nonstopmode", main_tex], check=False)
        # BibTeX if needed
        if any(Path('.').glob("*.bib")) or os.path.exists(f"{base}.bbl"):
            subprocess.run(["bibtex", base], check=False)
        subprocess.run(["pdflatex", "-interaction=nonstopmode", main_tex], check=False)
        print("✅ Compilation attempt finished. Review any .log or .pdf manually.")
    except FileNotFoundError:
        print("❌ LaTeX (pdflatex/bibtex) not found in PATH – compile manually.")
    finally:
        os.chdir(cwd)

def create_tarball(tmp_dir):
    """Create ax.tar exactly as the guide recommends."""
    cwd = os.getcwd()
    os.chdir(tmp_dir)
    tar_file = "ax.tar"
    subprocess.run(["tar", "-cvvf", tar_file, "."], check=True)
    print(f"✅ Tarball created: {tmp_dir}/{tar_file}")
    os.chdir(cwd)

def main():
    args = parse_args()
    tmp_dir = args.tmp_dir

    print("🚀 arXiv Preparation Script (embodies the full guide)")
    copy_to_tmp(args.source_dir, tmp_dir)

    print("\n" + "="*70)
    print("🚨 MANUAL STEPS REQUIRED BEFORE CONTINUING")
    print("• Merge appendix/supplement if split")
    print("• Edit style file(s) to remove journal-specific text (or ensure published style)")
    print("• Delete any other unused files you see")
    input("Press Enter once the manual steps above are complete...")

    flatten_subdirs(tmp_dir)
    clean_unnecessary_files(tmp_dir)
    remove_comments_from_tex(tmp_dir)
    add_arxiv_typeout(tmp_dir, args.main_tex)

    if args.compile:
        try_compile(tmp_dir, args.main_tex)
        print("\n⚠️  After compilation: keep only the new .bbl and delete .bib + other generated files.")

    create_tarball(tmp_dir)

    print("\n" + "="*70)
    print("✅ Preparation complete!")
    print(f"Upload file →  {tmp_dir}/ax.tar")
    print("Next manual actions:")
    print("• Inspect arXiv’s extracted file list and remove extras")
    print("• Clean abstract/title/authors in the metadata form")
    print("• Choose subject area (consult advisor)")
    print("• After posting, share the paper password with co-authors.")

if __name__ == "__main__":
    main()