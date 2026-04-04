import glob
wp_dir = "whitepapers"
fixed = 0
for fpath in glob.glob(wp_dir + "/PAPER_*.md"):
    with open(fpath, encoding="utf-8", errors="replace") as f:
        text = f.read()
    if "\\\\times" not in text and "\\\\bigl" not in text and "\\\\rho" not in text:
        continue
    orig = text
    text = text.replace("\\\\times", "\\times")
    text = text.replace("\\\\bigl", "\\bigl")
    text = text.replace("\\\\bigr", "\\bigr")
    text = text.replace("\\\\rho", "\\rho")
    if text != orig:
        with open(fpath, "w", encoding="utf-8") as f:
            f.write(text)
        fixed += 1
print(f"Fixed: {fixed} files")
