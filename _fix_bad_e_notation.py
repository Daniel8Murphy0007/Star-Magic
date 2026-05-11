"""Fix mangled scientific notation: 'X.X\1\times 10^{-1}Y' -> 'X.X\times 10^{-1Y}'.

Pattern signature is `\1\times 10^{-1}` followed by a single digit.
"""
import re, pathlib

PAT = re.compile(r"\\1\\times 10\^\{-1\}(\d)")

def fix_text(t: str) -> tuple[str, int]:
    new, n = PAT.subn(r"\\times 10^{-1\1}", t)
    return new, n

def main():
    root = pathlib.Path(__file__).parent / "whitepapers"
    total = 0
    files = 0
    for p in root.glob("*.md"):
        t = p.read_text(encoding="utf-8")
        new, n = fix_text(t)
        if n:
            p.write_text(new, encoding="utf-8")
            print(f"{p.name}: {n}")
            total += n
            files += 1
    print(f"Total replacements: {total} across {files} files")

if __name__ == "__main__":
    main()
