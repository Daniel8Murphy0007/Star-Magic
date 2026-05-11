"""Strip ALL combining diacritics on Latin/Greek letters and convert to LaTeX commands.

Covers: hat (U+0302), tilde (U+0303), macron (U+0304), dot above (U+0307),
breve (U+0306), acute (U+0301), grave (U+0300), diaeresis (U+0308).
"""
import pathlib

GREEK = [
    'alpha','beta','gamma','delta','epsilon','zeta','eta','theta','iota','kappa',
    'lambda','mu','nu','xi','omicron','pi','rho','sigma','tau','upsilon','phi','chi','psi','omega',
    'Alpha','Beta','Gamma','Delta','Epsilon','Zeta','Eta','Theta','Iota','Kappa',
    'Lambda','Mu','Nu','Xi','Omicron','Pi','Rho','Sigma','Tau','Upsilon','Phi','Chi','Psi','Omega',
]

ACCENTS = {
    0x0300: 'grave', 0x0301: 'acute', 0x0302: 'hat', 0x0303: 'tilde',
    0x0304: 'bar',   0x0306: 'breve', 0x0307: 'dot',  0x0308: 'ddot',
}

LATIN = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

def fix_text(t: str) -> int:
    total = 0
    for cp, cmd in ACCENTS.items():
        ch = chr(cp)
        for letter in LATIN:
            pat = letter + ch
            if pat in t:
                n = t.count(pat)
                t = t.replace(pat, f'\\{cmd}{{{letter}}}')
                total += n
        for g in GREEK:
            pat = g + ch
            if pat in t:
                n = t.count(pat)
                t = t.replace(pat, f'\\{cmd}{{\\{g}}}')
                total += n
    return t, total

def main():
    root = pathlib.Path(__file__).parent / "whitepapers"
    grand = 0
    for p in root.glob("*.md"):
        t = p.read_text(encoding="utf-8")
        new, n = fix_text(t)
        if n:
            p.write_text(new, encoding="utf-8")
            print(f"{p.name}: {n}")
            grand += n
    print(f"Total: {grand}")

if __name__ == "__main__":
    main()
