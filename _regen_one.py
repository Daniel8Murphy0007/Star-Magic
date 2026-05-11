"""
Single-paper rebuild helper.
Usage: python _regen_one.py PAPER_016_Quantum_Entanglement_UQFF_Nonlocal_Correlations
Reads whitepapers/<name>.md, preprocesses via same map as _regen_stale_207.py,
builds pdf/<name>.pdf. Prints full pdflatex error tail on failure.
"""
import os, re, subprocess, sys, tempfile

PANDOC = os.path.expandvars(r"%LOCALAPPDATA%\Pandoc\pandoc.exe")
PDFLATEX = r"C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
HEADER = "_pdf_unicode_header.tex"
PDF_DIR = "pdf"
WP_DIR = "whitepapers"

MAP = {
    r"\times": "×", r"\pm": "±", r"\mp": "∓", r"\cdot": "·",
    r"\to": "→", r"\rightarrow": "→", r"\leftarrow": "←",
    r"\Rightarrow": "⇒", r"\Leftarrow": "⇐",
    r"\approx": "≈", r"\sim": "∼", r"\propto": "∝",
    r"\neq": "≠", r"\leq": "≤", r"\geq": "≥", r"\equiv": "≡",
    r"\gtrsim": "≳", r"\lesssim": "≲", r"\ll": "≪", r"\gg": "≫",
    r"\infty": "∞", r"\partial": "∂", r"\nabla": "∇",
    r"\hbar": "ℏ", r"\ell": "ℓ", r"\degree": "°",
    r"\alpha":"α", r"\beta":"β", r"\gamma":"γ", r"\delta":"δ",
    r"\epsilon":"ε", r"\varepsilon":"ε", r"\zeta":"ζ", r"\eta":"η",
    r"\theta":"θ", r"\vartheta":"ϑ", r"\iota":"ι", r"\kappa":"κ",
    r"\lambda":"λ", r"\mu":"μ", r"\nu":"ν", r"\xi":"ξ",
    r"\pi":"π", r"\varpi":"ϖ", r"\rho":"ρ", r"\varrho":"ϱ",
    r"\sigma":"σ", r"\varsigma":"ς", r"\tau":"τ", r"\upsilon":"υ",
    r"\phi":"φ", r"\varphi":"φ", r"\chi":"χ", r"\psi":"ψ", r"\omega":"ω",
    r"\Gamma":"Γ", r"\Delta":"Δ", r"\Theta":"Θ", r"\Lambda":"Λ",
    r"\Xi":"Ξ", r"\Pi":"Π", r"\Sigma":"Σ", r"\Upsilon":"Υ",
    r"\Phi":"Φ", r"\Psi":"Ψ", r"\Omega":"Ω",
}
SORTED_KEYS = sorted(MAP.keys(), key=len, reverse=True)
INLINE = re.compile(r"\$(" + "|".join(re.escape(k) for k in SORTED_KEYS) + r")\$")

# Bare-command pattern (outside math + code) for replacement with unicode
BARE = re.compile(r"\\(" + "|".join(re.escape(k[1:]) for k in SORTED_KEYS) + r")(?![A-Za-z])")

def preprocess(text):
    text = INLINE.sub(lambda m: MAP[m.group(1)], text)
    # Now walk lines: skip inside fenced code blocks and $$..$$ blocks
    lines = text.split("\n")
    out = []
    in_code = False
    in_block = False
    for line in lines:
        stripped = line.lstrip()
        if stripped.startswith("```"):
            in_code = not in_code
            out.append(line); continue
        if in_code:
            out.append(line); continue
        dd = line.count("$$")
        if in_block:
            out.append(line)
            if dd % 2 == 1:
                in_block = False
            continue
        if dd >= 2:
            out.append(line); continue
        if dd == 1:
            out.append(line); in_block = True; continue
        # not in math block; replace bare commands outside $...$ regions
        # walk char by char tracking single $
        s = line
        out_chars = []
        j = 0
        inm = False
        while j < len(s):
            c = s[j]
            if c == "$":
                inm = not inm
                out_chars.append(c); j += 1; continue
            if not inm and c == "\\" and j+1 < len(s):
                m = BARE.match(s, j)
                if m:
                    cmd = "\\" + m.group(1)
                    out_chars.append(MAP[cmd])
                    j = m.end(); continue
            out_chars.append(c); j += 1
        out.append("".join(out_chars))
    return "\n".join(out)

def main():
    if len(sys.argv) != 2:
        print("Usage: python _regen_one.py <PAPER_NAME_no_extension>")
        sys.exit(2)
    name = sys.argv[1].replace(".md", "")
    md = os.path.join(WP_DIR, name + ".md")
    if not os.path.exists(md):
        print(f"NOT FOUND: {md}")
        sys.exit(2)
    with open(md, "r", encoding="utf-8") as f:
        content = f.read()
    fixed = preprocess(content)
    fd, tmp = tempfile.mkstemp(suffix=".md", prefix="_tmp_regen_", dir=WP_DIR)
    os.close(fd)
    pdf = os.path.join(PDF_DIR, name + ".pdf")
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(fixed)
        proc = subprocess.run(
            [PANDOC, tmp, f"--pdf-engine={PDFLATEX}", "-H", HEADER,
             "-V", "geometry:margin=1in", "-V", "fontsize=11pt",
             "-V", "colorlinks=true", "--highlight-style=tango", "-o", pdf],
            capture_output=True, text=True, encoding="utf-8",
            errors="replace", timeout=120,
        )
    finally:
        try: os.remove(tmp)
        except OSError: pass
    if proc.returncode == 0 and os.path.exists(pdf):
        size = os.path.getsize(pdf)
        print(f"OK   {name}.pdf  ({size//1024} KB)")
        sys.exit(0)
    else:
        err = (proc.stderr or "") + (proc.stdout or "")
        # Print tail of error
        tail = err.splitlines()
        for line in tail[-30:]:
            try:
                print(line)
            except UnicodeEncodeError:
                print(line.encode("ascii", "replace").decode("ascii"))
        print(f"\nFAIL  exit={proc.returncode}")
        sys.exit(1)

if __name__ == "__main__":
    main()
