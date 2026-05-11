"""
Preprocess + regenerate stale PDFs.
Replaces digit-adjacent inline LaTeX math escapes ($\\cmd$) with unicode chars
that pdflatex handles via _pdf_unicode_header.tex.
Original whitepapers/*.md files are NOT modified.
"""
import os, re, subprocess, sys, time, tempfile, shutil

PANDOC = os.path.expandvars(r"%LOCALAPPDATA%\Pandoc\pandoc.exe")
PDFLATEX = r"C:\Users\tmsjd\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe"
HEADER = "_pdf_unicode_header.tex"
PDF_DIR = "pdf"

# Map LaTeX commands -> unicode equivalents (newunicodechar handles these).
MAP = {
    r"\times": "×",
    r"\pm": "±",
    r"\mp": "∓",
    r"\cdot": "·",
    r"\to": "→",
    r"\rightarrow": "→",
    r"\leftarrow": "←",
    r"\Rightarrow": "⇒",
    r"\Leftarrow": "⇐",
    r"\approx": "≈",
    r"\sim": "∼",
    r"\propto": "∝",
    r"\neq": "≠",
    r"\leq": "≤",
    r"\geq": "≥",
    r"\equiv": "≡",
    r"\gtrsim": "≳",
    r"\lesssim": "≲",
    r"\ll": "≪",
    r"\gg": "≫",
    r"\infty": "∞",
    r"\partial": "∂",
    r"\nabla": "∇",
    r"\hbar": "ℏ",
    r"\ell": "ℓ",
    r"\degree": "°",
    # Greek lowercase
    r"\alpha": "α", r"\beta": "β", r"\gamma": "γ", r"\delta": "δ",
    r"\epsilon": "ε", r"\varepsilon": "ε", r"\zeta": "ζ", r"\eta": "η",
    r"\theta": "θ", r"\vartheta": "ϑ", r"\iota": "ι", r"\kappa": "κ",
    r"\lambda": "λ", r"\mu": "μ", r"\nu": "ν", r"\xi": "ξ",
    r"\pi": "π", r"\varpi": "ϖ", r"\rho": "ρ", r"\varrho": "ϱ",
    r"\sigma": "σ", r"\varsigma": "ς", r"\tau": "τ", r"\upsilon": "υ",
    r"\phi": "φ", r"\varphi": "φ", r"\chi": "χ", r"\psi": "ψ", r"\omega": "ω",
    # Greek uppercase
    r"\Gamma": "Γ", r"\Delta": "Δ", r"\Theta": "Θ", r"\Lambda": "Λ",
    r"\Xi": "Ξ", r"\Pi": "Π", r"\Sigma": "Σ", r"\Upsilon": "Υ",
    r"\Phi": "Φ", r"\Psi": "Ψ", r"\Omega": "Ω",
}

# Sort longest first so \Gamma matches before \G etc.
SORTED_KEYS = sorted(MAP.keys(), key=len, reverse=True)
# Build a single regex: $\cmd$ where cmd is one of the known names
INLINE_PATTERN = re.compile(
    r"\$(" + "|".join(re.escape(k) for k in SORTED_KEYS) + r")\$"
)

def preprocess(text: str) -> str:
    def sub(m):
        return MAP[m.group(1)]
    return INLINE_PATTERN.sub(sub, text)

def regen(md_path: str) -> tuple[bool, str]:
    with open(md_path, "r", encoding="utf-8") as f:
        content = f.read()
    fixed = preprocess(content)
    name = os.path.splitext(os.path.basename(md_path))[0]
    pdf_path = os.path.join(PDF_DIR, name + ".pdf")
    # Write temp .md next to original so relative paths still work
    tmp_dir = os.path.dirname(md_path)
    fd, tmp_md = tempfile.mkstemp(suffix=".md", prefix="_tmp_regen_", dir=tmp_dir)
    os.close(fd)
    try:
        with open(tmp_md, "w", encoding="utf-8") as f:
            f.write(fixed)
        cmd = [
            PANDOC, tmp_md,
            f"--pdf-engine={PDFLATEX}",
            "-H", HEADER,
            "-V", "geometry:margin=1in",
            "-V", "fontsize=11pt",
            "-V", "colorlinks=true",
            "--highlight-style=tango",
            "-o", pdf_path,
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True,
                              encoding="utf-8", errors="replace", timeout=120)
        if proc.returncode == 0 and os.path.exists(pdf_path):
            return True, ""
        err = (proc.stderr or "") + (proc.stdout or "")
        line = next((l for l in err.splitlines() if l.startswith("l.") or l.startswith("! ")), "")
        return False, line[:140]
    finally:
        try:
            os.remove(tmp_md)
        except OSError:
            pass

def main():
    with open("_stale_pdfs_to_regen.txt", "r", encoding="utf-8-sig") as f:
        files = [l.strip() for l in f if l.strip()]
    total = len(files)
    ok = 0
    fail = 0
    failures = []
    t0 = time.time()
    for i, md in enumerate(files, 1):
        success, err = regen(md)
        if success:
            ok += 1
            if i % 10 == 0 or i == total:
                elapsed = time.time() - t0
                print(f"[{i}/{total}] OK  ok={ok} fail={fail} elapsed={elapsed:.0f}s", flush=True)
        else:
            fail += 1
            name = os.path.splitext(os.path.basename(md))[0]
            failures.append((name, err))
            print(f"[{i}/{total}] FAIL {name} :: {err}", flush=True)
    elapsed = time.time() - t0
    print()
    print("=" * 50)
    print(f"Regenerated OK: {ok} / {total}")
    print(f"Failed:         {fail}")
    print(f"Elapsed:        {elapsed/60:.1f} min")
    with open("_regen_stale_207_failures.txt", "w", encoding="utf-8") as f:
        for name, err in failures:
            f.write(f"{name}\t{err}\n")
    print("Failure list -> _regen_stale_207_failures.txt")

if __name__ == "__main__":
    main()
