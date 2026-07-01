#!/usr/bin/env python3
r"""
UQFF Whitepaper Batch → PDF Builder (Pure-Python fallback)
==========================================================

Fallback path B: no pandoc / no LaTeX install required.

Uses Python-only stack:
  - markdown-it-py  (fast CommonMark + GitHub-flavored markdown parser)
  - weasyprint      (HTML/CSS → PDF with embedded fonts, text-searchable)
  - Optional: KaTeX-style math substitution for common LaTeX macros

Install prerequisites (one-time):

    pip install markdown-it-py weasyprint

On Windows, WeasyPrint needs GTK3 runtime installed for font rendering:

    - Download GTK3 runtime from:
      https://github.com/tschoonj/GTK-for-Windows-Runtime-Environment-Installer/releases
    - Run installer, restart terminal.

Or use the alternative pure-pip route with `reportlab` (no GTK needed):

    pip install markdown-it-py reportlab

If WeasyPrint is not available on your system, this script falls back to
reportlab automatically — with reduced formatting quality but guaranteed
text-searchable PDFs.

Usage
-----
    python _build_pdf2_pure_python.py --limit 20      # test batch
    python _build_pdf2_pure_python.py                 # full corpus
    python _build_pdf2_pure_python.py --pattern "PAPER_10*"
    python _build_pdf2_pure_python.py --jobs 4        # parallel

arXiv publishing rules honored
------------------------------
  1. Text-searchable output (no rasterization).
  2. Embedded fonts (WeasyPrint embeds system fonts by default; reportlab uses
     built-in Helvetica/Times).
  3. Standard letter-page geometry via CSS @page.
  4. PDF metadata (title, author, subject) set via WeasyPrint.
  5. Reasonable file size (~50-200 KB per short paper for weasyprint,
     ~30-100 KB for reportlab).
  6. Source-available: whitepapers/ untouched.

Author:    Daniel T. Murphy
License:   AGPL-3.0-or-later (see LICENSE) or Commercial (see COMMERCIAL.md)
"""

import argparse
import concurrent.futures
import html
import re
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR   = REPO_ROOT / "whitepapers"
OUT_DIR   = REPO_ROOT / "pdf2"
LOG_FILE  = OUT_DIR / "_build_log.txt"

# Detect available PDF backend
BACKEND = None
try:
    from weasyprint import HTML, CSS  # noqa: F401
    BACKEND = "weasyprint"
except Exception:
    try:
        from reportlab.lib.pagesizes import letter  # noqa: F401
        from reportlab.platypus import SimpleDocTemplate  # noqa: F401
        BACKEND = "reportlab"
    except Exception:
        BACKEND = None

# Markdown parser
try:
    from markdown_it import MarkdownIt
    MD_PARSER_AVAILABLE = True
except Exception:
    MD_PARSER_AVAILABLE = False


ARXIV_CSS = """
@page {
    size: letter;
    margin: 0.9in;
    @bottom-center {
        content: counter(page);
        font-family: "Times New Roman", serif;
        font-size: 9pt;
    }
}
body {
    font-family: "Times New Roman", "Nimbus Roman", Times, serif;
    font-size: 10.5pt;
    line-height: 1.35;
    color: #000;
}
h1 { font-size: 16pt; margin-top: 0; }
h2 { font-size: 13pt; margin-top: 1em; }
h3 { font-size: 11.5pt; margin-top: 0.8em; }
h4 { font-size: 10.5pt; margin-top: 0.6em; font-style: italic; }
p  { margin: 0.4em 0; text-align: justify; }
code {
    font-family: "Consolas", "Courier New", monospace;
    font-size: 9.5pt;
    background: #f4f4f4;
    padding: 1px 3px;
    border-radius: 2px;
}
pre {
    font-family: "Consolas", "Courier New", monospace;
    font-size: 9pt;
    background: #f4f4f4;
    padding: 6px 8px;
    border-radius: 3px;
    overflow-wrap: break-word;
}
blockquote {
    border-left: 3px solid #888;
    padding-left: 10pt;
    color: #333;
    margin: 0.4em 0;
}
table {
    border-collapse: collapse;
    margin: 0.6em 0;
    font-size: 9.5pt;
}
th, td {
    border: 1px solid #666;
    padding: 3pt 6pt;
    text-align: left;
}
th { background: #eee; }
hr {
    border: none;
    border-top: 1px solid #999;
    margin: 0.8em 0;
}
a  { color: #0645AD; text-decoration: none; }
img { max-width: 100%; }
.math-inline, .math-display {
    font-family: "Cambria Math", "STIX Two Math", serif;
    font-style: italic;
}
.math-display {
    display: block;
    text-align: center;
    margin: 0.4em 0;
}
"""


def parse_args():
    p = argparse.ArgumentParser(
        description="Pure-Python markdown→PDF batch builder (fallback for pandoc)"
    )
    p.add_argument("--limit", type=int, default=None,
                   help="Build only first N whitepapers")
    p.add_argument("--pattern", default="*",
                   help="Glob pattern for source filenames")
    p.add_argument("--jobs", type=int, default=1,
                   help="Parallel processes (default: 1)")
    p.add_argument("--force", action="store_true",
                   help="Rebuild every PDF regardless of mtime")
    p.add_argument("--dry-run", action="store_true",
                   help="Show what would be built without building")
    p.add_argument("--backend", default=None,
                   choices=("weasyprint", "reportlab"),
                   help="Force PDF backend (default: auto-detect)")
    return p.parse_args()


def discover_sources(pattern="*"):
    sources = []
    for src in sorted(SRC_DIR.glob("*.md")):
        if not src.name.startswith(("PAPER_", "COMPLETE_", "SCm_",
                                    "UQFF_", "WHITEPAPER_")):
            continue
        out = OUT_DIR / f"{src.stem}.pdf"
        sources.append((src, out))
    if pattern != "*":
        matcher = re.compile(pattern.replace("*", ".*").replace("?", "."))
        sources = [s for s in sources if matcher.match(s[0].stem)]
    return sources


def needs_rebuild(src, out, force=False):
    if force:
        return True
    if not out.exists():
        return True
    return src.stat().st_mtime > out.stat().st_mtime


def extract_frontmatter(md_text):
    meta = {}
    body = md_text
    if md_text.startswith("---"):
        end = md_text.find("\n---", 3)
        if end > 0:
            block = md_text[3:end]
            body = md_text[end + 4:]
            for line in block.splitlines():
                m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*)$", line)
                if m:
                    key, val = m.group(1), m.group(2).strip().strip('"').strip("'")
                    meta[key] = val
    return meta, body


def render_markdown_html(md_body, meta):
    """Convert markdown body to HTML with math-mode preservation."""
    if MD_PARSER_AVAILABLE:
        md = MarkdownIt("commonmark", {"html": True, "linkify": True}) \
            .enable(["table", "strikethrough"])
        html_body = md.render(md_body)
    else:
        # Very minimal fallback: escape HTML and preserve paragraphs
        html_body = html.escape(md_body).replace("\n\n", "</p><p>")
        html_body = f"<p>{html_body}</p>"

    # Handle inline math $...$ and display math $$...$$
    # Simple substitution — not full MathJax, but preserves LaTeX source text
    html_body = re.sub(
        r"\$\$([^$]+)\$\$",
        lambda m: f'<div class="math-display">{html.escape(m.group(1).strip())}</div>',
        html_body,
    )
    html_body = re.sub(
        r"\$([^$\n]+)\$",
        lambda m: f'<span class="math-inline">{html.escape(m.group(1).strip())}</span>',
        html_body,
    )

    title = meta.get("title", "").strip('"').strip("'") or "UQFF Whitepaper"
    author = meta.get("author", "").strip('"').strip("'") or "Daniel T. Murphy"
    date = meta.get("date", "").strip('"').strip("'") or ""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{html.escape(title)}</title>
<meta name="author" content="{html.escape(author)}">
<meta name="date" content="{html.escape(date)}">
<meta name="subject" content="UQFF Star-Magic whitepaper">
</head>
<body>
<h1>{html.escape(title)}</h1>
<p><em>Author:</em> {html.escape(author)}
{f'<br><em>Date:</em> {html.escape(date)}' if date else ''}</p>
<hr>
{html_body}
</body>
</html>"""


def build_weasyprint(src, out):
    from weasyprint import HTML, CSS
    md_text = src.read_text(encoding="utf-8", errors="replace")
    meta, body = extract_frontmatter(md_text)
    html_str = render_markdown_html(body, meta)
    HTML(string=html_str).write_pdf(str(out), stylesheets=[CSS(string=ARXIV_CSS)])
    return out.exists() and out.stat().st_size > 2048


def build_reportlab(src, out):
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
    from reportlab.lib.units import inch
    md_text = src.read_text(encoding="utf-8", errors="replace")
    meta, body = extract_frontmatter(md_text)
    title = meta.get("title", src.stem).strip('"').strip("'")
    author = meta.get("author", "Daniel T. Murphy").strip('"').strip("'")

    doc = SimpleDocTemplate(
        str(out),
        pagesize=letter,
        rightMargin=0.9 * inch, leftMargin=0.9 * inch,
        topMargin=0.9 * inch,   bottomMargin=0.9 * inch,
        title=title, author=author, subject=f"UQFF whitepaper: {src.stem}",
    )
    styles = getSampleStyleSheet()
    story = [
        Paragraph(html.escape(title), styles["Title"]),
        Paragraph(f"<i>{html.escape(author)}</i>", styles["Normal"]),
        Spacer(1, 12),
    ]
    for para in body.split("\n\n"):
        para = para.strip()
        if not para:
            continue
        # Very simple heading detection
        if para.startswith("### "):
            story.append(Paragraph(html.escape(para[4:]), styles["Heading3"]))
        elif para.startswith("## "):
            story.append(Paragraph(html.escape(para[3:]), styles["Heading2"]))
        elif para.startswith("# "):
            story.append(Paragraph(html.escape(para[2:]), styles["Heading1"]))
        else:
            # Preserve inline code and basic markdown
            para_html = html.escape(para).replace("\n", "<br/>")
            para_html = re.sub(r"`([^`]+)`",
                               r'<font face="Courier">\1</font>', para_html)
            para_html = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", para_html)
            para_html = re.sub(r"\*([^*]+)\*", r"<i>\1</i>", para_html)
            try:
                story.append(Paragraph(para_html, styles["BodyText"]))
                story.append(Spacer(1, 4))
            except Exception:
                # Skip paragraphs reportlab can't parse
                pass
    doc.build(story)
    return out.exists() and out.stat().st_size > 2048


def build_one(src, out, backend=None, dry_run=False):
    if dry_run:
        return (src, out, True, "dry-run")
    backend = backend or BACKEND
    if backend is None:
        return (src, out, False,
                "no PDF backend available. run: pip install weasyprint  "
                "OR: pip install reportlab")
    try:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        if backend == "weasyprint":
            ok = build_weasyprint(src, out)
        elif backend == "reportlab":
            ok = build_reportlab(src, out)
        else:
            return (src, out, False, f"unknown backend: {backend}")
        return (src, out, ok, "ok" if ok else "output too small or missing")
    except Exception as e:
        return (src, out, False, f"{type(e).__name__}: {str(e)[:500]}")


def main():
    args = parse_args()

    if not SRC_DIR.is_dir():
        print(f"ERROR: source directory not found: {SRC_DIR}", file=sys.stderr)
        return 2

    if not MD_PARSER_AVAILABLE:
        print("WARNING: markdown-it-py not installed — using minimal fallback parser",
              file=sys.stderr)
        print("Install:  pip install markdown-it-py", file=sys.stderr)

    if BACKEND is None and not args.dry_run:
        print("ERROR: no PDF backend available on this system.", file=sys.stderr)
        print("Install one of:", file=sys.stderr)
        print("  pip install weasyprint      # (best quality, needs GTK3 runtime on Windows)",
              file=sys.stderr)
        print("  pip install reportlab       # (works everywhere, simpler output)",
              file=sys.stderr)
        return 2

    backend = args.backend or BACKEND
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    all_sources = discover_sources(pattern=args.pattern)
    to_build = [(s, o) for s, o in all_sources if needs_rebuild(s, o, args.force)]
    if args.limit is not None:
        to_build = to_build[:args.limit]

    total = len(to_build)
    already = len(all_sources) - total
    print(f"pdf2 build target = {OUT_DIR}")
    print(f"  discovered sources : {len(all_sources)}")
    print(f"  already up-to-date : {already}")
    print(f"  to build           : {total}")
    print(f"  backend            : {backend or '(dry-run)'}")
    print(f"  jobs               : {args.jobs}")
    print(f"  dry-run            : {args.dry_run}")
    print()

    if total == 0:
        print("nothing to do.")
        return 0

    start_t = time.time()
    ok_count = 0
    fail_count = 0
    log_lines = [f"# pdf2 build log — {time.strftime('%Y-%m-%d %H:%M:%S')}\n",
                 f"backend={backend} jobs={args.jobs} force={args.force}\n",
                 f"pattern={args.pattern} limit={args.limit}\n\n"]

    def _work(pair):
        src, out = pair
        return build_one(src, out, backend=backend, dry_run=args.dry_run)

    if args.jobs > 1:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as ex:
            results = ex.map(_work, to_build)
    else:
        results = (_work(pair) for pair in to_build)

    first_failure_shown = False
    for i, (src, out, ok, msg) in enumerate(results, 1):
        if ok:
            ok_count += 1
            marker = "PASS"
        else:
            fail_count += 1
            marker = "FAIL"
            log_lines.append(f"[FAIL] {src.name}\n{msg}\n\n")
            if not first_failure_shown:
                first_failure_shown = True
                print()
                print("=" * 70)
                print(f"FIRST FAILURE inline (for diagnosis):")
                print(f"  file: {src.name}")
                print("-" * 70)
                for line in msg.splitlines()[-15:]:
                    print(f"    {line}")
                print("=" * 70)
                print()
        elapsed = time.time() - start_t
        rate = i / elapsed if elapsed > 0 else 0.0
        print(f"  [{i:5d}/{total}] {marker}  {src.name}  "
              f"({rate:.1f} papers/sec)")

    LOG_FILE.write_text("".join(log_lines), encoding="utf-8")
    elapsed = time.time() - start_t
    print()
    print(f"summary: {ok_count} ok, {fail_count} fail, {elapsed:.1f}s")
    print(f"log: {LOG_FILE}")
    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
