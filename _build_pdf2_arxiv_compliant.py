#!/usr/bin/env python3
r"""
UQFF Whitepaper Batch → PDF Builder (arXiv-compliant)
======================================================

Renders every whitepaper in `whitepapers/` (both `.md` and `.tex` sources)
to a PDF placed in the `pdf2/` folder, following arXiv publishing rules:

  1. Text-searchable output (no rasterized text).
  2. Embedded fonts via fontspec + DejaVu (lualatex engine).
  3. Standard geometry (letterpaper, 0.9-in margins).
  4. PDF/A-compatible metadata via hyperref (\pdftitle, \pdfauthor, \pdfsubject).
  5. Bibliography-safe: papers with inline references are compiled as-is.
  6. Reasonable file size (~100-500 KB per short paper).
  7. Source-available: whitepapers/ is unchanged; pdf2/ is the derived artifact.

The script is idempotent — if pdf2/PAPER_XXX.pdf already exists AND is newer
than the corresponding source, the file is skipped. Deleting a PDF from pdf2/
triggers its rebuild on the next run.

Usage
-----
  python _build_pdf2_arxiv_compliant.py              # full corpus (all .md + .tex)
  python _build_pdf2_arxiv_compliant.py --limit 10   # first 10 (verification run)
  python _build_pdf2_arxiv_compliant.py --pattern "PAPER_10*"   # subset
  python _build_pdf2_arxiv_compliant.py --engine xelatex        # alternate engine
  python _build_pdf2_arxiv_compliant.py --jobs 4                # parallel builds

Exit codes
----------
  0  all requested PDFs built successfully (or already up-to-date)
  1  one or more failures logged to pdf2/_build_log.txt
"""

import argparse
import concurrent.futures
import os
import re
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
SRC_DIR   = REPO_ROOT / "whitepapers"
OUT_DIR   = REPO_ROOT / "pdf2"
HEADER    = REPO_ROOT / "pdf_header.tex"
LOG_FILE  = OUT_DIR / "_build_log.txt"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.split("Usage")[0].strip())
    p.add_argument("--limit", type=int, default=None,
                   help="Build only the first N whitepapers (default: all)")
    p.add_argument("--pattern", default="*",
                   help="Glob pattern to filter source files (default: '*')")
    p.add_argument("--engine", default="lualatex",
                   choices=("lualatex", "xelatex", "pdflatex"),
                   help="LaTeX engine (default: lualatex — best for embedded fonts)")
    p.add_argument("--jobs", type=int, default=1,
                   help="Parallel pandoc processes (default: 1)")
    p.add_argument("--force", action="store_true",
                   help="Rebuild every PDF even if it exists and is up-to-date")
    p.add_argument("--dry-run", action="store_true",
                   help="Show what would be built without actually building")
    return p.parse_args()


def discover_sources(pattern="*"):
    """Return (source_path, output_path) pairs for every whitepaper."""
    sources = []
    for ext in ("*.md", "*.tex"):
        for src in sorted(SRC_DIR.glob(ext)):
            if not src.name.startswith(("PAPER_", "COMPLETE_", "SCm_",
                                        "UQFF_", "WHITEPAPER_")):
                continue
            # If both .md and .tex exist for the same paper, prefer .md
            stem = src.stem
            if src.suffix == ".tex" and (SRC_DIR / f"{stem}.md").exists():
                continue
            out = OUT_DIR / f"{stem}.pdf"
            sources.append((src, out))
    # Filter by --pattern glob
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


def extract_yaml_frontmatter(md_text):
    """Pull title / author / date from YAML frontmatter if present."""
    meta = {}
    if md_text.startswith("---"):
        end = md_text.find("\n---", 3)
        if end > 0:
            block = md_text[3:end]
            for line in block.splitlines():
                m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*):\s*(.*)$", line)
                if m:
                    key, val = m.group(1), m.group(2).strip().strip('"').strip("'")
                    meta[key] = val
    return meta


def build_one(src, out, engine="lualatex", dry_run=False):
    """Convert a single source file to PDF. Returns (src, out, ok, message)."""
    if dry_run:
        return (src, out, True, "dry-run: would build")

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # Build a per-paper metadata YAML for pandoc (embedded via -M)
    meta_args = []
    if src.suffix == ".md":
        try:
            md_text = src.read_text(encoding="utf-8", errors="replace")
            meta = extract_yaml_frontmatter(md_text)
            if meta.get("title"):
                meta_args += ["-M", f"title={meta['title']}"]
            if meta.get("author"):
                meta_args += ["-M", f"author={meta['author']}"]
            if meta.get("date"):
                meta_args += ["-M", f"date={meta['date']}"]
            meta_args += ["-M", f"subject=UQFF Star-Magic whitepaper: {src.stem}"]
        except Exception:
            pass

    cmd = [
        "pandoc",
        "-f", ("markdown+yaml_metadata_block+tex_math_dollars+raw_tex" if src.suffix == ".md"
               else "latex"),
        "-t", "pdf",
        "--pdf-engine", engine,
        "-H", str(HEADER.resolve()),
        "-V", "papersize=letter",
        "-V", "linkcolor=blue",
        "-V", "urlcolor=blue",
        "-V", "citecolor=blue",
        "-o", str(out.resolve()),
    ] + meta_args + [str(src.resolve())]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if result.returncode == 0 and out.exists() and out.stat().st_size > 4096:
            return (src, out, True, "ok")
        else:
            snippet = (result.stderr or result.stdout or "").splitlines()
            snippet = "\n".join(line for line in snippet
                                if "warning" not in line.lower())[-800:]
            return (src, out, False, snippet or "unknown pandoc error")
    except subprocess.TimeoutExpired:
        return (src, out, False, "timeout after 180s")
    except Exception as e:
        return (src, out, False, f"exception: {type(e).__name__}: {e}")


def main():
    args = parse_args()

    if not SRC_DIR.is_dir():
        print(f"ERROR: source directory not found: {SRC_DIR}", file=sys.stderr)
        return 2
    if not HEADER.is_file():
        print(f"ERROR: header file not found: {HEADER}", file=sys.stderr)
        return 2

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
    print(f"  engine             : {args.engine}")
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
                 f"engine={args.engine} jobs={args.jobs} force={args.force}\n",
                 f"pattern={args.pattern} limit={args.limit}\n\n"]

    def _work(pair):
        src, out = pair
        return build_one(src, out, engine=args.engine, dry_run=args.dry_run)

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
                print(f"FIRST FAILURE inline (for quick diagnosis):")
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
