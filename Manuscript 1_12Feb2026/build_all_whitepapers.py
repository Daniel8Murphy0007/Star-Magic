# -*- coding: utf-8 -*-
"""
build_all_whitepapers.py
Universal markdown → PDF converter for the Star-Magic UQFF whitepaper series.
Processes every PAPER_*.md in:
  - <repo>/whitepapers/
  - <repo>/  (root, excluding already-covered ones)
Outputs all PDFs to:  <repo>/Manuscript 1_12Feb2026/whitepapers_pdf/
Requires: pip install reportlab
"""

import os, re, sys

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import inch
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        HRFlowable, KeepTogether
    )
except ImportError:
    sys.exit("reportlab not found – run: pip install reportlab")

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_DIR   = os.path.dirname(SCRIPT_DIR)
WP_DIR     = os.path.join(REPO_DIR, "whitepapers")
OUT_DIR    = os.path.join(SCRIPT_DIR, "whitepapers_pdf")
os.makedirs(OUT_DIR, exist_ok=True)

# ── Styles ───────────────────────────────────────────────────────────────────
W, H = A4
TW   = W - 2.2 * inch

DARK  = colors.HexColor("#1a1a2e")
BLUE  = colors.HexColor("#16213e")
RULE  = colors.HexColor("#0f3460")
LGREY = colors.HexColor("#f5f5f5")
TBLHD = colors.HexColor("#d0e4f7")
HBOX  = colors.HexColor("#eef4ff")

def _sty(name, font="Times-Roman", size=10, lead=14, align=TA_LEFT,
         color=colors.black, sb=4, sa=3, li=0, fl=0, bg=None, bw=0, bp=0):
    kw = dict(fontName=font, fontSize=size, leading=lead, alignment=align,
              textColor=color, spaceBefore=sb, spaceAfter=sa,
              leftIndent=li, firstLineIndent=fl)
    if bg: kw["backColor"] = bg
    if bw: kw.update(borderColor=RULE, borderWidth=bw, borderPad=bp)
    return ParagraphStyle(name, **kw)

META  = _sty("ME", "Times-Italic",  9, 12, TA_CENTER, colors.grey,   sb=0, sa=2)
TITLE = _sty("TI", "Times-Bold",   16, 22, TA_CENTER, DARK,          sb=0, sa=4)
AUTH  = _sty("AU", "Times-Roman",  11, 15, TA_CENTER, DARK,          sb=0, sa=2)
ABS   = _sty("AB", "Times-Roman",  10, 14, TA_JUSTIFY,colors.black,  li=20, sa=4)
H1    = _sty("H1", "Times-Bold",   13, 17, TA_LEFT,   DARK,          sb=12, sa=3)
H2    = _sty("H2", "Times-Bold",   11, 15, TA_LEFT,   BLUE,          sb=8,  sa=2)
H3    = _sty("H3", "Times-BoldItalic",10,14,TA_LEFT,  DARK,          sb=6,  sa=2)
H4    = _sty("H4", "Times-Bold",   10, 13, TA_LEFT,   colors.black,  sb=4,  sa=1)
BODY  = _sty("BO", "Times-Roman",  10, 14, TA_JUSTIFY,colors.black,  sa=3)
CODE  = _sty("CO", "Courier",       8, 11, TA_LEFT,   colors.black,  li=10, sa=2, sb=2, bg=LGREY)
BULL  = _sty("BU", "Times-Roman",  10, 14, TA_LEFT,   colors.black,  li=18, fl=-9, sa=2)
NUM   = _sty("NU", "Times-Roman",  10, 14, TA_LEFT,   colors.black,  li=24, fl=-14, sa=2)
CAP   = _sty("CA", "Times-Italic",  8, 11, TA_CENTER, colors.grey,   sa=6)
FOOT  = _sty("FO", "Times-Italic",  8, 11, TA_CENTER, colors.grey,   sb=0, sa=0)

HSTYLES = [None, H1, H2, H3, H4]   # index by heading level 1-4

# ── Markdown helpers ─────────────────────────────────────────────────────────
_ESC = [("&", "&amp;"), ("<", "&lt;"), (">", "&gt;")]

def _escape(t):
    for a, b in _ESC:
        t = t.replace(a, b)
    return t

def _inline(text):
    """Convert inline markdown to ReportLab XML markup.

    Uses a placeholder strategy: code/math spans are extracted first so that
    bold/italic regexes cannot match across or inside them.
    """
    saved = {}
    counter = [0]

    def _save(markup):
        key = "\x00PH%d\x00" % counter[0]
        saved[key] = markup
        counter[0] += 1
        return key

    # ── Step 1: extract code/math spans before any other processing ─────────
    # $$...$$
    t = re.sub(r'\$\$(.+?)\$\$',
               lambda m: _save('<font face="Courier" size="8">%s</font>' % _escape(m.group(1))),
               text)
    # $...$
    t = re.sub(r'\$(.+?)\$',
               lambda m: _save('<font face="Courier" size="9">%s</font>' % _escape(m.group(1))),
               t)
    # `code`
    t = re.sub(r'`([^`]+)`',
               lambda m: _save('<font face="Courier" size="8">%s</font>' % _escape(m.group(1))),
               t)

    # ── Step 2: escape HTML in the remaining text (placeholders are safe) ───
    t = _escape(t)

    # ── Step 3: bold / italic ────────────────────────────────────────────────
    # ***bold-italic***
    t = re.sub(r'\*\*\*(.+?)\*\*\*', r'<b><i>\1</i></b>', t)
    # **bold**
    t = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', t)
    # __bold__
    t = re.sub(r'__(.+?)__', r'<b>\1</b>', t)
    # *italic* or _italic_  (won't touch text inside saved spans)
    t = re.sub(r'\*([^*\n]+?)\*', r'<i>\1</i>', t)
    t = re.sub(r'_([^_\n]+?)_',   r'<i>\1</i>', t)
    # Unicode -> HTML entities for common symbols
    t = t.replace('×', '&#215;').replace('→', '&#8594;').replace('←', '&#8592;')
    t = t.replace('−', '&#8722;').replace('≈', '&#8776;').replace('≤', '&#8804;')
    t = t.replace('≥', '&#8805;').replace('≠', '&#8800;').replace('∞', '&#8734;')
    t = t.replace('∑', '&#8721;').replace('∫', '&#8747;').replace('√', '&#8730;')
    t = t.replace('π', '&#960;').replace('Λ', '&#923;').replace('λ', '&#955;')
    t = t.replace('α', '&#945;').replace('β', '&#946;').replace('γ', '&#947;')
    t = t.replace('δ', '&#948;').replace('ε', '&#949;').replace('θ', '&#952;')
    t = t.replace('κ', '&#954;').replace('ω', '&#969;').replace('Ω', '&#937;')
    t = t.replace('ρ', '&#961;').replace('σ', '&#963;').replace('τ', '&#964;')
    t = t.replace('φ', '&#966;').replace('ψ', '&#968;').replace('η', '&#951;')
    t = t.replace('ℏ', '&#8463;').replace('ħ', '&#8463;')
    t = t.replace('☉', '&#9737;').replace('⊙', '&#8857;')
    t = t.replace('™', '&#8482;').replace('©', '&#169;').replace('®', '&#174;')
    t = t.replace('…', '...').replace('\u2014', '&#8212;').replace('\u2013', '&#8211;')
    t = t.replace('\u2019', "'").replace('\u2018', "'")
    t = t.replace('\u201c', '"').replace('\u201d', '"')
    # superscript: ² ³ etc.
    for raw, ent in [('²','<super>2</super>'),('³','<super>3</super>'),
                     ('⁻','<super>-</super>'),('⁰','<super>0</super>'),
                     ('¹','<super>1</super>')]:
        t = t.replace(raw, ent)
    # subscript notation like _{x}
    t = re.sub(r'_\{([^}]+)\}', lambda m: '<sub>%s</sub>' % m.group(1), t)
    t = re.sub(r'\^\{([^}]+)\}', lambda m: '<super>%s</super>' % m.group(1), t)

    # ── Step 5: restore saved code/math spans ───────────────────────────────
    for key, val in saved.items():
        t = t.replace(key, val)

    return t

def _p(text, style=None):
    try:
        return Paragraph(_inline(text), style or BODY)
    except Exception:
        return Paragraph(_escape(text), style or BODY)

def _rule():
    return HRFlowable(width="100%", thickness=0.5, color=RULE, spaceAfter=3, spaceBefore=3)

def _sp(n=6):
    return Spacer(1, n)

# ── Table parser ─────────────────────────────────────────────────────────────
def _parse_table(lines):
    """Parse markdown table lines → ReportLab Table."""
    rows = []
    for line in lines:
        if re.match(r'^\|[-: |]+\|?\s*$', line):
            continue  # separator row
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        row = []
        for c in cells:
            try:
                row.append(Paragraph(_inline(c), _sty("TC","Times-Roman",8,11)))
            except Exception:
                row.append(Paragraph(_escape(c), _sty("TC","Times-Roman",8,11)))
        if row:
            rows.append(row)
    if not rows:
        return []
    # Uniform column widths
    ncols = max(len(r) for r in rows)
    # Pad short rows
    for r in rows:
        while len(r) < ncols:
            r.append(Paragraph("", _sty("TC","Times-Roman",8,11)))
    cw = TW / ncols
    t = Table(rows, colWidths=[cw]*ncols, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",  (0,0), (-1,0),  TBLHD),
        ("FONTNAME",    (0,0), (-1,0),  "Times-Bold"),
        ("FONTSIZE",    (0,0), (-1,-1), 8),
        ("GRID",        (0,0), (-1,-1), 0.35, colors.grey),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white, LGREY]),
        ("ALIGN",       (0,0), (-1,-1), "CENTER"),
        ("VALIGN",      (0,0), (-1,-1), "MIDDLE"),
        ("TOPPADDING",  (0,0), (-1,-1), 3),
        ("BOTTOMPADDING",(0,0),(-1,-1), 3),
        ("LEFTPADDING", (0,0), (-1,-1), 4),
        ("RIGHTPADDING",(0,0),(-1,-1),  4),
    ]))
    return [_sp(3), t, _sp(5)]

# ── Main markdown parser ─────────────────────────────────────────────────────
def md_to_story(md_text, paper_id=""):
    """Convert full markdown text to a reportlab story list."""
    story = []
    lines = md_text.splitlines()
    i = 0
    in_code   = False
    code_buf  = []
    tbl_buf   = []
    in_table  = False

    def flush_code():
        nonlocal code_buf
        if code_buf:
            for cl in code_buf:
                story.append(Paragraph(_escape(cl) if cl.strip() else " ", CODE))
            story.append(_sp(4))
            code_buf = []

    def flush_table():
        nonlocal tbl_buf
        if tbl_buf:
            story.extend(_parse_table(tbl_buf))
            tbl_buf = []

    while i < len(lines):
        line = lines[i]
        raw  = line.rstrip()

        # ── Code fences ────────────────────────────────────────────────────
        if raw.startswith("```") or raw.startswith("~~~"):
            if in_code:
                flush_code()
                in_code = False
            else:
                flush_table(); in_table = False
                in_code = True
            i += 1; continue

        if in_code:
            code_buf.append(raw)
            i += 1; continue

        # ── Table rows ─────────────────────────────────────────────────────
        if raw.startswith("|"):
            in_table = True
            tbl_buf.append(raw)
            i += 1; continue
        else:
            if in_table:
                flush_table(); in_table = False

        # ── Blank line ─────────────────────────────────────────────────────
        if not raw.strip():
            story.append(_sp(4))
            i += 1; continue

        # ── Horizontal rule ────────────────────────────────────────────────
        if re.match(r'^[-*_]{3,}\s*$', raw):
            story.append(_rule())
            i += 1; continue

        # ── ATX Headings ───────────────────────────────────────────────────
        hm = re.match(r'^(#{1,4})\s+(.*)', raw)
        if hm:
            level = min(len(hm.group(1)), 4)
            text  = hm.group(2).strip()
            # strip trailing # marks
            text  = re.sub(r'\s+#+\s*$', '', text)
            style = HSTYLES[level]
            story.append(_p(text, style))
            i += 1; continue

        # ── Setext headings ────────────────────────────────────────────────
        if i+1 < len(lines):
            nxt = lines[i+1].rstrip()
            if re.match(r'^=+\s*$', nxt):
                story.append(_p(raw, H1))
                i += 2; continue
            if re.match(r'^-+\s*$', nxt) and len(nxt) >= 2:
                story.append(_p(raw, H2))
                i += 2; continue

        # ── Block-quote (treat as note box) ────────────────────────────────
        if raw.startswith(">"):
            text = re.sub(r'^>\s?', '', raw)
            bq = _sty("BQ","Times-Italic",9,13,TA_LEFT,DARK,li=14,sb=2,sa=2,bg=HBOX,bw=1,bp=5)
            story.append(Paragraph(_inline(text), bq))
            i += 1; continue

        # ── Unordered list ────────────────────────────────────────────────
        bm = re.match(r'^(\s*)[-*+]\s+(.*)', raw)
        if bm:
            text = bm.group(2)
            story.append(_p(text, BULL))
            i += 1; continue

        # ── Ordered list ─────────────────────────────────────────────────
        nm = re.match(r'^(\s*)\d+[.)]\s+(.*)', raw)
        if nm:
            text = nm.group(2)
            story.append(_p(text, NUM))
            i += 1; continue

        # ── LaTeX display math ($$...$$) multi-line ────────────────────────
        if raw.strip().startswith("$$"):
            math_lines = []
            rest = raw.strip()[2:]
            if rest.endswith("$$"):
                math_lines.append(rest[:-2])
                i += 1
            else:
                if rest.strip():
                    math_lines.append(rest)
                i += 1
                while i < len(lines):
                    ml = lines[i].rstrip()
                    if ml.strip().endswith("$$"):
                        math_lines.append(ml.strip().rstrip("$").rstrip())
                        i += 1; break
                    math_lines.append(ml)
                    i += 1
            for ml in math_lines:
                clean = ml.strip()
                if clean:
                    story.append(Paragraph(_escape(clean), CODE))
            story.append(_sp(3))
            continue

        # ── Regular paragraph ─────────────────────────────────────────────
        # Collect continuation lines
        para_lines = [raw]
        i += 1
        while i < len(lines):
            nline = lines[i].rstrip()
            if (not nline.strip() or nline.startswith("#") or
                    nline.startswith("|") or nline.startswith("```") or
                    nline.startswith("~~~") or nline.startswith(">") or
                    re.match(r'^(\s*)[-*+]\s+', nline) or
                    re.match(r'^(\s*)\d+[.)]\s+', nline) or
                    re.match(r'^[-*_]{3,}\s*$', nline) or
                    nline.strip().startswith("$$")):
                break
            para_lines.append(nline)
            i += 1

        full = " ".join(pl.strip() for pl in para_lines if pl.strip())
        if full:
            story.append(_p(full, BODY))

    flush_code()
    flush_table()
    return story


# ── Metadata extraction ──────────────────────────────────────────────────────
def extract_meta(md_text, filename):
    """Extract PAPER_XXX id, title, author, session from markdown."""
    paper_id = re.sub(r'\.md$', '', filename.split("_")[0] + "_" +
                      filename.split("_")[1] if "_" in filename else filename, flags=re.I)
    # Try to get paper number from filename
    pid_m  = re.match(r'(PAPER_\d+)', filename, re.I)
    paper_id = pid_m.group(1) if pid_m else filename.replace(".md","")

    lines = md_text.splitlines()
    title_line = ""
    author_line = ""
    session_line = ""
    for ln in lines[:20]:
        s = ln.strip().lstrip("#").strip()
        if not title_line and s and not s.startswith("Author") and not s.startswith("Session"):
            title_line = s
        if re.search(r'Author', ln, re.I):
            author_line = ln
        if re.search(r'Session', ln, re.I):
            session_line = ln
    return paper_id, title_line, author_line, session_line


# ── PDF builder ──────────────────────────────────────────────────────────────
def build_pdf(md_path, out_dir):
    filename = os.path.basename(md_path)
    stem     = re.sub(r'\.md$', '', filename, flags=re.I)
    # Strip "_whitepaper" suffix to normalize output names
    stem     = stem.replace("_whitepaper", "")
    out_path = os.path.join(out_dir, stem + ".pdf")

    try:
        with open(md_path, 'r', encoding='utf-8', errors='replace') as f:
            md_text = f.read()
    except Exception as e:
        print("  SKIP (read error): %s — %s" % (filename, e))
        return False

    paper_id, title, author, session = extract_meta(md_text, filename)

    story = []
    # Header block — use safe _p() fallbacks to avoid crashing on bad markup
    story.append(_sp(10))
    story.append(_p(paper_id, META))
    story.append(_p(_inline(title) if title else paper_id, TITLE))
    if author:
        story.append(_p(author, AUTH))
    if session:
        story.append(_p(session, META))
    story.append(_sp(6))
    story.append(_rule())

    # Content
    story.extend(md_to_story(md_text, paper_id))

    # Footer
    story.append(_sp(8))
    story.append(_rule())
    story.append(_p("%s | UQFF v4.74 | Star-Magic UQFF \xa9 2025 Daniel T. Murphy" % paper_id, FOOT))

    try:
        doc = SimpleDocTemplate(
            out_path, pagesize=A4,
            rightMargin=1.0*inch, leftMargin=1.0*inch,
            topMargin=0.9*inch,   bottomMargin=0.8*inch,
            title=paper_id, author="Daniel T. Murphy"
        )
        doc.build(story)
        return True
    except Exception as e:
        print("  ERROR building %s: %s" % (out_path, e))
        return False


# ── Collect all source files ─────────────────────────────────────────────────
def collect_sources():
    sources = []
    # 1. All whitepapers/ directory
    if os.path.isdir(WP_DIR):
        for fn in sorted(os.listdir(WP_DIR)):
            if fn.upper().endswith(".MD"):
                sources.append(os.path.join(WP_DIR, fn))
    # 2. Root PAPER_*.md files NOT already covered by whitepapers/
    wp_stems = set()
    for s in sources:
        stem = re.sub(r'\.md$', '', os.path.basename(s), flags=re.I)
        stem = stem.replace("_whitepaper", "")
        wp_stems.add(stem.upper())

    for fn in sorted(os.listdir(REPO_DIR)):
        if not fn.upper().endswith(".MD"):
            continue
        if not re.match(r'PAPER_\d+', fn, re.I):
            continue
        stem = re.sub(r'\.md$', '', fn, flags=re.I)
        stem = stem.replace("_whitepaper", "")
        if stem.upper() in wp_stems:
            continue   # already in whitepapers/
        sources.append(os.path.join(REPO_DIR, fn))

    return sources


# ── Main ─────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    sources = collect_sources()
    total  = len(sources)
    ok     = 0
    failed = 0

    print("=" * 60)
    print("Star-Magic UQFF Whitepaper PDF Builder")
    print("Output directory: %s" % OUT_DIR)
    print("Total papers to convert: %d" % total)
    print("=" * 60)

    for idx, src in enumerate(sources, 1):
        fn = os.path.basename(src)
        sys.stdout.write("[%d/%d] %s ... " % (idx, total, fn))
        sys.stdout.flush()
        try:
            result = build_pdf(src, OUT_DIR)
        except Exception as e:
            print("FAIL (%s)" % e)
            failed += 1
            continue
        if result:
            ok += 1
            print("OK")
        else:
            failed += 1

    print("=" * 60)
    print("Done.  Success: %d  |  Failed: %d  |  Total: %d" % (ok, failed, total))
    print("PDFs in: %s" % OUT_DIR)
