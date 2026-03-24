"""
generate_pdfs.py  —  UQFF PDF Generator (ReportLab edition)
============================================================
Converts:
  • Markdown papers     → nicely formatted PDF
  • C++ source files    → code-formatted PDF (monospace, syntax-highlighted)

Output directory: ./pdfs/
Usage: python generate_pdfs.py

Dependencies (all already installed):
  pip install reportlab markdown
"""

import os
import re
import sys
import textwrap
from pathlib import Path

import markdown

from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, HRFlowable,
    Preformatted, KeepTogether, PageBreak, Table, TableStyle
)
from reportlab.platypus.tableofcontents import TableOfContents
from reportlab.lib.colors import HexColor, black, white, Color

# ─── Color palette ────────────────────────────────────────────────────────────
UQFF_DARK       = HexColor('#0d1b2a')   # deep navy
UQFF_ACCENT     = HexColor('#1a6b8a')   # teal accent
UQFF_ACCENT2    = HexColor('#2d9cdb')   # lighter teal (h2)
UQFF_ACCENT3    = HexColor('#56b9da')   # h3
UQFF_BG_CODE    = HexColor('#eef4f8')   # code block background (light)
UQFF_TEXT       = HexColor('#1c2b36')   # body text
UQFF_MUTED      = HexColor('#4a6377')   # captions / footers
UQFF_RULE       = HexColor('#1a6b8a')   # horizontal rule
UQFF_TABLE_HDR  = HexColor('#1a6b8a')
UQFF_TABLE_ALT  = HexColor('#e8f4f8')
UQFF_CODE_TEXT  = HexColor('#1c2b36')   # monospace — dark on light bg

# ─── Page geometry ─────────────────────────────────────────────────────────────
PAGE_W, PAGE_H  = letter
MARGIN          = 0.85 * inch
INNER_W         = PAGE_W - 2 * MARGIN

# ─── Style factory ────────────────────────────────────────────────────────────

def make_styles():
    base = getSampleStyleSheet()

    def add(name, parent='Normal', **kw):
        s = ParagraphStyle(name, parent=base[parent], **kw)
        base.add(s)
        return s

    # Body
    add('UBody', fontName='Times-Roman', fontSize=10.5, leading=15,
        textColor=UQFF_TEXT, alignment=TA_JUSTIFY,
        spaceBefore=2, spaceAfter=4)

    # Headings
    add('UH1', fontName='Helvetica-Bold', fontSize=20, leading=26,
        textColor=UQFF_ACCENT, alignment=TA_LEFT,
        spaceBefore=18, spaceAfter=8)
    add('UH2', fontName='Helvetica-Bold', fontSize=15, leading=20,
        textColor=UQFF_ACCENT2, alignment=TA_LEFT,
        spaceBefore=14, spaceAfter=6)
    add('UH3', fontName='Helvetica-Bold', fontSize=12.5, leading=17,
        textColor=UQFF_ACCENT3, alignment=TA_LEFT,
        spaceBefore=10, spaceAfter=4)
    add('UH4', fontName='Helvetica-BoldOblique', fontSize=11, leading=15,
        textColor=UQFF_MUTED, alignment=TA_LEFT,
        spaceBefore=8, spaceAfter=3)

    # Title block
    add('UTitle', fontName='Helvetica-Bold', fontSize=24, leading=30,
        textColor=UQFF_ACCENT, alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=6)
    add('USubtitle', fontName='Helvetica-Oblique', fontSize=12, leading=16,
        textColor=UQFF_MUTED, alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=16)
    add('UMeta', fontName='Helvetica', fontSize=9.5, leading=13,
        textColor=UQFF_MUTED, alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=4)

    # Code / pre-formatted
    add('UCode', fontName='Courier', fontSize=8.5, leading=12,
        textColor=UQFF_CODE_TEXT, backColor=UQFF_BG_CODE,
        leftIndent=8, rightIndent=8, borderPadding=(5, 6, 5, 6),
        spaceBefore=6, spaceAfter=6)
    add('UCodeInline', fontName='Courier', fontSize=9,
        textColor=UQFF_ACCENT2, parent='UBody')

    # Lists
    add('UBullet', fontName='Times-Roman', fontSize=10.5, leading=14,
        textColor=UQFF_TEXT, leftIndent=18, firstLineIndent=-10,
        spaceBefore=1, spaceAfter=2)
    add('UNumbered', fontName='Times-Roman', fontSize=10.5, leading=14,
        textColor=UQFF_TEXT, leftIndent=22, firstLineIndent=-12,
        spaceBefore=1, spaceAfter=2)

    # Block quote / note
    add('UQuote', fontName='Times-Italic', fontSize=10, leading=14,
        textColor=UQFF_MUTED, leftIndent=24, rightIndent=8,
        spaceBefore=6, spaceAfter=6, borderPadding=(4, 0, 4, 8),
        borderColor=UQFF_ACCENT, borderWidth=2)

    # Footer
    add('UFooter', fontName='Helvetica', fontSize=8, leading=10,
        textColor=UQFF_MUTED, alignment=TA_CENTER)

    # C++ source file title
    add('SrcTitle', fontName='Helvetica-Bold', fontSize=18, leading=24,
        textColor=UQFF_ACCENT, alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=4)
    add('SrcMeta', fontName='Courier', fontSize=9, leading=13,
        textColor=UQFF_MUTED, alignment=TA_CENTER,
        spaceBefore=0, spaceAfter=12)

    return base


# ─── Page template (header / footer) ──────────────────────────────────────────

def make_doc(output_path: Path, title: str):
    def on_page(canvas, doc, title=title):
        canvas.saveState()
        # Header rule
        canvas.setStrokeColor(UQFF_RULE)
        canvas.setLineWidth(0.8)
        canvas.line(MARGIN, PAGE_H - 0.55*inch, PAGE_W - MARGIN, PAGE_H - 0.55*inch)
        # Header text
        canvas.setFont('Helvetica', 8)
        canvas.setFillColor(UQFF_MUTED)
        canvas.drawString(MARGIN, PAGE_H - 0.48*inch, 'UQFF — Daniel T. Murphy')
        safe_title = title[:80]
        canvas.drawRightString(PAGE_W - MARGIN, PAGE_H - 0.48*inch, safe_title)
        # Footer rule
        canvas.line(MARGIN, 0.55*inch, PAGE_W - MARGIN, 0.55*inch)
        canvas.drawCentredString(PAGE_W / 2, 0.38*inch,
            f'Page {doc.page}  •  Unified Quantum Field Framework (UQFF) v4.81')
        canvas.restoreState()

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=letter,
        leftMargin=MARGIN, rightMargin=MARGIN,
        topMargin=0.75*inch, bottomMargin=0.75*inch,
        title=title,
        author='Daniel T. Murphy — UQFF Research Consortium',
        subject='Unified Quantum Field Framework',
    )
    doc.build_kwargs = {'onFirstPage': on_page, 'onLaterPages': on_page}
    return doc


# ─── Markdown → ReportLab flowables ───────────────────────────────────────────

# Escape XML special chars for Paragraph
_XML_ESC = str.maketrans({'&': '&amp;', '<': '&lt;', '>': '&gt;',
                           '"': '&quot;', "'": '&#39;'})

def xml_escape(s: str) -> str:
    return s.translate(_XML_ESC)

def inline_markup(raw_text: str, styles) -> str:
    """Convert **bold**, *italic*, `code` inline to RL XML.

    Processing order:
      1. Extract backtick code spans (protect them from being italic-matched)
      2. XML-escape the non-code content
      3. Apply **bold** and *italic* (underscore-italic intentionally OMITTED to
         avoid false positives on technical identifiers, filenames, LaTeX)
      4. Restore code spans as coloured Courier font
    """
    # ── Step 1: save backtick spans before any other processing ───────────────
    code_spans: list[str] = []

    def save_code(m: re.Match) -> str:
        idx = len(code_spans)
        code_spans.append(m.group(1))
        return f'\x01CODE{idx:04d}\x01'   # non-XML placeholder

    text = re.sub(r'`([^`\n]+)`', save_code, raw_text)

    # ── Step 2: XML-escape everything else ────────────────────────────────────
    text = text.translate(_XML_ESC)

    # ── Step 3: bold+italic inline markup (NO underscore variants) ────────────
    text = re.sub(r'\*\*\*(.+?)\*\*\*', r'<b><i>\1</i></b>', text)
    text = re.sub(r'\*\*(.+?)\*\*',     r'<b>\1</b>', text)
    # *italic* — only when * is not adjacent to another * or a word character
    text = re.sub(r'(?<![*\w])\*([^*\n]+?)\*(?![*\w])', r'<i>\1</i>', text)

    # ── Step 4: restore code spans as coloured Courier ────────────────────────
    for idx, span in enumerate(code_spans):
        safe = span.translate(_XML_ESC)
        text = text.replace(
            f'\x01CODE{idx:04d}\x01',
            f'<font name="Courier" color="#2d9cdb">{safe}</font>'
        )

    return text


def md_to_flowables(md_text: str, styles) -> list:
    """Parse Markdown text into ReportLab flowables."""
    flowables = []
    lines = md_text.splitlines()
    i = 0
    in_code  = False
    code_buf = []
    table_buf = []
    in_table = False

    def flush_code():
        nonlocal in_code, code_buf
        if code_buf:
            raw = '\n'.join(code_buf)
            # Wrap long lines
            wrapped = []
            for ln in raw.splitlines():
                if len(ln) > 110:
                    for chunk in textwrap.wrap(ln, 110, subsequent_indent='    ', break_long_words=True):
                        wrapped.append(chunk)
                else:
                    wrapped.append(ln)
            wrapped_text = xml_escape('\n'.join(wrapped))
            flowables.append(Preformatted(wrapped_text, styles['UCode']))
        code_buf.clear()
        in_code = False

    def flush_table():
        nonlocal in_table, table_buf
        if not table_buf:
            in_table = False
            return
        # Parse header + separator + rows
        rows = []
        for r_idx, row_line in enumerate(table_buf):
            # Skip separator row (---|---)
            if re.match(r'^[\s|:-]+$', row_line.replace('|', '').strip()):
                continue
            cells = [c.strip() for c in row_line.strip().strip('|').split('|')]
            rows.append(cells)
        if not rows:
            table_buf.clear()
            in_table = False
            return
        max_cols = max(len(r) for r in rows)
        # Pad rows to equal cols
        padded = [r + [''] * (max_cols - len(r)) for r in rows]
        col_w = INNER_W / max_cols if max_cols else INNER_W

        tbl_data = []
        for r_idx, row in enumerate(padded):
            if r_idx == 0:
                hdr_style = ParagraphStyle('TblHdr', parent=styles['UBody'],
                                           fontName='Helvetica-Bold', textColor=white)
                tbl_data.append([
                    Paragraph(inline_markup(c, styles), hdr_style)
                    for c in row
                ])
            else:
                tbl_data.append([
                    Paragraph(inline_markup(c, styles), styles['UBody'])
                    for c in row
                ])
        col_widths = [col_w] * max_cols
        t = Table(tbl_data, colWidths=col_widths, repeatRows=1)
        t.setStyle(TableStyle([
            ('BACKGROUND',  (0,0), (-1,0), UQFF_TABLE_HDR),
            ('TEXTCOLOR',   (0,0), (-1,0), white),
            ('FONTNAME',    (0,0), (-1,0), 'Helvetica-Bold'),
            ('FONTSIZE',    (0,0), (-1,0), 9),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, UQFF_TABLE_ALT]),
            ('FONTSIZE',    (0,1), (-1,-1), 9),
            ('GRID',        (0,0), (-1,-1), 0.4, UQFF_RULE),
            ('VALIGN',      (0,0), (-1,-1), 'TOP'),
            ('LEFTPADDING', (0,0), (-1,-1), 5),
            ('RIGHTPADDING',(0,0), (-1,-1), 5),
            ('TOPPADDING',  (0,0), (-1,-1), 3),
            ('BOTTOMPADDING',(0,0),(-1,-1), 3),
        ]))
        flowables.append(t)
        flowables.append(Spacer(1, 6))
        table_buf.clear()
        in_table = False

    while i < len(lines):
        line = lines[i]

        # ── Fenced code block ─────────────────────────────────────────────────
        if line.strip().startswith('```'):
            if in_code:
                flush_code()
            else:
                if in_table: flush_table()
                in_code = True
            i += 1
            continue

        if in_code:
            code_buf.append(line)
            i += 1
            continue

        # ── Table detection ───────────────────────────────────────────────────
        if '|' in line and line.strip().startswith('|'):
            if not in_table:
                in_table = True
            table_buf.append(line)
            i += 1
            continue
        elif in_table:
            flush_table()

        # ── Blank line ────────────────────────────────────────────────────────
        if not line.strip():
            flowables.append(Spacer(1, 4))
            i += 1
            continue

        # ── Headings ──────────────────────────────────────────────────────────
        m = re.match(r'^(#{1,6})\s+(.*)', line)
        if m:
            level = len(m.group(1))
            text  = inline_markup(m.group(2), styles)
            style = {1:'UH1', 2:'UH2', 3:'UH3', 4:'UH4'}.get(level, 'UH4')
            if level == 1:
                flowables.append(HRFlowable(width='100%', thickness=1.2,
                                            color=UQFF_RULE, spaceAfter=4))
            flowables.append(Paragraph(text, styles[style]))
            if level <= 2:
                flowables.append(HRFlowable(width='100%', thickness=0.5,
                                            color=UQFF_RULE, spaceAfter=2))
            i += 1
            continue

        # ── Horizontal rule ───────────────────────────────────────────────────
        if re.match(r'^[-*_]{3,}\s*$', line):
            flowables.append(HRFlowable(width='100%', thickness=0.8,
                                        color=UQFF_RULE, spaceBefore=6, spaceAfter=6))
            i += 1
            continue

        # ── Block quote ───────────────────────────────────────────────────────
        if line.startswith('>'):
            text = line.lstrip('> ').strip()
            flowables.append(Paragraph(inline_markup(text, styles), styles['UQuote']))
            i += 1
            continue

        # ── Unordered list ────────────────────────────────────────────────────
        if re.match(r'^[\s]*[-*+]\s+', line):
            text = re.sub(r'^[\s]*[-*+]\s+', '', line)
            text = inline_markup(text, styles)
            flowables.append(Paragraph('•  ' + text, styles['UBullet']))
            i += 1
            continue

        # ── Ordered list ──────────────────────────────────────────────────────
        if re.match(r'^[\s]*\d+\.\s+', line):
            m2 = re.match(r'^[\s]*(\d+)\.\s+(.*)', line)
            if m2:
                text = inline_markup(m2.group(2), styles)
                flowables.append(Paragraph(f'{m2.group(1)}.  {text}', styles['UNumbered']))
            i += 1
            continue

        # ── Regular paragraph ─────────────────────────────────────────────────
        # Accumulate continuation lines
        para_lines = [line]
        while i + 1 < len(lines):
            nxt = lines[i + 1]
            if (not nxt.strip() or nxt.startswith('#') or nxt.startswith('```')
                    or nxt.startswith('>') or ('|' in nxt and nxt.strip().startswith('|'))
                    or re.match(r'^[\s]*[-*+]\s+', nxt)
                    or re.match(r'^[\s]*\d+\.\s+', nxt)
                    or re.match(r'^[-*_]{3,}\s*$', nxt)):
                break
            para_lines.append(nxt)
            i += 1
        text = ' '.join(para_lines)
        text = inline_markup(text, styles)
        flowables.append(Paragraph(text, styles['UBody']))
        i += 1

    if in_code:
        flush_code()
    if in_table:
        flush_table()

    return flowables


# ─── Title block extractor ─────────────────────────────────────────────────────

def extract_title_meta(md_text: str):
    """Try to extract H1 title from first heading or first line."""
    for line in md_text.splitlines():
        m = re.match(r'^#\s+(.*)', line)
        if m:
            return m.group(1).strip()
    return 'UQFF Research Paper'


# ─── Markdown → PDF ───────────────────────────────────────────────────────────

def md_file_to_pdf(md_path: Path, out_dir: Path, styles):
    md_text = md_path.read_text(encoding='utf-8')
    paper_title = extract_title_meta(md_text)
    out_path = out_dir / (md_path.stem + '.pdf')

    doc = make_doc(out_path, paper_title)

    # Build title block
    story = []
    story.append(Spacer(1, 0.3 * inch))
    story.append(Paragraph(xml_escape(paper_title), styles['UTitle']))
    story.append(Paragraph('Unified Quantum Field Framework (UQFF)', styles['USubtitle']))
    story.append(Paragraph('Daniel T. Murphy  •  UQFF Research Consortium  •  2026', styles['UMeta']))
    story.append(Paragraph(f'Source: {md_path.name}', styles['UMeta']))
    story.append(HRFlowable(width='80%', thickness=1.5, color=UQFF_ACCENT,
                             spaceBefore=8, spaceAfter=16, hAlign='CENTER'))

    # Convert body
    story += md_to_flowables(md_text, styles)
    story.append(Spacer(1, 0.5 * inch))

    doc.build(story, **doc.build_kwargs)
    return out_path


# ─── C++ source → PDF ─────────────────────────────────────────────────────────

# Very lightweight keyword highlighter — colour keywords without HTML
_CPP_KW = re.compile(
    r'\b(auto|bool|break|case|catch|class|const|constexpr|continue|default|'
    r'delete|do|double|else|enum|explicit|extern|false|float|for|friend|if|'
    r'inline|int|long|namespace|new|nullptr|operator|private|protected|'
    r'public|return|short|signed|sizeof|static|struct|switch|template|this|'
    r'throw|true|try|typedef|typename|union|unsigned|using|virtual|void|'
    r'volatile|while|#include|#ifndef|#define|#endif|#pragma|std::)\b'
)

def cpp_line_markup(line: str) -> tuple[str, str]:
    """Return (safe_line, colour) for line, wrapping keywords in font tags."""
    safe = xml_escape(line)
    # Colour single-line comments
    if '//' in safe:
        idx = safe.index('//')
        comment = safe[idx:]
        code    = safe[:idx]
        code    = _CPP_KW.sub(r'<font color="#1a6b8a"><b>\1</b></font>', code)
        safe    = code + f'<font color="#4a7a60">{comment}</font>'
    else:
        safe = _CPP_KW.sub(r'<font color="#1a6b8a"><b>\1</b></font>', safe)
    return safe


def cpp_file_to_pdf(cpp_path: Path, out_dir: Path, styles):
    src_text = cpp_path.read_text(encoding='utf-8', errors='replace')
    out_path = out_dir / (cpp_path.stem + '.pdf')

    # Extract description comment block from top of file
    desc_lines = []
    for ln in src_text.splitlines()[:60]:
        if ln.strip().startswith('//') or ln.strip().startswith('*') or ln.strip().startswith('/*'):
            cl = re.sub(r'^\s*(/\*+|\*+/?|//+)\s?', '', ln).strip()
            if cl:
                desc_lines.append(cl)
        elif ln.strip() and not ln.strip().startswith('#') and desc_lines:
            break

    doc = make_doc(out_path, cpp_path.name)
    story = []

    # Title block
    story.append(Spacer(1, 0.3 * inch))
    story.append(Paragraph(xml_escape(cpp_path.name), styles['SrcTitle']))
    story.append(Paragraph('UQFF C++ Source Module', styles['USubtitle']))
    if desc_lines:
        desc_text = xml_escape(' '.join(desc_lines[:8]))
        story.append(Paragraph(desc_text, styles['SrcMeta']))
    story.append(Paragraph('Daniel T. Murphy  •  UQFF Research Consortium  •  2026', styles['UMeta']))
    story.append(Paragraph(f'Lines: {len(src_text.splitlines()):,}  '
                            f'Size: {len(src_text.encode()):,} bytes', styles['UMeta']))
    story.append(HRFlowable(width='80%', thickness=1.5, color=UQFF_ACCENT,
                             spaceBefore=8, spaceAfter=14, hAlign='CENTER'))

    # Description summary if present
    if desc_lines:
        story.append(Paragraph('<b>Module Description</b>', styles['UH2']))
        for dl in desc_lines:
            escaped = xml_escape(dl)
            story.append(Paragraph(escaped, styles['UBody']))
        story.append(Spacer(1, 10))
        story.append(HRFlowable(width='100%', thickness=0.5, color=UQFF_RULE, spaceAfter=6))

    # Source code — split into pages of ~60 lines each to avoid memory issues
    all_lines  = src_text.splitlines()
    total_lines = len(all_lines)
    CHUNK = 55   # lines per code block

    story.append(Paragraph('<b>Source Code</b>', styles['UH2']))

    for chunk_start in range(0, total_lines, CHUNK):
        chunk = all_lines[chunk_start:chunk_start + CHUNK]
        # Build a pre-formatted block with line numbers
        numbered = []
        for rel, ln in enumerate(chunk, 1):
            abs_ln = chunk_start + rel
            # Truncate very long lines
            if len(ln) > 115:
                ln = ln[:112] + '...'
            numbered.append(f'{abs_ln:5d}  {ln}')
        block_text = xml_escape('\n'.join(numbered))
        story.append(Preformatted(block_text, styles['UCode']))

    story.append(Spacer(1, 0.4 * inch))
    story.append(Paragraph(f'— End of {cpp_path.name} —', styles['UMeta']))

    doc.build(story, **doc.build_kwargs)
    return out_path


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    repo = Path(__file__).parent
    out_dir = repo / 'pdf'
    out_dir.mkdir(parents=True, exist_ok=True)

    styles = make_styles()

    # Files to convert
    md_targets = [
        'PAPER_422_UQFF_29System_CrossValidation_Matrix.md',
        'PAPER_423_Um_Complete_SSq_Vacuum_Thermal_Damping.md',
        'PAPER_424_FUBii_Um_Universal_Companion_Catalog.md',
        'PAPER_425_DPM_Four_Component_Correlation.md',
        'PAPER_426_UA_SCm_JWST_ALMA_CERN_Validation_Table.md',
        'PAPER_427_26D_Resonance_Layer_Amplitude_Frequency.md',
        'PAPER_428_HRes_Periodic_Table_Universal_Nuclear.md',
        'PAPER_429_Three_New_Number_Systems_Vacuum_Dipole_Buoyancy.md',
    ]

    cpp_targets = [
        'UQFFLearningAssessment_001.cpp',
        'UQFFLearningAssessment_002.cpp',
    ]

    print('=' * 70)
    print('  UQFF PDF Generator  —  ReportLab edition')
    print(f'  Output directory: {out_dir}')
    print('=' * 70)

    errors = []

    print('\n── Markdown papers ──────────────────────────────────────────────────')
    for fname in md_targets:
        src = repo / 'whitepapers' / fname
        if not src.exists():
            print(f'  [SKIP]  {fname}  (not found)')
            errors.append(f'NOT FOUND: {fname}')
            continue
        try:
            out = md_file_to_pdf(src, out_dir, styles)
            size_kb = out.stat().st_size / 1024
            print(f'  [OK]    {out.name:60s}  ({size_kb:6.1f} KB)')
        except Exception as exc:
            print(f'  [FAIL]  {fname}: {exc}')
            errors.append(f'FAILED: {fname}: {exc}')

    print('\n── C++ learning assessments ─────────────────────────────────────────')
    for fname in cpp_targets:
        src = repo / fname
        if not src.exists():
            print(f'  [SKIP]  {fname}  (not found)')
            errors.append(f'NOT FOUND: {fname}')
            continue
        try:
            out = cpp_file_to_pdf(src, out_dir, styles)
            size_kb = out.stat().st_size / 1024
            print(f'  [OK]    {out.name:60s}  ({size_kb:6.1f} KB)')
        except Exception as exc:
            print(f'  [FAIL]  {fname}: {exc}')
            errors.append(f'FAILED: {fname}: {exc}')

    print('\n' + '=' * 70)
    if errors:
        print(f'  {len(errors)} error(s):')
        for e in errors:
            print(f'    {e}')
    else:
        total = len(md_targets) + len(cpp_targets)
        print(f'  All {total} PDFs generated successfully.')
    print('=' * 70)

    # List generated PDFs
    pdfs = sorted(out_dir.glob('*.pdf'))
    print(f'\nGenerated PDFs in {out_dir}:\n')
    total_size = 0
    for p in pdfs:
        sz = p.stat().st_size
        total_size += sz
        print(f'  {p.name:65s}  {sz/1024:7.1f} KB')
    print(f'\n  Total: {len(pdfs)} files, {total_size/1024:.1f} KB')
    return 0 if not errors else 1


if __name__ == '__main__':
    sys.exit(main())
