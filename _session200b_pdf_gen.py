#!/usr/bin/env python3
"""Generate PDFs for Session 200 PAPER_863-869 whitepapers using reportlab."""
import os, sys, re

from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer

BASE = os.path.dirname(os.path.abspath(__file__))
WP_DIR = os.path.join(BASE, "whitepapers")
PDF_DIR = os.path.join(BASE, "pdf")

PAPERS = [
    "PAPER_863_Water_Reactor_Birkeland_H2_Electrolysis_Efficiency",
    "PAPER_864_LRC_PseudoMonopole_SparkGap_Resonance",
    "PAPER_865_FieldGenerator_Spooky_NonLocal_TempDrop",
    "PAPER_866_DCEACE_Reversal_NdFeB_Caduceus_Motor",
    "PAPER_867_Mosquito_BioThermal_Efficiency_Benchmark",
    "PAPER_868_Topoconductor_Quantum_Cooling_Comparison",
    "PAPER_869_MilkyWay_82Day_StarTracking_UFT_Analysis",
]

styles = getSampleStyleSheet()
title_style = ParagraphStyle("Title2", parent=styles["Title"], fontSize=14, leading=18, spaceAfter=6)
h2_style = ParagraphStyle("H2", parent=styles["Heading2"], fontSize=12, leading=14, spaceAfter=4, spaceBefore=10)
body_style = ParagraphStyle("Body2", parent=styles["Normal"], fontSize=10, leading=13, spaceAfter=4)
code_style = ParagraphStyle("Code2", parent=styles["Code"], fontSize=9, leading=11, spaceAfter=2,
                            backColor="#f0f0f0", leftIndent=12)

def md_to_flowables(text):
    """Convert markdown text to reportlab flowables (simple parser)."""
    flowables = []
    for line in text.split("\n"):
        stripped = line.strip()
        if not stripped:
            flowables.append(Spacer(1, 6))
            continue
        # Escape XML entities
        safe = stripped.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
        if stripped.startswith("# "):
            flowables.append(Paragraph(safe[2:], title_style))
        elif stripped.startswith("## "):
            flowables.append(Paragraph(safe[3:], h2_style))
        elif stripped.startswith("- "):
            # Code lines inside list items
            content = safe[2:]
            if content.startswith("`") and content.endswith("`"):
                flowables.append(Paragraph(content[1:-1], code_style))
            else:
                flowables.append(Paragraph("&bull; " + content, body_style))
        elif stripped.startswith("**") or stripped.startswith("*"):
            flowables.append(Paragraph(safe, body_style))
        elif stripped.startswith("---"):
            flowables.append(Spacer(1, 8))
        elif stripped.startswith("```"):
            pass  # skip code fences
        else:
            flowables.append(Paragraph(safe, body_style))
    return flowables

ok = 0
for name in PAPERS:
    md_path = os.path.join(WP_DIR, name + ".md")
    pdf_path = os.path.join(PDF_DIR, name + ".pdf")
    try:
        with open(md_path, "r", encoding="utf-8") as f:
            txt = f.read()
        doc = SimpleDocTemplate(pdf_path, pagesize=letter,
                                leftMargin=0.75*inch, rightMargin=0.75*inch,
                                topMargin=0.75*inch, bottomMargin=0.75*inch)
        doc.build(md_to_flowables(txt))
        print(f"  OK: {name}.pdf")
        ok += 1
    except Exception as e:
        print(f"  FAIL: {name} -> {e}")

print(f"\n{ok}/7 PDFs generated")
