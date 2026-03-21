# -*- coding: utf-8 -*-
"""
build_papers_371_376.py
Generates PDFs for PAPER_371 through PAPER_376
Requires: pip install reportlab
"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import inch
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, HRFlowable
)
import os

DIR = os.path.dirname(os.path.abspath(__file__))
W, H = A4
TW = W - 2.2 * inch

DARK  = colors.HexColor("#1a1a2e")
BLUE  = colors.HexColor("#16213e")
RULE  = colors.HexColor("#0f3460")
LGREY = colors.HexColor("#f5f5f5")
TBLHD = colors.HexColor("#d0e4f7")
HBOX  = colors.HexColor("#eef4ff")

def sty(name, font="Times-Roman", size=11, lead=16, align=TA_LEFT,
        color=colors.black, sb=6, sa=4, li=0, fl=0, bg=None, bw=0, bp=0):
    kw = dict(fontName=font, fontSize=size, leading=lead, alignment=align,
              textColor=color, spaceBefore=sb, spaceAfter=sa,
              leftIndent=li, firstLineIndent=fl)
    if bg: kw["backColor"] = bg
    if bw:
        kw["borderColor"] = RULE; kw["borderWidth"] = bw; kw["borderPad"] = bp
    return ParagraphStyle(name, **kw)

TTS  = sty("TT","Times-Bold",18,24,TA_CENTER,DARK,sb=0,sa=6)
AUS  = sty("AU","Times-Roman",12,16,TA_CENTER,DARK,sb=0,sa=3)
MES  = sty("ME","Times-Italic",10,14,TA_CENTER,colors.grey,sb=0,sa=2)
ALS  = sty("AL","Times-Bold",11,14,TA_CENTER,DARK,sb=6,sa=4)
ABS  = sty("AB","Times-Roman",10,14,TA_JUSTIFY,colors.black,li=28,sa=6)
H1S  = sty("H1","Times-Bold",14,18,TA_LEFT,DARK,sb=14,sa=4)
H2S  = sty("H2","Times-Bold",12,16,TA_LEFT,BLUE,sb=10,sa=3)
H3S  = sty("H3","Times-BoldItalic",11,14,TA_LEFT,DARK,sb=8,sa=3)
BOS  = sty("BO","Times-Roman",11,16,TA_JUSTIFY,colors.black,sa=6)
EQS  = sty("EQ","Courier",9,12,TA_LEFT,colors.black,li=20,sa=3,sb=3,bg=LGREY)
BUS  = sty("BU","Times-Roman",11,15,TA_LEFT,colors.black,li=24,fl=-12,sa=3)
RES  = sty("RE","Times-Roman",10,14,TA_LEFT,colors.black,li=24,fl=-24,sa=3)
CAS  = sty("CA","Times-Italic",9,12,TA_CENTER,colors.grey,sa=8)

def B(t):   return "<b>%s</b>" % t
def I(t):   return "<i>%s</i>" % t
def TT(t):  return "<font face='Courier' size='9'>%s</font>" % t
def SUP(t): return "<super>%s</super>" % t
def SUB(t): return "<sub>%s</sub>" % t
def p(text, style=None): return Paragraph(text, style or BOS)
def rule(): return HRFlowable(width="100%", thickness=0.6, color=RULE, spaceAfter=4, spaceBefore=4)
def eq(text): return Paragraph(text, EQS)
def SP(n=6): return Spacer(1, n)

def tbl(rows, widths=None, caption=None):
    t = Table(rows, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),TBLHD),("FONTNAME",(0,0),(-1,0),"Times-Bold"),
        ("FONTSIZE",(0,0),(-1,-1),9),("FONTNAME",(0,1),(-1,-1),"Times-Roman"),
        ("GRID",(0,0),(-1,-1),0.4,colors.grey),
        ("ROWBACKGROUNDS",(0,1),(-1,-1),[colors.white,LGREY]),
        ("ALIGN",(0,0),(-1,-1),"CENTER"),("VALIGN",(0,0),(-1,-1),"MIDDLE"),
        ("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4),
        ("LEFTPADDING",(0,0),(-1,-1),6),("RIGHTPADDING",(0,0),(-1,-1),6),
    ]))
    out = [SP(4), t]
    if caption: out.append(p(caption, CAS))
    out.append(SP(8))
    return out

def title_block(num, title, session, src):
    return [
        SP(15),
        p(num, MES),
        p(title, TTS),
        SP(6),
        p("Daniel T. Murphy", AUS),
        p("Independent Research &#8226; Star-Magic / UQFF Project", MES),
        p("Session: " + session + " &#8226; Source: " + src + " &#8226; UQFF v4.74", MES),
        SP(8), rule()
    ]

def footer(num, session):
    return [
        SP(10), rule(),
        p(num + " | UQFF v4.74 (Session 115) | Commit: 870cb4f | Session: "
          + session + " | Author: Daniel T. Murphy | &#169;2025", MES)
    ]

def build(out_path, story):
    doc = SimpleDocTemplate(out_path, pagesize=A4,
        rightMargin=1.1*inch, leftMargin=1.1*inch,
        topMargin=1.0*inch, bottomMargin=1.0*inch)
    doc.build(story)
    print("PDF written ->", out_path)


# ═══════════════════════════════════════════════════════
# PAPER_371 — MUGE 12-Term Superconductive Resonance
# ═══════════════════════════════════════════════════════
def build_371():
    out = os.path.join(DIR, "PAPER_371_MUGE_12Term_Superconductive_Resonance.pdf")
    s = title_block("PAPER_371",
        "MUGE 12-Term Superconductive Resonance Framework",
        "101", "grok_share_11254865.txt (lines 2000-2700)")

    s += [p("Abstract", ALS),
          p("This paper presents the complete 12-term MUGE (Modified Unified Gravity Equation) "
            "Superconductive Resonance Framework derived from the Star Magic UQFF formalism. "
            "The framework extends the standard UQFF by introducing twelve independently computable "
            "acceleration terms rooted in vacuum energy, THz resonance phenomena, aether coupling, "
            "and superconductive quantum effects. Validated against seven astrophysical systems "
            "including magnetar SGR1745-2900 and Sagittarius A*.", ABS),
          SP(6), rule()]

    s.append(p("1.  Master Equation", H1S))
    s.append(p("The MUGE Resonance master acceleration sums twelve independent terms:"))
    s.append(eq(
        "g(r,t) = a_DPM + a_THz + a_vac,diff + a_super,freq + a_aether,res + Ug4i"
    ))
    s.append(eq(
        "       + a_quantum,freq + a_Aether,freq + a_fluid,freq + a_osc + a_exp,freq + f_TRZ"
    ))
    s.append(SP(8))

    s.append(p("2.  Term Definitions", H1S))

    terms = [
        ("2.1", "DPM Acceleration", "a_DPM",
         "F_DPM = I * A * (omega_1 - omega_2)\n"
         "a_DPM = F_DPM * f_DPM * E_vac,neb * c * V_sys",
         "Differential Pressure Motor term coupling EM frequency differences to vacuum energy and system volume."),
        ("2.2", "THz Coupling", "a_THz",
         "a_THz = f_THz * E_vac,neb * v_exp * a_DPM / (E_vac,ISM * c)",
         "THz-frequency coupling between nebular and ISM vacuum energy mediating energy transfer."),
        ("2.3", "Vacuum Energy Differential", "a_vac,diff",
         "a_vac,diff = dE_vac * v_exp^2 * a_DPM / (E_vac,neb * c^2)",
         "Second-order velocity coupling to the vacuum energy differential between nebular and ISM."),
        ("2.4", "Superconductive Frequency", "a_super,freq",
         "a_super,freq = F_super * f_THz * a_DPM / (E_vac,neb * c)",
         "Superconductive force (F_super = 6.287e-19 N) coupling at THz resonance frequency."),
        ("2.5", "Aether Resonance", "a_aether,res",
         "a_aether,res = U_A,SCM * omega_i * f_THz * a_DPM * (1 + f_TRZ)",
         "Aether-superconducting manifold coupling modulated by THz frequency and TRZ correction."),
        ("2.6", "Reactive Ug4 Coupling", "Ug4i",
         "Ug4i = k4_res * E_react(t) * f_react * a_DPM / E_vac,neb * c",
         "Reactive gravitational coupling with time-dependent reaction energy E_react(t)."),
        ("2.7", "Quantum Frequency", "a_quantum,freq",
         "a_quantum,freq = f_quantum * E_vac,neb * a_DPM / (E_vac,ISM * c)",
         "Resonance at f_quantum = 1.445e-17 Hz = 2pi/t_Hubble (Hubble time resonance frequency)."),
        ("2.8", "Aether Frequency", "a_Aether,freq",
         "a_Aether,freq = f_Aether * E_vac,neb * a_DPM / (E_vac,ISM * c)",
         "Aether coupling at f_Aether = 1.576e-35 Hz (sub-Planckian aether frequency)."),
        ("2.9", "Fluid Frequency", "a_fluid,freq",
         "a_fluid,freq = f_fluid * E_vac,neb * V_sys / (E_vac,ISM * c)",
         "Fluid dynamical term coupling system volume to vacuum energy. "
         "Dominant for SGR1745-2900: a_fluid,freq = 1.773e-9 m/s2."),
        ("2.10", "Oscillation Term", "a_osc",
         "a_osc = f_osc * cos(2 * pi * f_osc * t)",
         "Temporal oscillation at f_osc = 4.57e14 Hz (visible-light frequency range)."),
        ("2.11", "Expansion Frequency", "a_exp,freq",
         "a_exp,freq = 2*pi * H_z * t * E_vac,neb * a_DPM / (E_vac,ISM * c)",
         "Time-dependent expansion coupling proportional to Hubble parameter H_z."),
        ("2.12", "Time-Reversal Correction", "f_TRZ",
         "f_TRZ = 0.1   (constant additive term)",
         "Topological resonance zone correction providing a baseline 0.1 m/s2 additive contribution."),
    ]
    for (num, name, symbol, formula, desc) in terms:
        s.append(p("%s  %s  (%s)" % (num, B(name), I(symbol)), H2S))
        for line in formula.split("\n"):
            s.append(eq(line))
        s.append(SP(3))
        s.append(p(desc))

    s.append(p("3.  Canonical Parameter Values (ResonanceParams)", H1S))
    s += tbl([
        ["Symbol", "Value", "Units", "Description"],
        ["f_DPM", "1e12", "Hz", "DPM frequency"],
        ["f_THz", "1e12", "Hz", "THz coupling frequency"],
        ["E_vac,neb", "7.09e-36", "J", "Nebular vacuum energy"],
        ["E_vac,ISM", "7.09e-37", "J", "ISM vacuum energy"],
        ["dE_vac", "6.381e-36", "J", "Vacuum energy differential"],
        ["F_super", "6.287e-19", "N", "Superconductive force"],
        ["U_A,SCM", "10", "--", "Aether SCm coupling"],
        ["omega_i", "1e-8", "rad/s", "Intrinsic angular frequency"],
        ["k4_res", "1.0", "--", "Resonance Ug4 coupling"],
        ["f_react", "1e10", "Hz", "Reactive frequency"],
        ["f_quantum", "1.445e-17", "Hz", "Hubble time resonance = 2pi/t_H"],
        ["f_Aether", "1.576e-35", "Hz", "Sub-Planckian aether frequency"],
        ["f_osc", "4.57e14", "Hz", "Oscillation (visible-light) frequency"],
        ["f_TRZ", "0.1", "--", "Time-reversal correction (constant)"],
    ], widths=[TW*0.18, TW*0.18, TW*0.14, TW*0.50],
       caption="Table 1. ResonanceParams canonical default values")

    s.append(p("4.  Validation: Expected Unit Test Values", H1S))
    s += tbl([
        ["Test", "System", "Expected Value"],
        ["a_fluid_freq", "SGR1745-2900", "1.773e-9 m/s" + SUP("2")],
        ["resonance_MUGE", "SGR1745-2900", "1.773e-9 m/s" + SUP("2")],
        ["a_THz (a_DPM=3.545e-42, v_exp=1e3)", "--", "1.182e-33"],
        ["a_vac_diff", "--", "3.545e-53"],
        ["a_super_freq", "--", "1.048e-21"],
        ["a_aether_res", "--", "3.900e-38"],
        ["a_quantum_freq", "--", "1.708e-66"],
        ["a_Aether_freq", "--", "1.863e-84"],
        ["a_exp_freq (t=3.799e10)", "--", "1.623e-57"],
    ], widths=[TW*0.48, TW*0.22, TW*0.30],
       caption="Table 2. Unit test expected values for MUGE resonance terms")

    s.append(p("5.  Implementation", H1S))
    s.append(p(B("C++:") + " " + TT("STAR_MAGIC_09SEPT_UQFF_MODULE.cpp")
               + ", namespace " + TT("StarMagic09Sept_Session101")))
    for fn in ["compute_aDPM()", "compute_aTHz()", "compute_avac_diff()",
               "compute_asuper_freq()", "compute_aaether_res()", "compute_Ug4i()",
               "compute_aquantum_freq()", "compute_aAether_freq()", "compute_afluid_freq()",
               "compute_Osc_term()", "compute_aexp_freq()", "compute_resonance_MUGE()"]:
        s.append(p("&#8226;  " + TT(fn), BUS))
    s.append(p(B("Python:") + " " + TT("CondensedPhysics4.py")
               + ", class " + TT("MUGESuperconductive12TermResonanceCalculator") + " (CP4 #19)"))
    s.append(p(B("Wolfram macro:") + " " + TT("WOLFRAM_TERM_MUGE_RESONANCE")))

    s += footer("PAPER_371", "101")
    build(out, s)


# ═══════════════════════════════════════════════════════
# PAPER_372 — Compressed UQFF with B/Bcrit Superconductivity
# ═══════════════════════════════════════════════════════
def build_372():
    out = os.path.join(DIR, "PAPER_372_Compressed_UQFF_Bcrit.pdf")
    s = title_block("PAPER_372",
        "Compressed UQFF with B/B" + SUB("crit") + " Superconductivity",
        "101", "grok_share_11254865.txt (lines 2700-3400)")

    s += [p("Abstract", ALS),
          p("This paper presents the Compressed UQFF formulation incorporating Newtonian gravity, "
            "Hubble expansion, superconductivity via the B/B" + SUB("crit") + " flux-quenching "
            "factor, environmental coupling, cosmological constant contribution, quantum coherence, "
            "fluid dynamics, and dark matter perturbation. Parameterised for seven astrophysical "
            "systems and validated against SGR1745-2900. This is the first UQFF formulation to "
            "explicitly incorporate Bardeen-Cooper-Schrieffer (BCS) superconductivity quenching "
            "into the gravitational coupling (Linear Meissner form; see PAPER_375 for the "
            "exponential form).", ABS),
          SP(6), rule()]

    s.append(p("1.  Master Equation", H1S))
    s.append(p("The Compressed UQFF master gravity equation combines eight modular terms:"))
    s.append(eq("g_UQFF(r,t) = GM(t)/r^2 * [1 + H(t,z)] * [1 - B/B_crit] * [1 + F_env]"))
    s.append(eq("            + (Ug1 + Ug2 + Ug3' + Ug4) + Lambda*c^2/3"))
    s.append(eq("            + (hbar/dxdp) * integral(psi* H_hat psi dV) * 2pi/t_Hubble"))
    s.append(eq("            + rho_fluid*V*g + (M_vis+M_DM)*(drho/rho + 3GM/r^3)"))
    s.append(p("where H(t,z) = H" + SUB("0") + "t, H" + SUB("0")
               + " = 2.269&#215;10" + SUP("&#8722;18") + " s" + SUP("&#8722;1")
               + " (Newtonian cosmological expansion approximation)."))

    s.append(p("2.  Modular Component Functions", H1S))
    s += tbl([
        ["Function", "Formula", "Key Constants"],
        ["compressed_base", "G*M/r^2", "G = 6.674e-11"],
        ["compressed_expansion", "1 + H0*t", "H0 = 2.269e-18 s-1"],
        ["compressed_super_adj", "1 - B/B_crit", "Linear Meissner form"],
        ["compressed_env", "1.0", "Default (no environment coupling)"],
        ["compressed_cosm", "Lambda*c^2/3", "Lambda = 1.1e-52 m-2"],
        ["compressed_quantum", "(hbar/1e-68)*2.176e-18*(2pi/t_H)", "t_H = 4.35e17 s"],
        ["compressed_fluid", "rho_f * V * g_local", "from MUGESystem"],
        ["compressed_perturbation", "(M+M_DM)*(drho/rho + 3GM/r^3)", "drho/rho = 1e-5"],
    ], widths=[TW*0.30, TW*0.42, TW*0.28],
       caption="Table 1. Compressed UQFF modular component functions")

    s.append(p("3.  Seven Astrophysical System Parameters", H1S))
    s += tbl([
        ["System", "M (kg)", "r (m)", "B (T)", "B_crit (T)", "V_sys (m3)"],
        ["SGR1745-2900", "2.984e30", "1e4", "1e10", "1e11", "4.189e12"],
        ["Sagittarius A*", "8.155e36", "1e12", "1e-5", "1e-4", "3.552e45"],
        ["Tapestry Starbirth", "1.989e35", "3.086e17", "1e-4", "1e-3", "1e53"],
        ["Westerlund 2", "1.989e35", "3.086e17", "1e-4", "1e-3", "1e53"],
        ["Pillars of Creation", "1.989e32", "9.46e15", "1e-4", "1e-3", "3.552e48"],
        ["Rings of Relativity", "1.989e36", "3.086e17", "1e-5", "1e-4", "1e54"],
        ["Student's Guide Universe", "1e53", "1e26", "1e-10", "1e-9", "1e80"],
    ], widths=[TW*0.26, TW*0.14, TW*0.12, TW*0.12, TW*0.13, TW*0.13] + [TW*0.10],
       caption="Table 2. MUGESystem parameters for seven astrophysical systems")

    s.append(p("4.  Validation", H1S))
    s.append(p(B("Unit test:") + " " + TT("test_compute_compressed_MUGE(SGR1745-2900)")))
    s.append(p("Expected: " + B("1.782 &#215; 10" + SUP("39") + " m/s" + SUP("2"))))
    s.append(p("Dominated by compressed_base &#215; expansion; B/B" + SUB("crit")
               + " = 0.1 &#8594; 90% retention under linear Meissner form."))

    s.append(p("5.  Physical Interpretation", H1S))
    s.append(p("The [1 &#8722; B/B" + SUB("crit") + "] factor models the Meissner effect: "
               "as the magnetic field B approaches the critical field B" + SUB("crit")
               + ", gravitational coupling is quenched &#8212; mirroring how a Type-I "
               "superconductor expels magnetic flux below B" + SUB("crit")
               + ". The Compressed UQFF thus unifies gravitomagnetic quenching with "
               "cosmological expansion and quantum corrections in a single framework."))
    s.append(p("For Type-II exponential treatment (London penetration depth model) see PAPER_375."))

    s.append(p("6.  Implementation", H1S))
    s.append(p(B("C++:") + " " + TT("STAR_MAGIC_09SEPT_UQFF_MODULE.cpp")
               + ", namespace " + TT("CompressedUQFF")))
    s.append(p(B("Python:") + " " + TT("CondensedPhysics4.py")
               + ", class " + TT("CompressedUQFFBcritSuperconductivityCalculator") + " (CP4 #20)"))
    s.append(p(B("Wolfram macro:") + " " + TT("WOLFRAM_TERM_COMPRESSED_UQFF_BCRIT")))

    s += footer("PAPER_372", "101")
    build(out, s)


# ═══════════════════════════════════════════════════════
# PAPER_373 — Morris-Thorne Wormhole Null Geodesics
# ═══════════════════════════════════════════════════════
def build_373():
    out = os.path.join(DIR, "PAPER_373_MorrisThorne_Wormhole_Null_Geodesics.pdf")
    s = title_block("PAPER_373",
        "Morris-Thorne Wormhole Null Geodesics",
        "101", "grok_share_11254865.txt (lines 2700-2800)")

    s += [p("Significance", ALS),
          p("FIRST wormhole physics integration in the Star Magic / UQFF CP pipeline (CP1&#8211;CP4). "
            "The wormhole coupling to MUGE resonance gravity via the a" + SUB("worm")
            + " term is derived in PAPER_375.", ABS),
          p("Abstract", ALS),
          p("This paper presents the implementation of Morris-Thorne traversable wormhole null "
            "geodesics within the Star Magic UQFF framework. The wormhole metric, geodesic equations, "
            "traversal and reflection conditions, and embedding functions are fully specified. "
            "This constitutes the first wormhole physics to appear in the Condensed Physics "
            "calculator pipeline.", ABS),
          SP(6), rule()]

    s.append(p("1.  Wormhole Metric (Morris-Thorne 1988)", H1S))
    s.append(p("The static, spherically symmetric traversable wormhole metric:"))
    s.append(eq("ds^2 = -dt^2 + dr^2 + (b^2 + r^2)(d_theta^2 + sin^2(theta) d_phi^2)"))
    s.append(p("where b = 1.0 m is the throat radius. This is the simplest traversable wormhole "
               "shape function b(r) = b" + SUB("0") + " = const, corresponding to zero radial "
               "tidal forces at the throat."))

    s.append(p("2.  Null Geodesic Equations", H1S))
    s.append(p("For null geodesics (ds" + SUP("2") + " = 0) with conserved energy "
               "E = dt/d&#955; and angular momentum L = (b" + SUP("2")
               + "+ r" + SUP("2") + ") d&#966;/d&#955;:"))
    s.append(eq("dr/d_lambda   = +/- sqrt( E^2 - L^2/(b^2 + r^2) )"))
    s.append(eq("dphi/d_lambda = L / (b^2 + r^2)"))
    s.append(eq("dt/d_lambda   = E"))
    s.append(p("where &#955; is the affine parameter."))

    s.append(p("3.  Traversal and Reflection Cases", H1S))
    s.append(p(B("Case 1: Traversal") + " (L = 0.5, E = 1.0)"))
    s.append(p("V" + SUB("eff") + "(r) = L" + SUP("2")
               + "/(b" + SUP("2") + "+r" + SUP("2") + ") = 0.25/(1+r" + SUP("2")
               + "). Since V" + SUB("eff") + " < E" + SUP("2") + " = 1 for all r, "
               "the null ray crosses the throat r = 0 and continues to the second universe."))
    s.append(p(B("Case 2: Reflection") + " (L = 1.5, E = 1.0)"))
    s.append(p("V" + SUB("eff") + "(0) = L" + SUP("2") + "/b" + SUP("2")
               + " = 2.25 > E" + SUP("2") + " = 1. The turning radius:"))
    s.append(eq("r_min = sqrt( L^2/E^2 - b^2 ) = sqrt(2.25 - 1) = sqrt(1.25) = 1.118 m"))
    s += tbl([
        ["Case", "L", "E", "Outcome", "r_min (m)"],
        ["Traversal", "0.5", "1.0", "Passes through throat", "0"],
        ["Reflection", "1.5", "1.0", "Reflects at r_min", "~1.118"],
    ], widths=[TW*0.18, TW*0.10, TW*0.10, TW*0.42, TW*0.20],
       caption="Table 1. Null geodesic traversal and reflection cases")

    s.append(p("4.  Embedding Functions", H1S))
    s.append(p("The wormhole embedding in Euclidean 3-space (equatorial slice, &#952; = &#960;/2):"))
    s.append(eq("z_embed(r)   = b * arcsinh(r/b)"))
    s.append(eq("rho_embed(r) = sqrt(b^2 + r^2)"))
    s.append(p("These define the \"funnel\" shape visible in standard wormhole visualisations."))

    s.append(p("5.  Numerical Integration Scheme", H1S))
    s.append(p("The C++ propagator uses first-order Euler (upgradeable to RK4):"))
    s.append(eq("for step = 0..n_steps:"))
    s.append(eq("    dr   = (+/-sqrt(E^2 - L^2/(b^2+r^2))) * d_lambda"))
    s.append(eq("    dphi = (L/(b^2+r^2)) * d_lambda"))
    s.append(eq("    dt   = E * d_lambda"))
    s.append(eq("    r  += dr,  phi += dphi,  t += dt"))
    s.append(p("Defaults: d&#955; = 0.1 m, n_steps = 100."))

    s.append(p("6.  Connection to UQFF", H1S))
    s.append(p("The wormhole throat geometry enters the MUGE resonance model via the "
               + B("wormhole-MUGE coupling term") + " (PAPER_375):"))
    s.append(eq("a_worm(r) = f_worm * E_vac,neb / (b^2 + r^2)"))
    s.append(p("where f" + SUB("worm") + " = 10" + SUP("&#8722;10") + " (wormhole coupling constant), "
               "E" + SUB("vac,neb") + " = 7.09&#215;10" + SUP("&#8722;36") + " J, b = 1.0 m, "
               "r = evaluation radius. At r = b = 1 m: "
               "a" + SUB("worm") + " = 3.545&#215;10" + SUP("&#8722;46") + " m/s" + SUP("2") + "."))

    s.append(p("7.  Implementation", H1S))
    s.append(p(B("C++:") + " " + TT("STAR_MAGIC_09SEPT_UQFF_MODULE.cpp")
               + ", namespace " + TT("WormholeGeodesics")))
    for item in ["Struct: GeodesicState {r, phi, t_coord}",
                 "Functions: drdt(), dphidt(), z_embed(), rho_embed(), propagate(), selftest()"]:
        s.append(p("&#8226;  " + item, BUS))
    s.append(p(B("Python:") + " " + TT("CondensedPhysics4.py")
               + ", class " + TT("MorrisThorneWormholeNullGeodesicsCalculator") + " (CP4 #21)"))
    s.append(p(B("Wolfram macro:") + " " + TT("WOLFRAM_TERM_WORMHOLE_GEODESIC")))

    s.append(p("References", H1S))
    s.append(p("[1] Morris, M.S. &amp; Thorne, K.S. (1988). Wormholes in spacetime and their use for "
               "interstellar travel. " + I("Am. J. Phys.") + " 56(5): 395-412.", RES))
    s.append(p("[2] Visser, M. (1995). " + I("Lorentzian Wormholes.") + " Springer-Verlag.", RES))

    s += footer("PAPER_373", "101")
    build(out, s)


# ═══════════════════════════════════════════════════════
# PAPER_374 — J1610+1811 Relativistic Quasar Jet UQFF-NS Coupling
# ═══════════════════════════════════════════════════════
def build_374():
    out = os.path.join(DIR, "PAPER_374_J1610_Relativistic_Quasar_Jet.pdf")
    s = title_block("PAPER_374",
        "J1610+1811 Relativistic Quasar Jet UQFF&#8211;NS Coupling",
        "101", "grok_share_11254865.txt (lines 5100-5200)")

    s += [p("Abstract", ALS),
          p("This paper presents the coupling of UQFF resonance gravity (from the 12-term MUGE "
            "model, PAPER_371) into a Navier-Stokes quasar jet simulation constrained by the "
            "relativistic parameters of the high-redshift quasar J1610+1811. While PAPER_360 "
            "computed FU and Bi properties for J1610+1811 at z=6.5, this paper addresses the "
            "same object at z=3.122 with a physically motivated relativistic jet velocity "
            "v_SCm = 0.99c derived from observed jet power and luminosity. This represents the "
            "first NS-UQFF coupling simulation driven by actual high-redshift quasar observational "
            "data.", ABS),
          SP(6), rule()]

    s.append(p("1.  Observational Data  (J1610+1811)", H1S))
    s += tbl([
        ["Parameter", "Symbol", "Value", "Source"],
        ["Redshift", "z", "3.122", "Optical / spectroscopic"],
        ["Jet power", "P_jet", "~4e45 W", "X-ray observation"],
        ["Total luminosity", "L", "~2e46 W", "Multi-band"],
        ["Derived jet velocity", "v_SCm", "0.99c = 2.97e8 m/s", "P_jet/L constraint"],
        ["Lorentz factor", "gamma", "~7.09", "Special relativity"],
    ], widths=[TW*0.28, TW*0.14, TW*0.28, TW*0.30],
       caption="Table 1. J1610+1811 observational parameters at z=3.122")
    s.append(p(B("Derivation of v_SCm:") + " The jet velocity is constrained such that the "
               "relativistic kinetic-power ratio P" + SUB("jet") + "/L ≈ 0.2 is consistent "
               "with v" + SUB("SCm") + "/c ≈ 0.99 for a relativistic jet."))
    s.append(p(B("Note on z:") + " PAPER_360 used z=6.5 (the FU/Bi high-z quasar paper); "
               "PAPER_374 uses z=3.122 which appears in the standalone C++ annotated section "
               "of grok_share_11254865.txt (lines 5100+). These represent two distinct epochs."))

    s.append(p("2.  UQFF&#8211;NS Coupling Algorithm", H1S))
    s.append(p("The coupling algorithm proceeds in six steps:"))
    for step in [
        "Instantiate SGR1745-2900 proxy SMBH &#8594; MUGESystem sagA (SgrA* as quasar host).",
        "Compute g_UQFF = compute_resonance_MUGE(sagA, ResonanceParams{}) "
        "&#8594; 12-term MUGE acceleration.",
        "Instantiate FluidSolver (N=32 grid, visc=0.0001, dt=0.1) using Jos Stam Stable Fluids method (PAPER_369).",
        "Set jet_force = v_SCm_rel / 10.0 = 2.97e7 m/s" + SUP("2") + ".",
        "For step=1..10: inject jet force into central column; add UQFF body force "
        "g_UQFF/1e30 uniformly; execute NS step (diffuse &#8594; advect &#8594; project).",
        "Compute mean |v| = (1/N" + SUP("2") + ") &#8721; sqrt(u" + SUB("ij")
        + SUP("2") + "+v" + SUB("ij") + SUP("2") + "); render ASCII velocity field.",
    ]:
        s.append(p(step, BUS))

    s.append(eq("ASCII field:  '#' = |v|>1.0   '+' = |v|>0.5   '.' = |v|>0.1   ' ' = |v|<=0.1"))

    s.append(p("3.  Physical Interpretation", H1S))
    s.append(p("The coupling of g" + SUB("UQFF") + "/1e30 as a body force in the NS grid models "
               "the effect of the vacuum-energy-mediated gravitational field from a SMBH "
               "(SgrA* proxy, M=8.155&#215;10" + SUP("36") + " kg) on the fluid dynamics of a "
               "relativistic jet. The factor 1e30 normalises the UQFF acceleration to the NS grid scale."))

    s.append(p("4.  Distinction from PAPER_360 and PAPER_369", H1S))
    s += tbl([
        ["Feature", "PAPER_360", "PAPER_369", "PAPER_374"],
        ["Object", "J1610+1811", "Generic SgrA*", "J1610+1811"],
        ["Redshift z", "6.5", "--", "3.122"],
        ["Physics", "FU, Bi calculations", "NS Stable Fluids", "NS + UQFF resonance force"],
        ["v_jet", "--", "1e8 m/s (generic)", "0.99c = 2.97e8 m/s"],
        ["UQFF coupling", "None", "Velocity only", "Full 12-term MUGE body force"],
    ], widths=[TW*0.22, TW*0.22, TW*0.22, TW*0.34],
       caption="Table 2. Comparison with related papers")

    s.append(p("5.  Implementation", H1S))
    s.append(p(B("C++:") + " " + TT("STAR_MAGIC_09SEPT_UQFF_MODULE.cpp")
               + ", namespace " + TT("J1610QuasarJet")))
    s.append(p("Constants: z_redshift=3.122, P_jet=4e45, L_luminosity=2e46, v_SCm_rel=0.99c"))
    s.append(p("Function: " + TT("simulate_relativistic_quasar_jet(os, NS_steps=10)")))
    s.append(p(B("Python:") + " " + TT("CondensedPhysics4.py")
               + ", class " + TT("J1610RelativisticQuasarJetUQFFNSCalculator") + " (CP4 #22)"))

    s += footer("PAPER_374", "101")
    build(out, s)


# ═══════════════════════════════════════════════════════
# PAPER_375 — UQFF Advanced Integration
# ═══════════════════════════════════════════════════════
def build_375():
    out = os.path.join(DIR, "PAPER_375_UQFF_Advanced_Integration.pdf")
    s = title_block("PAPER_375",
        "UQFF Advanced Integration:<br/>Wormhole-MUGE | Meissner Exponential"
        "<br/>Relativistic &#947; | Error Propagation",
        "101", "grok_share_11254865.txt (lines 7500-8800)")

    s += [p("Abstract", ALS),
          p("This paper presents four new mathematical formulations extending the UQFF framework: "
            "(1) a wormhole-MUGE coupling term from the Morris-Thorne metric; (2) the exponential "
            "Meissner superconductivity model replacing the linear B/B" + SUB("crit")
            + " quenching (PAPER_372); (3) a special-relativistic Lorentz correction applied to "
            "the DPM acceleration for high-velocity systems; and (4) a formal error propagation "
            "formalism for uncertainty quantification. All four are combined into the Unified UQFF "
            "master equation and validated for J1610+1811 at v=0.99c.", ABS),
          SP(6), rule()]

    s.append(p("1.  Wormhole-MUGE Coupling Term", H1S))
    s.append(p("The Morris-Thorne wormhole geometry (PAPER_373) introduces a new acceleration term:"))
    s.append(eq("a_worm(r) = f_worm * E_vac,neb / (b^2 + r^2)"))
    s.append(p("where f" + SUB("worm") + " = 10" + SUP("&#8722;10") + " (dimensionless coupling), "
               "E" + SUB("vac,neb") + " = 7.09&#215;10" + SUP("&#8722;36") + " J, b = 1.0 m (throat). "
               "This encodes the gravitational contribution of a wormhole throat modulated by local "
               "vacuum energy density."))
    s.append(p("At r = b = 1 m:"))
    s.append(eq("a_worm(1) = 1e-10 * 7.09e-36 / (1 + 1) = 3.545e-46 m/s^2"))

    s.append(p("2.  Meissner Exponential Superconductivity", H1S))
    s.append(p("PAPER_372 uses the linear Meissner approximation (1 &#8722; B/B" + SUB("crit")
               + "). This paper introduces the physically more accurate "
               + B("exponential form") + " applicable to Type-II superconductors "
               "(London penetration depth model):"))
    s.append(eq("Linear:       (1 - B/B_crit)"))
    s.append(eq("Exponential:  exp(-B/B_crit)"))
    s.append(p("The exponential form ensures monotone decay from e" + SUP("0") + " = 1.0 "
               "(no field) to 0 (field well above B" + SUB("crit")
               + "), without the unphysical negative values that arise from the linear form when "
               "B > B" + SUB("crit") + "."))
    s += tbl([
        ["System", "B/B_crit", "Linear factor", "Exponential factor"],
        ["SGR1745-2900", "0.1", "0.9", "0.905"],
        ["Sagittarius A*", "0.1", "0.9", "0.905"],
        ["Student's Guide", "0.1", "0.9", "0.905"],
    ], widths=[TW*0.35, TW*0.20, TW*0.23, TW*0.22],
       caption="Table 1. Linear vs exponential Meissner factor comparison")

    s.append(p("3.  Relativistic Lorentz Correction", H1S))
    s.append(p("For high-velocity systems (J1610+1811 jet at v = 0.99c), the DPM acceleration "
               "undergoes Lorentz suppression:"))
    s.append(eq("gamma = 1 / sqrt(1 - v^2/c^2)"))
    s.append(eq("a_DPM  -->  a_DPM / gamma"))
    s.append(p("For v = 0.99c:"))
    s.append(eq("gamma = 1 / sqrt(1 - 0.9801) = 1 / sqrt(0.0199) = 7.089"))
    s.append(p("This suppresses a" + SUB("DPM") + " by factor ~7.09, reflecting that the DPM "
               "force (electromagnetic in origin) is frame-dependent in the relativistic regime. "
               "All other resonance terms retain their coordinate-frame values."))

    s.append(p("4.  Error Propagation Formalism", H1S))
    s.append(p("Uncertainties in individual MUGE terms propagate in quadrature:"))
    s.append(eq("delta_g = sqrt( sum_i (delta_a_i)^2 )"))
    s.append(p("For a fractional error f = 1%:"))
    s.append(eq("delta_a_i = f * |a_i|  -->  delta_g = f * sqrt( sum_i a_i^2 )"))
    s.append(p("This provides a rigorous uncertainty bound enabling comparison with "
               "observational error bars."))

    s.append(p("5.  Unified UQFF Master Equation  (Complete Form)", H1S))
    s.append(p("Combining PAPER_371&#8211;375:"))
    s.append(eq("[Compressed UQFF — PAPER_372, Meissner exp]:"))
    s.append(eq("  GM/r^2 * (1+H0*t) * exp(-B/B_crit) * (1+F_env)"))
    s.append(eq("+ Ug_sum + Lambda*c^2/3 + quantum_integral + fluid + perturbation"))
    s.append(eq("[MUGE Resonance — PAPER_371, Lorentz corrected]:"))
    s.append(eq("+ a_DPM/gamma + a_THz + a_vac,diff + a_super,freq + a_aether,res + Ug4i"))
    s.append(eq("+ a_quantum,freq + a_Aether,freq + a_fluid,freq + a_osc + a_exp,freq + f_TRZ"))
    s.append(eq("[Wormhole — PAPER_373]:"))
    s.append(eq("+ a_worm = f_worm * E_vac,neb / (b^2 + r^2)"))
    s.append(eq("[Total uncertainty — PAPER_375]:"))
    s.append(eq("+/- delta_g = sqrt( sum_i (delta_a_i)^2 )"))

    s.append(p("6.  Canonical Parameter Summary", H1S))
    s += tbl([
        ["Symbol", "Value", "Source Paper"],
        ["f_worm", "1e-10", "PAPER_375"],
        ["Meissner form", "exp(-B/B_crit)", "PAPER_375 (vs linear PAPER_372)"],
        ["gamma (v=0.99c)", "~7.089", "PAPER_375"],
        ["delta_g (1% error, SGR1745)", "~1e-11", "PAPER_375"],
    ], widths=[TW*0.25, TW*0.35, TW*0.40],
       caption="Table 2. New parameters introduced in PAPER_375")

    s.append(p("7.  Implementation", H1S))
    s.append(p(B("C++:") + " " + TT("STAR_MAGIC_09SEPT_UQFF_MODULE.cpp")
               + ", namespace " + TT("UQFFAdvanced")))
    for fn in ["compute_a_wormhole(Evac_neb, b, r)",
               "meissner_exp(B, Bcrit)",
               "lorentz_gamma(v)",
               "apply_lorentz(aDPM, v)",
               "error_propagation(delta_terms)",
               "compute_unified_UQFF(sys, res, t, v_jet, b_worm, r_worm)",
               "compute_total_uncertainty(sys, p, frac_error)"]:
        s.append(p("&#8226;  " + TT(fn), BUS))
    s.append(p(B("Python:") + " " + TT("CondensedPhysics4.py")
               + ", class " + TT("UQFFWormholeMeissnerRelativisticGammaCalculator") + " (CP4 #23)"))
    s.append(p(B("Wolfram macro:") + " " + TT("WOLFRAM_TERM_UQFF_ADVANCED")))

    s += footer("PAPER_375", "101")
    build(out, s)


# ═══════════════════════════════════════════════════════
# PAPER_376 — UQFF Resonance Superconductive Formal Proof Set
# ═══════════════════════════════════════════════════════
def build_376():
    out = os.path.join(DIR, "PAPER_376_UQFF_Resonance_Formal_Proof_Set.pdf")
    s = title_block("PAPER_376",
        "UQFF Resonance Superconductive Formal Proof Set",
        "102", "grok_share_11254865.txt (lines 6001-10322)")

    s += [p("Overview", ALS),
          p("This paper formalizes the mathematical proof set validating the UQFF Resonance "
            "Superconductive model. The proof set document (May 15, 2025) provides dimensional "
            "consistency checks, boundary condition tests, resonance amplification proofs, "
            "superconductivity proofs, and empirical validation against astrophysical observations. "
            "CP4 class: " + TT("UQFFResonanceFormalProofSetCalculator") + " (CP4 #25).", ABS),
          SP(6), rule()]

    s.append(p("1.  Overview", H1S))
    s.append(p("Source documents: "
               + I("UQFF_Resonance Superconductive Universal Gravity Equation system proof set._15May2025.docx") + "; "
               + I("Compressed UQFF Equation_14May2025.docx") + "; "
               + I("Master UQFF Resonance Equation_14May2025.docx") + "."))
    s.append(p("Five proofs are presented: (1) dimensional consistency, (2) boundary conditions, "
               "(3) resonance amplification, (4) superconductivity, (5) empirical validation."))

    s.append(p("2.  Proof 1: Dimensional Consistency", H1S))
    s.append(p("All MUGE terms yield units of m/s" + SUP("2") + " (acceleration)."))
    s += tbl([
        ["Term", "Dimensional Form", "Unit Check"],
        ["Base", "G*M/r^2", "m3/(kg*s2) * kg / m2 = m/s2  OK"],
        ["Expansion", "(1 + H0*t) [dimensionless]", "multiplier  OK"],
        ["Super_adj", "(1 - B/B_crit) [dimensionless]", "multiplier  OK"],
        ["Cosm", "Lambda*c^2/3", "m-2 * (m/s)2 = s-2/m [contextual scaling]"],
        ["Quantum", "(hbar/dxdp)*integral*2pi/t_H", "scaled to m/s2 via normalisation"],
        ["Fluid", "rho_fluid*V*g_local", "kg/m3 * m3 * m/s2 [scaled]"],
        ["Perturbation", "(M+M_DM)*(drho/rho+3GM/r^3)", "m/s2 via mass normalisation"],
    ], widths=[TW*0.18, TW*0.42, TW*0.40],
       caption="Table 1. Dimensional consistency check for Compressed MUGE terms")
    s.append(p("All 12 Resonance MUGE terms (a_DPM through f_TRZ) likewise reduce to m/s" + SUP("2")
               + " through the E_vac,neb / E_vac,ISM / c normalisation chain."))

    s.append(p("3.  Proof 2: Boundary Conditions", H1S))
    s.append(eq("As r --> inf:   g_UQFF --> Lambda*c^2/3 ~ 3.3e-36 m/s^2  (cosmological floor)"))
    s.append(eq("As t --> 0:     g_UQFF --> G*M/r^2  (Newtonian gravity recovered)"))
    s.append(eq("As B --> B_crit (linear):  g_UQFF * (1 - B/B_crit) --> 0  (full quench)"))
    s.append(eq("As B --> B_crit (exp):     g_UQFF * exp(-1) ~ 0.368*g  (gradual suppression)"))
    s.append(eq("As B >> B_crit:            Meissner complete expulsion, g --> 0"))

    s.append(p("4.  Proof 3: Resonance Amplification", H1S))
    s.append(p("The quantum coherence integral amplifies at the cosmic resonance frequency:"))
    s.append(eq("omega_res = 2*pi / t_Hubble = 2*pi / 4.35e17 s ~ 1.445e-17 rad/s"))
    s.append(p("Key value: f" + SUB("quantum") + " = 1.445&#215;10" + SUP("&#8722;17")
               + " in ResonanceParams matches this frequency exactly. The a" + SUB("quantum,freq")
               + " term therefore ensures the Hubble time resonance is present in every MUGE computation."))
    s.append(eq("a_quantum,freq = f_quantum * E_vac,neb * a_DPM / E_vac,ISM / c"))
    s.append(eq("               = 1.445e-17 * 7.09e-36 * a_DPM / 7.09e-37 / 3e8"))

    s.append(p("5.  Proof 4: Superconductivity Proof", H1S))
    s.append(p(B("Linear Meissner (PAPER_372):") + "  g_adj = g_base &#215; (1 &#8722; B/B_crit)"))
    s.append(p(B("Exponential Meissner (PAPER_375/376):") + "  g_adj = g_base &#215; exp(&#8722;B/B_crit)"))
    s.append(p("Physical basis: London penetration depth &#955;" + SUB("L")
               + " &#8733; 1/&#8730;n" + SUB("s")
               + ", where n" + SUB("s") + " is superfluid carrier density. "
               "The exponential form applies for Type-II superconductors (Abrikosov vortex lattice)."))
    s += tbl([
        ["Condition", "Linear factor", "Exponential factor", "Physical meaning"],
        ["B = 0", "1.0", "1.0", "No superconductivity"],
        ["B = 0.1*B_crit", "0.9", "0.905", "Partial Meissner"],
        ["B = B_crit", "0.0 (exact quench)", "0.368 = e-1", "Exponential preferred"],
        ["B >> B_crit", "< 0 (unphysical)", "~0 (physical)", "Exponential necessary"],
    ], widths=[TW*0.22, TW*0.22, TW*0.22, TW*0.34],
       caption="Table 2. Meissner linear vs exponential comparison at key field values")

    s.append(p("6.  Proof 5: Empirical Validation", H1S))
    s.append(p("6.1  " + B("Magnetar SGR 1745-2900"), H2S))
    s.append(p(B("Observed:") + " X-ray flare timescales 10&#8211;100 days (Chandra, NuSTAR)"))
    s.append(p(B("UQFF Prediction:") + " E" + SUB("react")
               + "(t) = 10" + SUP("46") + " &#215; exp(&#8722;0.0005 &#215; t)"))
    s.append(eq("At t=10 days:   E_react ~ 1046 * exp(-5e-3)  ~ 1041 J  (flare active)"))
    s.append(eq("At t=100 days:  E_react ~ 1046 * exp(-0.05)  ~  995 J  (flare active)"))
    s.append(p("Flare active while E" + SUB("react") + " > threshold: &#8776; 10&#8211;100 day window &#10003;"))
    s.append(p("6.2  " + B("Sagittarius A* (Sgr A*)"), H2S))
    s.append(p(B("Observed:") + " Accretion rate &#8776; 10" + SUP("&#8722;8")
               + " M" + SUB("&#9737;") + "/yr (Event Horizon Telescope)"))
    s.append(p(B("UQFF Prediction:") + " resonance_MUGE(Sgr A*) &#8776; 4.105&#215;10" + SUP("29")
               + " m/s" + SUP("2") + " &#8212; consistent with high-luminosity EHT flares 2022&#8211;2025."))
    s.append(p("6.3  " + B("Newtonian Baseline"), H2S))
    s.append(eq("test_compute_compressed_base() at 1 AU:"))
    s.append(eq("  Expected: G*M_sun / (1 AU)^2 ~ 0.00593 m/s^2   Computed: OK (passes)"))

    s.append(p("7.  Unified Proof Equation  (Combined Form)", H1S))
    s.append(eq("[Compressed UQFF]:"))
    s.append(eq("  GM(t)/r^2*(1+H0*t)*exp(-B(t)/Bcrit)*(1+Fenv(t))"))
    s.append(eq("  + Ug_sum + Lambda*c^2/3 + quantum + fluid + perturbation"))
    s.append(eq("[MUGE Resonance]:"))
    s.append(eq("  + a_DPM/gamma + a_THz + a_vac,diff + a_super,freq + a_aether,res + Ug4i"))
    s.append(eq("  + a_quantum,freq + a_Aether,freq + a_fluid,freq + a_osc + a_exp,freq + f_TRZ"))
    s.append(eq("[Wormhole + Error]:"))
    s.append(eq("  + a_worm +/- delta_g = sqrt(sum_i (delta_a_i)^2)"))
    s.append(p("where &#947; = 1/&#8730;(1&#8722;v" + SUP("2")
               + "/c" + SUP("2") + ") (Lorentz factor, &#947; &#8776; 7.09 at v=0.99c)."))

    s.append(p("8.  Key Validated Constants", H1S))
    s += tbl([
        ["Parameter", "Value", "Proof Context"],
        ["H0", "2.269e-18 s-1", "Expansion baseline (Planck 2018 consistent)"],
        ["Lambda", "1.1e-52 m-2", "Cosmological constant (LCDM)"],
        ["hbar", "1.0546e-34 J*s", "Quantum coherence integral"],
        ["t_Hubble", "4.35e17 s", "Resonance amplification timescale"],
        ["B_crit", "1e11 T", "Magnetar critical field (PAPER_372)"],
        ["f_quantum", "1.445e-17 Hz", "= 2pi/t_Hubble (Hubble resonance)"],
        ["E_react(t=0)", "1e46 J", "Magnetar flare energy seed"],
        ["kappa", "0.0005 day-1", "SCm reactivity decay; matches 10-100 day flares"],
    ], widths=[TW*0.18, TW*0.22, TW*0.60],
       caption="Table 3. Key validated constants in the UQFF formal proof set")

    s.append(p("9.  CP4 Class", H1S))
    s.append(p(B("Class:") + " " + TT("UQFFResonanceFormalProofSetCalculator")))
    s.append(p(B("Category:") + " Formal Validation"))
    s.append(p(B("References:") + " PAPER_372, PAPER_373, PAPER_374, PAPER_375"))

    s += footer("PAPER_376", "102")
    build(out, s)


# ═══════════════════════════════════════════════════════
if __name__ == "__main__":
    build_371()
    build_372()
    build_373()
    build_374()
    build_375()
    build_376()
    print("\nAll 5 PDFs complete.")
