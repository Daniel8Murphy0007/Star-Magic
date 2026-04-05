#!/usr/bin/env python3
"""Session 200C — generate 8 whitepapers (.md) + 8 PDFs (reportlab)."""
import os, textwrap

ROOT = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
PDF_DIR = os.path.join(ROOT, "pdf")

PAPERS = [
    {
        "num": 870, "cls": "#454", "name": "DPMExtendedPeriodicTableProportionCalc",
        "title": "DPM Extended Periodic Table Proportion System",
        "abstract": (
            "Develops a system of unique DPM proportions (fUA', fSCm) for every atom "
            "in the universe, extending the Periodic Table to Z_max=10,000. "
            "fUA'=(Z_max-Z)/Z_max; fSCm=Z/Z_max. All atoms start radioactive with "
            "decay rate lambda=k_lambda*f_SCm and stabilize over time. "
            "Reactivity gradient R_EB=k_R*Z. SM_magnetic/non-magnetic alternation."
        ),
        "equations": [
            "f_UA' = (Z_max - Z) / Z_max",
            "f_SCm = Z / Z_max",
            "f_UA' + f_SCm = 1   (DPM fully defines nucleus)",
            "R_EB = k_R * Z   (reactivity gradient)",
            "lambda = k_lambda * f_SCm   (radioactive decay rate)",
            "L_quant proportional to f_UA' * f_SCm * R_EB",
        ],
    },
    {
        "num": 871, "cls": "#455", "name": "UniversalSpeedRangeCosmicPhotonDecelerationCalc",
        "title": "Universal Speed Range c^26*i^-26 and Cosmic Photon Deceleration",
        "abstract": (
            "Defines the universal speed range as c^26*i^-26, where c=2.998e8 m/s. "
            "E=c^2 is redefined as the speed of visible light in the cosmic vacuum (UA). "
            "Cosmic photons originate as heavy metal ions emitted from nuclear shells, "
            "decelerating from c^26*i^-26 to c^2 through 26 quantum layers."
        ),
        "equations": [
            "v_range = c^26 * i^{-26}   (universal speed range)",
            "E = c^2 = (2.998e8)^2 ~ 8.988e16 m^2/s^2   (light speed in UA)",
            "v_photon(layer) = c^{26-layer+1}",
            "Cosmic photon: heavy metal ion -> slows c^26*i^-26 -> c^2",
        ],
    },
    {
        "num": 872, "cls": "#456", "name": "ProtoIronProtoSiliconNuclearIdentityCalc",
        "title": "Proto-Iron / Proto-Silicon Nuclear Identity Mapping",
        "abstract": (
            "Maps proto-atomic nuclei to their heavier-element identity via DPM proportions. "
            "Proto-hydrogen nucleus = proto-iron (Z_identity=26, SM_magnetic, durable "
            "strong-force shell). Proto-helium nucleus = proto-silicon (Z_identity=14, "
            "SM_non-magnetic). Odd Z -> SM_magnetic; even Z -> SM_non-magnetic."
        ),
        "equations": [
            "Proto-H nucleus = Proto-Fe (Z_id=26, SM_magnetic)",
            "Proto-He nucleus = Proto-Si (Z_id=14, SM_non-magnetic)",
            "U_m = f_SCm * rho_SCm * c^2   (SCm-only influence)",
            "SM_property: odd Z -> magnetic, even Z -> non-magnetic",
        ],
    },
    {
        "num": 873, "cls": "#457", "name": "Ug1DPMGeophysicalGeometrySummationCalc",
        "title": "U_g1=DPM Geophysical Geometry Summation",
        "abstract": (
            "U_g1 is equivalent to the DPM; the total force is a summation over DPM "
            "variable forms reflecting unique forces with geophysical geometries. "
            "Components: SM_gravity (spherical), U_b buoyancy (spherical, counter-force), "
            "resonance (toroidal). Geometry functions G_k enable force diagramming."
        ),
        "equations": [
            "F_Ug1 = Sum_k [k_k * (fUA1*fSCm1*REB1)*(fUA2*fSCm2*REB2)/r^2 * G_k]",
            "F_SM_gravity = k_grav * DPM1*DPM2/r^2 * sin(theta)",
            "F_Ub = k_buoy * (fUA'*fSCm)/r^2 * sin(theta)",
            "F_resonance = k_res * nu_THz * fSCm * cos(phi)*f(nu_THz)",
            "F_Ug1_total = F_grav - F_buoy + F_resonance",
        ],
    },
    {
        "num": 874, "cls": "#458", "name": "Ug3ElectronTaggingTHzCirculationCalc",
        "title": "U_g3=U_i+U_m Electron Tagging and THz Circulation",
        "abstract": (
            "U_g3 = U_i + U_m in motion. U_i individually tags each electron via the "
            "THz hole pipeline from nucleus (Point A) to electron (Point B), counteracting "
            "the strong-force to position electrons at projected quantum energy shell balance. "
            "Electron DPM uses limited U_i for coherent circulation in U_g2 imaginary orbit shell. "
            "U_g2 monitors and adjusts shell position."
        ),
        "equations": [
            "F_Ug3 = (F_Ui + F_Um + F_DPM_e) * G_geo / r_shell^2",
            "F_Ui = k_i * fUA'_nuc * nu_THz * R_EB   (repulsive tagging)",
            "F_Um = k_m * fSCm_nuc * nu_res   (SCm-driven motion)",
            "F_DPM_e = k_e * (fUA'_e * fSCm_e) * nu_THz   (electron circulation)",
            "G_geo = sin(theta)*cos(phi)*f(nu_THz)",
            "E_shell = c * nu_res * hbar * fSCm * sin(theta)   (U_g2 monitor)",
        ],
    },
    {
        "num": 875, "cls": "#459", "name": "SMMagSurfaceConductionFragmentAssemblyCalc",
        "title": "SM_mag Surface Conduction and Fragment Assembly",
        "abstract": (
            "SM_magnetic surface moments conduct from many surface points in a chaos "
            "pattern but spatially coherent separation, induced by the internal DPM "
            "through the semi-solid shell. SM_mag arranges brittle layered string fragments "
            "on the durable proto-nucleus. Vacuum energy density becomes a fixed constant "
            "(capacitance) at proto-nucleus formation. ULF quantum ripples ULF_quantum^{-1..-26}."
        ),
        "equations": [
            "C_vac = rho_vac * r   (vacuum energy density = capacitance)",
            "E_crack = Sum_{i=1}^{26} (hbar*omega/i) * C_vac   (ULF quantum ripples)",
            "SM_mag = fSCm * rho_SCm * c^2   (surface conduction strength)",
            "coherence = 1 - e^{-fSCm*N_frag}   (chaos -> coherent transition)",
            "g_atomic_SM proportional to 1/(rho_vac * r^2)   (inversely proportional to vacuum density)",
        ],
    },
    {
        "num": 876, "cls": "#460", "name": "DPMCoherentConsciousnessSpookyActionCalc",
        "title": "DPM Coherent Consciousness and Spooky Action at a Distance",
        "abstract": (
            "The DPM (UA'+SCm) is coherent and produces shared consciousness both locally "
            "and at spooky distances. THz hole synchronization produces paired resonance "
            "that drives systems to collide and/or repel from galactic distances. "
            "Directly responsible for galactic red/blue shifting both locally and at "
            "spooky distance."
        ),
        "equations": [
            "coherence = fSCm_A * fSCm_B",
            "E_pair = rho_SCm * c^2 * nu_THz * coherence",
            "Delta_nu = k_shift * nu_THz * (fSCm_A + fSCm_B)   (red/blue shift)",
            "z = Delta_nu / nu_emit",
            "C_consciousness = coherence * (1 - e^{-nu_THz/10^12})",
            "E_surplus = E_pair / d^2   (signal bundle at distance)",
        ],
    },
    {
        "num": 877, "cls": "#461", "name": "ThreeAssumptionUQFFCosmogenesisCalc",
        "title": "Three-Assumption UQFF Cosmogenesis Master Integration",
        "abstract": (
            "Master integration of the complete UQFF cosmogenesis model: "
            "(1) Three reactive quantum fundamentals (electrostatic barrier, UA, SCm) "
            "form proto-nuclear shells via DPM. "
            "(2) Proto-shells evolve through EM bang and 2 expansion/contraction cycles "
            "to produce proto-atoms (proto-H=proto-Fe, proto-He=proto-Si). "
            "(3) Four U_g forces govern all interactions: U_g1=DPM, U_g2=shells, "
            "U_g3=U_i+U_m, U_g4i=control. "
            "26 quantum atomic states before mass; quantum-to-mass gradient at 7-10 U_mag degrees."
        ),
        "equations": [
            "=== Assumption 1 ===",
            "f_UA' = (Z_max-Z)/Z_max; f_SCm = Z/Z_max; R_EB = Z",
            "=== Assumption 2 (ACP) ===",
            "U_i = k*(rho_SCm - rho_UA/10)*omega*cos(pi*t)",
            "Psi_proto = Sum_{i=1}^{26} U_m,i",
            "C_vac = rho_vac*r   (capacitance = vacuum energy density)",
            "=== Assumption 3 (Forces) ===",
            "F_Ug1 = DPM summation with geometry",
            "E_Ug2 = c*nu*hbar*fSCm",
            "F_Ug3 = (U_i + U_m)/r^2",
            "E_Ug4i = fSCm*nu*rho_SCm",
        ],
    },
]


def create_whitepapers():
    """Create markdown whitepapers."""
    created = []
    for p in PAPERS:
        fname = f"PAPER_{p['num']:03d}_{p['name']}.md"
        fpath = os.path.join(ROOT, fname)
        eqs = "\n".join(f"  {i+1}. `{e}`" for i, e in enumerate(p["equations"]))
        content = textwrap.dedent(f"""\
            # PAPER_{p['num']:03d} — {p['title']}

            **CP4 Class:** {p['cls']} `{p['name']}`
            **Session:** 200C v5.61
            **Source:** `describe mass without using weight.txt` (3,094 lines)

            ## Abstract

            {p['abstract']}

            ## Key Equations

            {eqs}

            ## UQFF Context

            This calculator implements parameterized physics equations per the
            CondensedPhysics architecture (pure calculator, no hardcoded data).
            Inputs arrive from source2.cpp PRINCIPAL GUI via dataset dict.
            Outputs include primary_equations, available_equations, and simulation_set.

            ## Source Thread

            Grok 3 dialogue, Daniel Murphy, June 03-04 2025.
            Copyright Daniel T. Murphy, daniel.murphy00@gmail.com.

            ---
            *Star-Magic UQFF Codebase — Session 200C*
        """)
        with open(fpath, "w", encoding="utf-8") as f:
            f.write(content)
        created.append(fname)
        print(f"  [WP] {fname}")
    return created


def create_pdfs():
    """Create PDF whitepapers using reportlab."""
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import mm
        from reportlab.pdfgen import canvas
    except ImportError:
        print("  reportlab not available — skipping PDF generation")
        return []

    os.makedirs(PDF_DIR, exist_ok=True)
    created = []
    for p in PAPERS:
        fname = f"PAPER_{p['num']:03d}_{p['name']}.pdf"
        fpath = os.path.join(PDF_DIR, fname)
        c = canvas.Canvas(fpath, pagesize=A4)
        w, h = A4
        y = h - 40 * mm

        c.setFont("Helvetica-Bold", 16)
        c.drawString(20 * mm, y, f"PAPER_{p['num']:03d}")
        y -= 8 * mm
        c.setFont("Helvetica-Bold", 12)
        # Wrap title if needed
        title = p["title"]
        if len(title) > 70:
            c.drawString(20 * mm, y, title[:70])
            y -= 6 * mm
            c.drawString(20 * mm, y, title[70:])
        else:
            c.drawString(20 * mm, y, title)
        y -= 10 * mm

        c.setFont("Helvetica", 9)
        c.drawString(20 * mm, y, f"CP4 Class: {p['cls']}  {p['name']}")
        y -= 5 * mm
        c.drawString(20 * mm, y, "Session 200C v5.61 | Source: describe mass without using weight.txt")
        y -= 10 * mm

        c.setFont("Helvetica-Bold", 11)
        c.drawString(20 * mm, y, "Abstract")
        y -= 6 * mm
        c.setFont("Helvetica", 9)
        # Word-wrap abstract
        words = p["abstract"].split()
        line = ""
        for word in words:
            test = line + " " + word if line else word
            if len(test) > 90:
                c.drawString(20 * mm, y, line)
                y -= 4.5 * mm
                line = word
                if y < 30 * mm:
                    c.showPage()
                    y = h - 30 * mm
                    c.setFont("Helvetica", 9)
            else:
                line = test
        if line:
            c.drawString(20 * mm, y, line)
            y -= 6 * mm

        c.setFont("Helvetica-Bold", 11)
        c.drawString(20 * mm, y, "Key Equations")
        y -= 6 * mm
        c.setFont("Courier", 8)
        for eq in p["equations"]:
            if y < 25 * mm:
                c.showPage()
                y = h - 30 * mm
                c.setFont("Courier", 8)
            c.drawString(22 * mm, y, eq[:95])
            y -= 4.5 * mm

        y -= 6 * mm
        c.setFont("Helvetica", 8)
        c.drawString(20 * mm, y, "Star-Magic UQFF Codebase - Session 200C - Copyright Daniel T. Murphy")

        c.save()
        created.append(fname)
        print(f"  [PDF] {fname}")

    return created


def main():
    print("Creating whitepapers...")
    wps = create_whitepapers()
    print(f"\n{len(wps)} whitepapers created.\n")

    print("Creating PDFs...")
    pdfs = create_pdfs()
    print(f"\n{len(pdfs)} PDFs created.")

    # Count total PDFs
    if os.path.isdir(PDF_DIR):
        total = len([f for f in os.listdir(PDF_DIR) if f.endswith(".pdf")])
        print(f"Total PDFs in pdf/: {total}")


if __name__ == "__main__":
    main()
