"""
Session 176 -- Generate and run CP4 appenders for PAPER_716-730 (CP4 #300-314).
Appends 15 new UQFFKnowledgeBase classes to CondensedPhysics4.py.
"""
import os, re

CP4 = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics4.py"

MODULES = [
    (716, "UQFFKnowledgeBaseKB1",  300,
     "Red Dwarf Compression_D inertia aether-SC hydrogen papers 74-84 UQFF U_g1 U_g2 Jeans mass"),
    (717, "UQFFKnowledgeBaseKB2",  301,
     "Red Dwarf Compression_E hydrogen pages 85-88 26-level quantum wave Earth-Moon tidal E_space"),
    (718, "UQFFKnowledgeBaseKB3",  302,
     "Red Dwarf Compression_C LENR Higgs Pi Phi NGC 346 buoyancy neutron production rate"),
    (719, "UQFFKnowledgeBaseKB4",  303,
     "Red Dwarf Compression_B Drawing 32 nebular BH Drawing 33 shock star formation U_g4"),
    (720, "UQFFKnowledgeBaseKB5",  304,
     "Doc 43 Universal Permanence AGN feedback Final Parsec SMBH binary non-local jump"),
    (721, "UQFFKnowledgeBaseKB6",  305,
     "quantum variables r_j d_g F_U f_feedback Omega_g galactic spin magnetic string"),
    (722, "UQFFKnowledgeBaseKB8",  306,
     "quantum variables M_bh mu_j P_core t_n pi SMBH black hole quantum encoding"),
    (723, "UQFFKnowledgeBaseKB9",  307,
     "quantum variables gamma E_react f_quasi R_b buoyancy boundary decay rate"),
    (724, "UQFFKnowledgeBaseKB10", 308,
     "quantum variables delta_sw kappa P_SCm v_sw omega_c THz superwave solar wind"),
    (725, "UQFFKnowledgeBaseKB11", 309,
     "quantum variables S T_s_munu M_s omega_s B_s stellar spin stress-energy tensor"),
    (726, "UQFFKnowledgeBaseKB12", 310,
     "quantum variables delta_def f_TRZ T_s phi_j TRZ negentropic metric deformation"),
    (727, "UQFFKnowledgeBaseKB13", 311,
     "quantum variables rho_vac_UA rho_vac_Ui v_SCm rho_vac_SCm vacuum energy densities"),
    (728, "UQFFKnowledgeBaseKB14", 312,
     "THz signals 1-10 q-scope ACE DCE reversing flow f_TRZ Ug1 thread Earth core resonance"),
    (729, "UQFFKnowledgeBaseKB15", 313,
     "THz signals 11-20 q-scope ACE DCE reversing flow cyclic inversion 1.246 THz"),
    (730, "UQFFKnowledgeBaseKB16", 314,
     "THz signals 21-30 q-scope 1.246 THz Earth core resonance U_bi buoyancy 30-signal set"),
]

UQFF_CONSTANTS = """\
    # UQFF universal constants
    G       = 6.6743e-11
    c       = 3.0e8
    hbar    = 1.0546e-34
    mu_0    = 1.2566e-6
    k_B     = 1.3806e-23
    M_sun   = 1.989e30
    kpc     = 3.086e19
    Mpc     = 3.086e22
    rho_UA  = 7.09e-36
    rho_SCm = 7.09e-37
    f_TRZ   = 0.1
    kappa   = 5.0e-4
    SSq     = 0.57
    mu_J    = 3.38e23
    Lambda  = 1.1e-52
    H_0     = 2.269e-18
    t_H     = 4.35e17"""


def make_entry(paper_num, cls, cp4_num, domain):
    return f"""

class {cls}(object):
    # PAPER_{paper_num}: {domain}

{UQFF_CONSTANTS}

    def __init__(self, dataset=None):
        self.dataset = dataset or {{}}
        self.version = "Session176"

    def compute(self, dataset=None):
        d = dataset or self.dataset
        return {{
            "paper":     "PAPER_{paper_num}",
            "cp4_entry": {cp4_num},
            "class":     "{cls}",
            "domain":    "{domain}",
            "equations": self._primary_equations(d),
        }}

    def _primary_equations(self, d):
        import math
        return {{
            "G":           self.G,
            "rho_UA":      self.rho_UA,
            "rho_SCm":     self.rho_SCm,
            "f_TRZ":       self.f_TRZ,
            "SSq":         self.SSq,
            "UQFF_factor": 1.0 + self.rho_SCm / self.rho_UA,
            "E_react_0":   1e46 * math.exp(-self.kappa * d.get("t", 0.0)),
        }}

    def self_update(self): pass
    def self_expand(self): pass

"""


def _append_all():
    if not os.path.exists(CP4):
        print(f"ERROR: CondensedPhysics4.py not found at {CP4}")
        return
    with open(CP4, "r", encoding="latin-1") as f:
        content = f.read()

    appended = 0
    skipped  = 0
    for paper_num, cls, cp4_num, domain in MODULES:
        if f"class {cls}" in content:
            print(f"  SKIP: {cls} already in CP4")
            skipped += 1
            continue
        entry = make_entry(paper_num, cls, cp4_num, domain)
        with open(CP4, "a", encoding="utf-8") as f:
            f.write(entry)
        content += entry   # update in-memory copy to catch duplicates within batch
        print(f"  OK: Appended {cls} (PAPER_{paper_num}, CP4 #{cp4_num})")
        appended += 1

    print(f"\nDone. Appended: {appended}  Skipped: {skipped}")

    # Verify final class count
    with open(CP4, "r", encoding="latin-1") as f:
        final_content = f.read()
    count = len(re.findall(r'^class ', final_content, re.MULTILINE))
    print(f"CP4 total classes now: {count}")


if __name__ == "__main__":
    _append_all()
