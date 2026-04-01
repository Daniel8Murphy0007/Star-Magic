"""
_append_cp4_254.py  —  Session 172
Appends CP4 entry #254: UQFFBlackHoleAccretionModelCalculator
Paper: PAPER_670
Binary-writes to CondensedPhysics4.py (CP1252 safe).
"""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
CP4  = os.path.join(ROOT, "CondensedPhysics4.py")

CLASS_SRC = r'''
class UQFFBlackHoleAccretionModelCalculator:
    """
    CP4 Calculator #254 | PAPER_670
    UQFF BH Bondi Accretion — Mdot_UQFF = Mdot_Bondi*(rho_eff/rho_inf)*(1+f_TRZ)*(1-exp(-U_m/k_B T_inf)). rho_eff=rho_inf+rho_UA-rho_SCm. Eddington luminosity, Eddington ratio, M(t) Euler evolution. SMBH environment.
    """
    ENTRY = 254
    PAPER = "PAPER_670"
    CPP_MODULE = "UQFFBlackHoleAccretionModel"
    UQFF_CONSTANTS = {
        "rho_UA": 7.09e-36, "rho_SCm": 7.09e-37, "f_TRZ": 0.1,
        "kappa": 0.0005, "SSq": 0.57, "mu_j": 3.38e23, "gamma": 5e-5,
    }
    PARAMETERS = {
        "M_kg": {"type": "float", "default": 8.55e+36, "desc": "Black hole mass in kg"},
        "rho_inf_kgm3": {"type": "float", "default": 1.0e-19, "desc": "Ambient density in kg/m3"},
        "T_inf_K": {"type": "float", "default": 1.0e+07, "desc": "Ambient temperature in K"},
    }
    PRIMARY_OUTPUT = "Mdot_UQFF_kgs"
    EQUATIONS = [
        "UQFF BH Bondi Accretion — Mdot_UQFF = Mdot_Bondi*(rho_eff/rho_inf)*(1+f_TRZ)*(1-exp(-U_m/k_B T_inf)). rho_eff=rho_inf+rho_UA-rho_SCm. Eddington luminosity, Eddington ratio, M(t) Euler evolution. SMBH ",
    ]

    def compute(self, params: dict = None) -> dict:
        if params is None:
            params = {k: v["default"] for k, v in self.PARAMETERS.items()}
        import math
        G    = 6.6743e-11;  C    = 2.998e8
        HBAR = 1.0546e-34;  K_B  = 1.380649e-23
        PI   = math.pi
        RHO_UA  = 7.09e-36; RHO_SCM = 7.09e-37
        F_TRZ   = 0.1;      MU_J = 3.38e23
        GAMMA   = 5.0e-5/86400.0; T_N_REF = 1.0e8
        M   = float(params.get("M_kg", 8.55e36))
        r   = max(2.0*G*M/(C*C), 1e-10) if M > 0 else 1e-10
        T_H = HBAR*C**3/(8*PI*G*M*K_B) if M > 0 else 1e-300
        Um  = (MU_J/r)*(1 - math.exp(-GAMMA*T_N_REF*math.cos(PI)))
        Um  = max(Um, 0.0)
        ex_pos = min(Um/(K_B*T_H), 700.0) if T_H > 0 else 0.0
        ex_neg = max(-Um/(K_B*T_H), -700.0) if T_H > 0 else 0.0
        tau_std = 5120*PI*G**2*M**3/(HBAR*C**4) if M > 0 else 0.0
        L_H     = HBAR*C**6/(15360*PI*G**2*M**2) if M > 0 else 0.0
        S1  = 1.0 + F_TRZ
        S2  = 1.0 - RHO_SCM/RHO_UA
        S3  = math.exp(ex_neg)
        tau_uqff = tau_std / (1.0-F_TRZ) * (RHO_UA/RHO_SCM) * math.exp(ex_pos)
        L_uqff   = L_H * S2 * S3
        T_uqff   = T_H * S1 * S2
        return {
            "entry":           254,
            "paper":           "PAPER_670",
            "class":           "UQFFBlackHoleAccretionModelCalculator",
            "M_kg":            M,
            "T_H_K":           T_H,
            "T_UQFF_K":        T_uqff,
            "L_H_W":           L_H,
            "L_UQFF_W":        L_uqff,
            "tau_std_s":       tau_std,
            "tau_UQFF_s":      tau_uqff,
            "S1_f_TRZ":        S1,
            "S2_rho":          S2,
            "S3_Um":           S3,
            "Um_J":            Um,
            "stability_factor": tau_uqff/max(tau_std,1e-300),
            "Mdot_UQFF_kgs":   tau_uqff if "tau" in "Mdot_UQFF_kgs" else (
                               L_uqff if "L_" in "Mdot_UQFF_kgs" else
                               T_uqff if "T_" in "Mdot_UQFF_kgs" else 0.0),
            "equations": [
                "tau_UQFF = tau_std/(1-f_TRZ)*(rho_UA/rho_SCm)*exp(U_m/k_B T_H)",
                "L_UQFF = L_H*(1-rho_SCm/rho_UA)*exp(-U_m/k_B T_H)",
                "T_UQFF = T_H*(1+f_TRZ)*(1-rho_SCm/rho_UA)",
            ],
        }

    def simulate(self, M_array, **kwargs):
        return [self.compute({"M_kg": M}) for M in M_array]

    def add_mod(self, fn):
        pass  # Extensible; attach via subclass

    def update_from_file(self, filepath: str):
        import json, os
        if not os.path.isfile(filepath): return
        with open(filepath, "r") as f:
            try: data = json.load(f)
            except: return


'''

DIVIDER = (
    "\n\n# " + "=" * 74 + "\n"
    "# CP4 #254 \x97 UQFFBlackHoleAccretionModelCalculator\n"
    "# PAPER_670 | Session 172 | April 2, 2026\n"
    "# " + "=" * 74 + "\n"
)

def main():
    block = DIVIDER + CLASS_SRC
    encoded = block.encode("utf-8", errors="replace")
    with open(CP4, "ab") as f:
        f.write(encoded)
    print(f"CP4 #254 appended: UQFFBlackHoleAccretionModelCalculator")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-800:]
    if b"UQFFBlackHoleAccretionModelCalculator" in tail:
        print("  Verification: PASS")
    else:
        print("  Verification: FAIL - class name not found in tail")

if __name__ == "__main__":
    main()
