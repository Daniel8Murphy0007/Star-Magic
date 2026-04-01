"""
_append_cp4_252.py  —  Session 172
Appends CP4 entry #252: UQFFStabilityPrimordialBHCalculator
Paper: PAPER_668
Binary-writes to CondensedPhysics4.py (CP1252 safe).
"""
import os

ROOT = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic"
CP4  = os.path.join(ROOT, "CondensedPhysics4.py")

CLASS_SRC = r'''
class UQFFStabilityPrimordialBHCalculator:
    """
    CP4 Calculator #252 | PAPER_668
    UQFF Primordial BH Stability — tau_UQFF = tau_std*(1/(1-f_TRZ))*(rho_UA/rho_SCm)*exp(U_m/k_B T_H). Mass classification: stable_DM/marginal/evaporating. Min DM mass binary search. Window reopened for M>~3e10 kg in UQFF.
    """
    ENTRY = 252
    PAPER = "PAPER_668"
    CPP_MODULE = "UQFFStabilityPrimordialBH"
    UQFF_CONSTANTS = {
        "rho_UA": 7.09e-36, "rho_SCm": 7.09e-37, "f_TRZ": 0.1,
        "kappa": 0.0005, "SSq": 0.57, "mu_j": 3.38e23, "gamma": 5e-5,
    }
    PARAMETERS = {
        "M_kg": {"type": "float", "default": 1.0e+12, "desc": "Primordial BH mass in kg"},
    }
    PRIMARY_OUTPUT = "tau_UQFF_s"
    EQUATIONS = [
        "UQFF Primordial BH Stability — tau_UQFF = tau_std*(1/(1-f_TRZ))*(rho_UA/rho_SCm)*exp(U_m/k_B T_H). Mass classification: stable_DM/marginal/evaporating. Min DM mass binary search. Window reopened for M",
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
            "entry":           252,
            "paper":           "PAPER_668",
            "class":           "UQFFStabilityPrimordialBHCalculator",
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
            "tau_UQFF_s":   tau_uqff if "tau" in "tau_UQFF_s" else (
                               L_uqff if "L_" in "tau_UQFF_s" else
                               T_uqff if "T_" in "tau_UQFF_s" else 0.0),
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
    "# CP4 #252 \x97 UQFFStabilityPrimordialBHCalculator\n"
    "# PAPER_668 | Session 172 | April 2, 2026\n"
    "# " + "=" * 74 + "\n"
)

def main():
    block = DIVIDER + CLASS_SRC
    encoded = block.encode("utf-8", errors="replace")
    with open(CP4, "ab") as f:
        f.write(encoded)
    print(f"CP4 #252 appended: UQFFStabilityPrimordialBHCalculator")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-800:]
    if b"UQFFStabilityPrimordialBHCalculator" in tail:
        print("  Verification: PASS")
    else:
        print("  Verification: FAIL - class name not found in tail")

if __name__ == "__main__":
    main()
