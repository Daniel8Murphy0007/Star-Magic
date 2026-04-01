#!/usr/bin/env python3
"""Append CP4 #268 (UQFFPrimordialBHEvaporationCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #268 — UQFFPrimordialBHEvaporationCalculator\n# PAPER_684 | Session 173 | appended by _append_cp4_268.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFPrimordialBHEvaporationCalculator:\n    """\n    CP4 #268 — PAPER_684 | UQFF primordial BH evaporation — dM/dt suppressed, tau_UQFF, M(t_form), RK4 evolution\n    CPP_MODULE: UQFFPrimordialBHEvaporation\n    """\n    ENTRY   = 268\n    PAPER   = "PAPER_684"\n    CPP_MODULE = "UQFFPrimordialBHEvaporation"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 1000000000000.0, \'PBH mass kg\'),\n        (\'t_form_s\', \'float\', 1e-23, \'formation time s\'),\n        (\'dt_s\', \'float\', 1e+50, \'RK4 time step s\'),\n    ]\n    PRIMARY_OUTPUT = "tau_uqff_s"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 1000000000000.0, \'t_form_s\': 1e-23, \'dt_s\': 1e+50}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 1000000000000.0)\n        t_form_s = params.get(\'t_form_s\', 1e-23)\n        dt_s = params.get(\'dt_s\', 1e+50)\n        M = mass_kg\n        TH = T_H(M); r = r_s(M); tau = tau_std(M)\n        Um = U_m(r, 1e8)\n        dM_std = -HBAR*C**4/(15360.0*PI*G**2*M**2) if M>0 else 0\n        arg = Um/(K_B*TH) if TH>0 else 0\n        dM_UQFF = dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))\n        tau_UQFF_s = tau/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)*math.exp(min(arg,700))\n        # M at formation from t_form\n        r_f = 2.998e8*t_form_s/2.0\n        rho_rad = 3.0/(32.0*PI*G*t_form_s**2)\n        M_form = (4.0/3.0)*PI*r_f**3*rho_rad\n        result = {\'tau_uqff_s\': tau_UQFF_s, \'tau_std_s\': tau,\n                   \'dM_dt_UQFF\': dM_UQFF, \'M_formation_kg\': M_form,\n                   \'suppression_factor\': (1.0-F_TRZ)*RHO_SCM/RHO_UA}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFPrimordialBHEvaporationCalculator" in existing:
        print(f"SKIP: UQFFPrimordialBHEvaporationCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #268 UQFFPrimordialBHEvaporationCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFPrimordialBHEvaporationCalculator" in tail or "UQFFPrimordialBHEvaporationCalculator" in existing:
        print(f"VERIFIED: UQFFPrimordialBHEvaporationCalculator present")
    else:
        print(f"WARNING: Could not verify UQFFPrimordialBHEvaporationCalculator — check file")

if __name__ == "__main__":
    main()
