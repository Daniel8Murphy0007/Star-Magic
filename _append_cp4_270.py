#!/usr/bin/env python3
"""Append CP4 #270 (UQFFModulationForM87Calculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #270 — UQFFModulationForM87Calculator\n# PAPER_686 | Session 173 | appended by _append_cp4_270.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFModulationForM87Calculator:\n    """\n    CP4 #270 — PAPER_686 | UQFF modulation for M87* — T_UQFF, shadow_UQFF, jet power, ring brightness, T_acc\n    CPP_MODULE: UQFFModulationForM87\n    """\n    ENTRY   = 270\n    PAPER   = "PAPER_686"\n    CPP_MODULE = "UQFFModulationForM87"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 1.293e+40, \'M87* mass kg (6.5e9 Msun)\'),\n        (\'T_acc_K\', \'float\', 10000000.0, \'accretion disk temperature K\'),\n    ]\n    PRIMARY_OUTPUT = "T_uqff_K"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 1.293e+40, \'T_acc_K\': 10000000.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 1.293e+40)\n        T_acc_K = params.get(\'T_acc_K\', 10000000.0)\n        M = mass_kg\n        TH    = T_H(M)\n        T_uqff_K = TH*(1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)\n        r_sh  = 3.0*math.sqrt(3.0)*G*M/C**2\n        r_sh_UQFF = r_sh*math.sqrt(1.0+F_TRZ*(RHO_UA/RHO_SCM))\n        P_JET_GR = 1.0e37\n        P_jet_UQFF = P_JET_GR*(1.0+F_TRZ)*math.sqrt(RHO_UA/RHO_SCM)\n        brightness_ratio = (RHO_UA/RHO_SCM)**(F_TRZ/4.0)\n        r = r_s(M); Um = U_m(r, 1e8)\n        T_acc_UQFF = T_acc_K*(1.0+Um/(K_B*TH)) if TH>0 else T_acc_K\n        result = {\'T_uqff_K\': T_uqff_K, \'T_H_K\': TH,\n                   \'shadow_GR_m\': r_sh, \'shadow_UQFF_m\': r_sh_UQFF,\n                   \'jet_power_UQFF_W\': P_jet_UQFF, \'brightness_ratio\': brightness_ratio,\n                   \'T_acc_UQFF_K\': T_acc_UQFF}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFModulationForM87Calculator" in existing:
        print(f"SKIP: UQFFModulationForM87Calculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #270 UQFFModulationForM87Calculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFModulationForM87Calculator" in tail or "UQFFModulationForM87Calculator" in existing:
        print(f"VERIFIED: UQFFModulationForM87Calculator present")
    else:
        print(f"WARNING: Could not verify UQFFModulationForM87Calculator — check file")

if __name__ == "__main__":
    main()
