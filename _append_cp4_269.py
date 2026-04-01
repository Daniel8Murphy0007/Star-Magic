#!/usr/bin/env python3
"""Append CP4 #269 (UQFFPBHDarkMatterImplicationsCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #269 — UQFFPBHDarkMatterImplicationsCalculator\n# PAPER_685 | Session 173 | appended by _append_cp4_269.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFPBHDarkMatterImplicationsCalculator:\n    """\n    CP4 #269 — PAPER_685 | UQFF PBH dark matter — M_crit shift, f_PBH boost, viable mass window scan\n    CPP_MODULE: UQFFPBHDarkMatterImplications\n    """\n    ENTRY   = 269\n    PAPER   = "PAPER_685"\n    CPP_MODULE = "UQFFPBHDarkMatterImplications"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 1000000000000.0, \'PBH mass kg\'),\n        (\'n_pbh_m3\', \'float\', 1000000.0, \'PBH number density m^-3\'),\n    ]\n    PRIMARY_OUTPUT = "f_pbh_boost"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 1000000000000.0, \'n_pbh_m3\': 1000000.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 1000000000000.0)\n        n_pbh_m3 = params.get(\'n_pbh_m3\', 1000000.0)\n        T_age = 4.34e17\n        M_crit_std  = (HBAR*C**4*T_age/(5120.0*PI*G**2))**(1.0/3.0)\n        tau_ratio   = 1.0/(1.0-F_TRZ)*(RHO_UA/RHO_SCM)\n        M_crit_UQFF = M_crit_std / tau_ratio**(1.0/3.0)\n        f_pbh_boost = tau_ratio**(2.0/3.0)\n        rho_PBH     = mass_kg * n_pbh_m3 * f_pbh_boost\n        viable_std  = mass_kg > M_crit_std\n        viable_UQFF = mass_kg > M_crit_UQFF\n        result = {\'f_pbh_boost\': f_pbh_boost, \'M_crit_std_kg\': M_crit_std,\n                   \'M_crit_UQFF_kg\': M_crit_UQFF, \'rho_PBH_UQFF\': rho_PBH,\n                   \'viable_std\': viable_std, \'viable_UQFF\': viable_UQFF}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFPBHDarkMatterImplicationsCalculator" in existing:
        print(f"SKIP: UQFFPBHDarkMatterImplicationsCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #269 UQFFPBHDarkMatterImplicationsCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFPBHDarkMatterImplicationsCalculator" in tail or "UQFFPBHDarkMatterImplicationsCalculator" in existing:
        print(f"VERIFIED: UQFFPBHDarkMatterImplicationsCalculator present")
    else:
        print(f"WARNING: Could not verify UQFFPBHDarkMatterImplicationsCalculator — check file")

if __name__ == "__main__":
    main()
