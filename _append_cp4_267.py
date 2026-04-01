#!/usr/bin/env python3
"""Append CP4 #267 (UQFFHawkingTemperatureModulationCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #267 — UQFFHawkingTemperatureModulationCalculator\n# PAPER_683 | Session 173 | appended by _append_cp4_267.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFHawkingTemperatureModulationCalculator:\n    """\n    CP4 #267 — PAPER_683 | UQFF Hawking T modulation — T_UQFF, modulation factor, Wien peak, Planck spectrum shift\n    CPP_MODULE: UQFFHawkingTemperatureModulation\n    """\n    ENTRY   = 267\n    PAPER   = "PAPER_683"\n    CPP_MODULE = "UQFFHawkingTemperatureModulation"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 1000000000000.0, \'BH mass kg\'),\n        (\'radius_m\', \'float\', 0.0, \'evaluation radius m (0=r_s)\'),\n        (\'time_s\', \'float\', 100000000.0, \'evaluation time s\'),\n    ]\n    PRIMARY_OUTPUT = "T_uqff_K"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 1000000000000.0, \'radius_m\': 0.0, \'time_s\': 100000000.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 1000000000000.0)\n        radius_m = params.get(\'radius_m\', 0.0)\n        time_s = params.get(\'time_s\', 100000000.0)\n        M = mass_kg\n        r = r_s(M) if radius_m<=0 else radius_m\n        TH = T_H(M)\n        Um = U_m(r, time_s)\n        fac = (1.0+F_TRZ)*(1.0-RHO_SCM/RHO_UA)*(1.0+Um/(K_B*TH)) if TH>0 else 1.0\n        T_uqff_K = TH * fac\n        wien_GR   = HBAR*C/(2.82*K_B*TH) if TH>0 else 0\n        wien_UQFF = HBAR*C/(2.82*K_B*T_uqff_K) if T_uqff_K>0 else 0\n        result = {\'T_uqff_K\': T_uqff_K, \'T_H_K\': TH, \'mod_factor\': fac,\n                   \'wien_GR_m\': wien_GR, \'wien_UQFF_m\': wien_UQFF}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFHawkingTemperatureModulationCalculator" in existing:
        print(f"SKIP: UQFFHawkingTemperatureModulationCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #267 UQFFHawkingTemperatureModulationCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFHawkingTemperatureModulationCalculator" in tail or "UQFFHawkingTemperatureModulationCalculator" in existing:
        print(f"VERIFIED: UQFFHawkingTemperatureModulationCalculator present")
    else:
        print(f"WARNING: Could not verify UQFFHawkingTemperatureModulationCalculator — check file")

if __name__ == "__main__":
    main()
