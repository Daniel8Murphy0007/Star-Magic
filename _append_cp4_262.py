#!/usr/bin/env python3
"""Append CP4 #262 (LISAVsLIGOComparisonsCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #262 — LISAVsLIGOComparisonsCalculator\n# PAPER_678 | Session 173 | appended by _append_cp4_262.py\n# ========================================================================\n'
CLASS_SRC = 'class LISAVsLIGOComparisonsCalculator:\n    """\n    CP4 #262 — PAPER_678 | UQFF LISA vs LIGO suppression comparison — R_supp, crossover frequency, Sn_UQFF\n    CPP_MODULE: LISAVsLIGOComparisons\n    """\n    ENTRY   = 262\n    PAPER   = "PAPER_678"\n    CPP_MODULE = "LISAVsLIGOComparisons"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 3.978e+31, \'BH mass kg (20 Msun)\'),\n        (\'frequency_hz\', \'float\', 150.0, \'GW frequency Hz\'),\n        (\'T_eff_K\', \'float\', 2.73, \'effective temperature K\'),\n    ]\n    PRIMARY_OUTPUT = "R_supp_ligo"\n    PRIMARY_INPUT  = "frequency_hz"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 3.978e+31, \'frequency_hz\': 150.0, \'T_eff_K\': 2.73}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 3.978e+31)\n        frequency_hz = params.get(\'frequency_hz\', 150.0)\n        T_eff_K = params.get(\'T_eff_K\', 2.73)\n        L_LISA = 2.5e9\n        r  = r_s(mass_kg); TH = T_H(mass_kg)\n        S1 = S_SCm(mass_kg, TH, r)\n        S2 = S_Um_f(mass_kg, r, 1e8, frequency_hz)\n        R_supp_ligo = (1.0-F_TRZ)*S1*S2\n        S_UA = max(0.0, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))\n        R_supp_lisa = (1.0-F_TRZ)*S_UA*S1\n        Um = U_m(r, 1e8)\n        S_UA2 = max(1e-300, 1.0 - RHO_UA*L_LISA/(K_B*T_eff_K))\n        if Um>0 and S_UA2>0:\n            f_cross = -math.log(S_UA2)*C**2/(2.0*PI*Um)\n        else:\n            f_cross = 0.0\n        result = {\'R_supp_ligo\': R_supp_ligo, \'R_supp_lisa\': R_supp_lisa,\n                   \'crossover_freq_Hz\': f_cross}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class LISAVsLIGOComparisonsCalculator" in existing:
        print(f"SKIP: LISAVsLIGOComparisonsCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #262 LISAVsLIGOComparisonsCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class LISAVsLIGOComparisonsCalculator" in tail or "LISAVsLIGOComparisonsCalculator" in existing:
        print(f"VERIFIED: LISAVsLIGOComparisonsCalculator present")
    else:
        print(f"WARNING: Could not verify LISAVsLIGOComparisonsCalculator — check file")

if __name__ == "__main__":
    main()
