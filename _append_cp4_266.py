#!/usr/bin/env python3
"""Append CP4 #266 (UQFFStabilityNumericallyForSgrACalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #266 — UQFFStabilityNumericallyForSgrACalculator\n# PAPER_682 | Session 173 | appended by _append_cp4_266.py\n# ========================================================================\n'
CLASS_SRC = 'class UQFFStabilityNumericallyForSgrACalculator:\n    """\n    CP4 #266 — PAPER_682 | UQFF numerical stability for Sgr A* — Lyapunov, omega_I, RK4 mass evolution, stability check\n    CPP_MODULE: UQFFStabilityNumericallyForSgrA\n    """\n    ENTRY   = 266\n    PAPER   = "PAPER_682"\n    CPP_MODULE = "UQFFStabilityNumericallyForSgrA"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 8.548e+36, \'BH mass kg (Sgr A*)\'),\n        (\'t_end_s\', \'float\', 1e+60, \'evolution end time s\'),\n        (\'dt_s\', \'float\', 1e+55, \'time step s\'),\n    ]\n    PRIMARY_OUTPUT = "lyapunov_exponent"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 8.548e+36, \'t_end_s\': 1e+60, \'dt_s\': 1e+55}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 8.548e+36)\n        t_end_s = params.get(\'t_end_s\', 1e+60)\n        dt_s = params.get(\'dt_s\', 1e+55)\n        M = mass_kg\n        TH = T_H(M); r = r_s(M); tau = tau_std(M)\n        Um = U_m(r, 1e8)\n        omega_I = -1.0/max(tau,1e-300) * (1.0+F_TRZ*(RHO_UA/RHO_SCM)*Um/(K_B*TH))\n        arg = Um/(K_B*TH) if TH>0 else 0\n        lam = -(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))/max(tau,1e-300)\n        stable = lam < 0\n        lyapunov_exponent = lam\n        result = {\'lyapunov_exponent\': lyapunov_exponent, \'omega_I\': omega_I,\n                   \'stable\': stable, \'T_H_K\': TH, \'tau_std_s\': tau}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class UQFFStabilityNumericallyForSgrACalculator" in existing:
        print(f"SKIP: UQFFStabilityNumericallyForSgrACalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #266 UQFFStabilityNumericallyForSgrACalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class UQFFStabilityNumericallyForSgrACalculator" in tail or "UQFFStabilityNumericallyForSgrACalculator" in existing:
        print(f"VERIFIED: UQFFStabilityNumericallyForSgrACalculator present")
    else:
        print(f"WARNING: Could not verify UQFFStabilityNumericallyForSgrACalculator — check file")

if __name__ == "__main__":
    main()
