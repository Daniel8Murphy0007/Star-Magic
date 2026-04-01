#!/usr/bin/env python3
"""Append CP4 #271 (M87MassEvolutionSimulationCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #271 — M87MassEvolutionSimulationCalculator\n# PAPER_687 | Session 173 | appended by _append_cp4_271.py\n# ========================================================================\n'
CLASS_SRC = 'class M87MassEvolutionSimulationCalculator:\n    """\n    CP4 #271 — PAPER_687 | M87* mass evolution — Bondi-UQFF accretion, Hawking evaporation, jet mass loss, Hubble-time RK4\n    CPP_MODULE: M87MassEvolutionSimulation\n    """\n    ENTRY   = 271\n    PAPER   = "PAPER_687"\n    CPP_MODULE = "M87MassEvolutionSimulation"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'mass_kg\', \'float\', 1.293e+40, \'M87* mass kg\'),\n        (\'rho_ism_kg_m3\', \'float\', 1.67e-25, \'ISM density kg/m^3\'),\n        (\'T_ism_K\', \'float\', 10000000.0, \'ISM temperature K\'),\n        (\'dt_s\', \'float\', 10000000000000.0, \'RK4 time step s\'),\n    ]\n    PRIMARY_OUTPUT = "dM_dt_total_kg_s"\n    PRIMARY_INPUT  = "mass_kg"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'mass_kg\': 1.293e+40, \'rho_ism_kg_m3\': 1.67e-25, \'T_ism_K\': 10000000.0, \'dt_s\': 10000000000000.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        mass_kg = params.get(\'mass_kg\', 1.293e+40)\n        rho_ism_kg_m3 = params.get(\'rho_ism_kg_m3\', 1.67e-25)\n        T_ism_K = params.get(\'T_ism_K\', 10000000.0)\n        dt_s = params.get(\'dt_s\', 10000000000000.0)\n        M = mass_kg\n        # Bondi-UQFF accretion\n        c_s = math.sqrt(5.0/3.0*K_B*T_ism_K/1.673e-27)\n        Mdot_Bondi = 4.0*PI*G**2*M**2*rho_ism_kg_m3/c_s**3\n        rho_eff = rho_ism_kg_m3 + RHO_UA - RHO_SCM\n        Mdot_UQFF = Mdot_Bondi*(rho_eff/rho_ism_kg_m3)*(1.0+F_TRZ)\n        # Evaporation\n        TH = T_H(M); r = r_s(M); Um = U_m(r, 1e8)\n        dM_std = -HBAR*C**4/(15360.0*PI*G**2*M**2) if M>0 else 0\n        arg = Um/(K_B*TH) if TH>0 else 0\n        dM_evap = dM_std*(1.0-F_TRZ)*(RHO_SCM/RHO_UA)*math.exp(-min(arg,700))\n        # Jet loss\n        rho_disk = rho_ism_kg_m3*1e6\n        B_eq = math.sqrt(8.0*PI*rho_disk*C**2)\n        if B_eq>1e8: B_eq=1e8\n        Omega_H = C/(4.0*G*M/C**2)*C if M>0 else 0\n        rs_m = r_s(M)\n        P_BZ = 0.01*B_eq**2*rs_m**2*C*Omega_H**2/(4.0*PI*C)\n        P_jet_UQFF = P_BZ*(1.0+F_TRZ)*math.sqrt(RHO_UA/RHO_SCM)\n        dM_jet = -P_jet_UQFF/C**2\n        dM_dt_total_kg_s = Mdot_UQFF + dM_evap + dM_jet\n        result = {\'dM_dt_total_kg_s\': dM_dt_total_kg_s, \'Mdot_Bondi_UQFF\': Mdot_UQFF,\n                   \'dM_evap\': dM_evap, \'dM_jet\': dM_jet, \'P_jet_W\': P_jet_UQFF}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class M87MassEvolutionSimulationCalculator" in existing:
        print(f"SKIP: M87MassEvolutionSimulationCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #271 M87MassEvolutionSimulationCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class M87MassEvolutionSimulationCalculator" in tail or "M87MassEvolutionSimulationCalculator" in existing:
        print(f"VERIFIED: M87MassEvolutionSimulationCalculator present")
    else:
        print(f"WARNING: Could not verify M87MassEvolutionSimulationCalculator — check file")

if __name__ == "__main__":
    main()
