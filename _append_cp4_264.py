#!/usr/bin/env python3
"""Append CP4 #264 (VortexQuantizationCalculator) to CondensedPhysics2.py"""
import os, sys

CP4  = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics2.py"
DIVIDER = '\n\n# ========================================================================\n# CP4 #264 — VortexQuantizationCalculator\n# PAPER_680 | Session 173 | appended by _append_cp4_264.py\n# ========================================================================\n'
CLASS_SRC = 'class VortexQuantizationCalculator:\n    """\n    CP4 #264 — PAPER_680 | UQFF vortex quantization — circulation kappa_v, core radius a_v, energy E_v, Magnus force\n    CPP_MODULE: VortexQuantization\n    """\n    ENTRY   = 264\n    PAPER   = "PAPER_680"\n    CPP_MODULE = "VortexQuantization"\n    UQFF_CONSTANTS = {\n        "RHO_UA":  7.09e-36,\n        "RHO_SCM": 7.09e-37,\n        "F_TRZ":   0.1,\n        "KAPPA":   0.0005,\n        "SSQ":     0.57,\n        "MU_J":    3.38e23,\n        "GAMMA":   5.0e-5 / 86400.0,\n    }\n    PARAMETERS = [\n        (\'winding_number\', \'int\', 1, \'vortex winding number n\'),\n        (\'R_outer_m\', \'float\', 1000000.0, \'outer boundary radius m\'),\n        (\'v_s_m_s\', \'float\', 1000.0, \'superfluid velocity m/s\'),\n    ]\n    PRIMARY_OUTPUT = "vortex_energy"\n    PRIMARY_INPUT  = "winding_number"\n\n    def compute(self, params: dict = None) -> dict:\n        if params is None:\n            params = {}\n        defaults = {\'winding_number\': 1, \'R_outer_m\': 1000000.0, \'v_s_m_s\': 1000.0}\n        for k, v in defaults.items():\n            params.setdefault(k, v)\n        winding_number = params.get(\'winding_number\', 1)\n        R_outer_m = params.get(\'R_outer_m\', 1000000.0)\n        v_s_m_s = params.get(\'v_s_m_s\', 1000.0)\n        H0_SI = 2.184e-18\n        M_UA  = HBAR*H0_SI/C**2\n        G_UA  = 1.0e-10\n        N_UA  = RHO_UA/(HBAR*H0_SI/C**2)\n        n = int(winding_number)\n        kv    = n * 2.0*PI*HBAR/M_UA\n        xi    = HBAR/math.sqrt(2.0*M_UA*G_UA*N_UA)\n        a_v   = xi*math.exp(-n*PI)\n        if a_v>0 and R_outer_m>a_v:\n            Ev = RHO_UA*kv**2/(4.0*PI)*math.log(R_outer_m/a_v)*(RHO_UA/RHO_SCM)\n        else:\n            Ev = 0.0\n        c_UA  = math.sqrt(G_UA*N_UA/M_UA)\n        omega_v_dummy = c_UA/R_outer_m  # rough estimate\n        omega_v_UQFF = omega_v_dummy*(1.0+F_TRZ*c_UA/C)\n        F_Magnus = RHO_UA*kv*v_s_m_s*(1.0+F_TRZ*(RHO_UA/RHO_SCM))\n        vortex_energy = Ev\n        result = {\'vortex_energy\': vortex_energy, \'circulation_kv\': kv,\n                   \'core_radius_m\': a_v, \'magnus_force_UQFF\': F_Magnus,\n                   \'omega_v_UQFF\': omega_v_UQFF}\n        self._last = result\n        return result\n\n    def simulate(self, sweep: list = None, sweep_param: str = None) -> list:\n        if sweep is None:\n            sweep = [{}]\n        results = []\n        for p in sweep:\n            results.append(self.compute(dict(p)))\n        return results\n\n    def add_mod(self, fn) -> None:\n        if not hasattr(self, "_mods"):\n            self._mods = []\n        self._mods.append(fn)\n\n    def update_from_file(self, filepath: str) -> None:\n        import os\n        if not os.path.isfile(filepath):\n            return\n        with open(filepath, "r", encoding="utf-8") as f:\n            for line in f:\n                line = line.strip()\n                if "=" in line:\n                    k, v = line.split("=", 1)\n                    try:\n                        object.__setattr__(self, k.strip(), float(v.strip()))\n                    except ValueError:\n                        pass\n'

def main():
    if not os.path.isfile(CP4):
        print(f"ERROR: Cannot find {CP4}", file=sys.stderr); sys.exit(1)
    # Check if already appended
    with open(CP4, "rb") as f:
        existing = f.read().decode("utf-8", errors="replace")
    if "class VortexQuantizationCalculator" in existing:
        print(f"SKIP: VortexQuantizationCalculator already present in CondensedPhysics2.py")
        return
    block = DIVIDER + CLASS_SRC
    with open(CP4, "ab") as f:
        f.write(block.encode("utf-8"))
    print(f"APPENDED: CP4 #264 VortexQuantizationCalculator -> CondensedPhysics2.py")
    # Verify
    with open(CP4, "rb") as f:
        tail = f.read()[-2000:].decode("utf-8", errors="replace")
    if "class VortexQuantizationCalculator" in tail or "VortexQuantizationCalculator" in existing:
        print(f"VERIFIED: VortexQuantizationCalculator present")
    else:
        print(f"WARNING: Could not verify VortexQuantizationCalculator — check file")

if __name__ == "__main__":
    main()
