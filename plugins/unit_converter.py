"""
Sample UQFF Plugin - Physics Unit Converter
Converts between different unit systems
"""

def process(equation: str) -> dict:
    """Convert units in equation."""
    # Physical constants
    c = 299792458  # m/s
    G = 6.67430e-11  # m^3/kg/s^2
    hbar = 1.054571817e-34  # J*s
    
    # Natural unit conversions
    conversions = {
        "eV_to_J": 1.602176634e-19,
        "J_to_eV": 6.241509074e18,
        "ly_to_m": 9.4607304725808e15,
        "pc_to_m": 3.0856775814913673e16,
        "au_to_m": 1.495978707e11,
        "solar_mass_kg": 1.98892e30,
    }
    
    return {
        "original": equation,
        "conversions_available": list(conversions.keys()),
        "constants": {"c": c, "G": G, "hbar": hbar}
    }
