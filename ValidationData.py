#!/usr/bin/env python3
"""
ValidationData.py - Validation Dataset Wrapper
==============================================

Links to CondensedPhysics_Validation.py for complete validation URLs and references.

PURPOSE:
    - Provide simplified interface to validation data
    - Import validation URLs from CondensedPhysics_Validation.py
    - Organize validation datasets by physics domain

DATA FLOW:
    CondensedPhysics_Validation.py (6000+ validation URLs)
                    ↓
    ValidationData.py (this file - simplified access)
                    ↓
    CondensedPhysics.py (calculator uses for citations)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Created: February 11, 2026
"""

# Import all validation data from CondensedPhysics_Validation.py
try:
    from CondensedPhysics_Validation import *
    print("✅ Loaded validation data from CondensedPhysics_Validation.py")
    print(f"   - Nuclear Binding: {len(NUCLEAR_BINDING_VALIDATION)} sources")
    print(f"   - UQFF Validation: Available")
    print(f"   - Observational Data: Available")
except ImportError as e:
    print(f"⚠️  Warning: Could not import CondensedPhysics_Validation.py: {e}")
    NUCLEAR_BINDING_VALIDATION = {}


# ═══════════════════════════════════════════════════════════════════════════════
# SIMPLIFIED VALIDATION ACCESS
# ═══════════════════════════════════════════════════════════════════════════════

def get_validation_url(domain: str, reference: str) -> str:
    """
    Get validation URL for a specific physics domain and reference.
    
    Args:
        domain: Physics domain (e.g., 'nuclear', 'cosmology', 'magnetar')
        reference: Reference name (e.g., 'pdg_2025', 'arxiv_2408_04231')
    
    Returns:
        URL string or empty string if not found
    """
    if domain == 'nuclear' and reference in NUCLEAR_BINDING_VALIDATION:
        return NUCLEAR_BINDING_VALIDATION[reference].get('url', '')
    
    # Add more domains as needed
    return ''


def list_validation_sources(domain: str = None) -> list:
    """
    List all available validation sources.
    
    Args:
        domain: Optional domain filter
    
    Returns:
        List of source names
    """
    if domain == 'nuclear':
        return list(NUCLEAR_BINDING_VALIDATION.keys())
    
    # Return all if no domain specified
    all_sources = []
    all_sources.extend(NUCLEAR_BINDING_VALIDATION.keys())
    return all_sources


# ═══════════════════════════════════════════════════════════════════════════════
# QUICK REFERENCE DATASETS
# ═══════════════════════════════════════════════════════════════════════════════

QUICK_REFERENCES = {
    'pdg_2025': {
        'url': 'https://pdg.lbl.gov/2025/reviews/rpp2024-rev-passage-particles-matter.pdf',
        'description': 'Particle Data Group 2025 - Nuclear binding energies',
        'category': 'Nuclear Physics'
    },
    'ensdf_nndc': {
        'url': 'https://www.nndc.bnl.gov/ensdf/',
        'description': 'Evaluated Nuclear Structure Data File',
        'category': 'Nuclear Database'
    },
    'gw170817_ligo': {
        'url': 'https://arxiv.org/abs/1710.05832',
        'description': 'LIGO/Virgo GW170817 detection paper',
        'category': 'Gravitational Waves'
    },
    'simbad': {
        'url': 'http://simbad.u-strasbg.fr/simbad/',
        'description': 'SIMBAD Astronomical Database',
        'category': 'Astrophysics'
    },
    'ned': {
        'url': 'https://ned.ipac.caltech.edu/',
        'description': 'NASA/IPAC Extragalactic Database',
        'category': 'Astrophysics'
    }
}


# ═══════════════════════════════════════════════════════════════════════════════
# VALIDATION BENCHMARKS
# ═══════════════════════════════════════════════════════════════════════════════

VALIDATION_BENCHMARKS = {
    'nuclear_binding': {
        'Pb_206_ground_state': 1.6 * 10**-12,  # J (10 MeV)
        'Fe_56_binding_per_nucleon': 8.79 * 1.602e-13,  # J
        'U_238_fission_energy': 200 * 1.602e-13,  # J
    },
    'cosmology': {
        'Lambda': 1.1e-52,  # m^-2
        'H_0': 70e3 / 3.086e22,  # s^-1 (70 km/s/Mpc)
        'Omega_m': 0.315,
        'Omega_Lambda': 0.685,
    },
    'magnetar': {
        'B_crit': 4.4e13,  # T (QED critical field)
        'B_sgr1745': 1e15,  # T (SGR 1745-2900)
        'P_typical': 2-12,  # seconds (spin period range)
    },
    'gw170817': {
        'M_total': 2.8 * 1.989e30,  # kg
        'M_ej_total': 0.05 * 1.989e30,  # kg
        'v_ej_c': 0.1,  # fraction of speed of light
        'distance_Mpc': 40,
    }
}


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("ValidationData.py - Validation Dataset Wrapper")
    print("=" * 70)
    
    # List available sources
    print("\nNuclear physics validation sources:")
    sources = list_validation_sources('nuclear')
    for source in sources[:5]:  # Show first 5
        print(f"  - {source}")
    
    print(f"\n  ... and {len(sources) - 5} more sources")
    
    # Get specific URL
    url = get_validation_url('nuclear', 'pdg_2025')
    print(f"\nPDG 2025 URL: {url}")
    
    # Show quick references
    print("\nQuick references:")
    for ref_id, ref_data in QUICK_REFERENCES.items():
        print(f"  {ref_id}: {ref_data['description']}")
    
    # Show validation benchmarks
    print("\nValidation benchmarks:")
    print(f"  Pb-206 ground state: {VALIDATION_BENCHMARKS['nuclear_binding']['Pb_206_ground_state']:.2e} J")
    print(f"  Hubble constant: {VALIDATION_BENCHMARKS['cosmology']['H_0']:.2e} s^-1")
    print(f"  QED critical field: {VALIDATION_BENCHMARKS['magnetar']['B_crit']:.2e} T")
    print(f"  GW170817 total mass: {VALIDATION_BENCHMARKS['gw170817']['M_total']:.2e} kg")
