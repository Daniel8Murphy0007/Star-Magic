#!/usr/bin/env python3
"""
IPData.py - UQFF Input Parameter Data Storage
==============================================

Stores fetched API data and user manual inputs for QCalc.py processing.

ARCHITECTURE:
    APIFetch.py → IPData.py (this file) → QCalc.py → OPData.py
    
STORAGE FORMAT:
    - Dictionary-based storage with query IDs
    - Automatic timestamping
    - Source tracking (which API provided data)
    - Parameter validation before storage

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

from dataclasses import dataclass, field, asdict
from typing import Dict, List, Any, Optional
from datetime import datetime
import json
import os

# ═══════════════════════════════════════════════════════════════════════════════
# INPUT PARAMETER SCHEMA
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class InputParameters:
    """
    Standardized input parameters for UQFF calculations.
    
    ALL parameters are optional - QCalc.py determines which equations
    are solvable based on available parameters.
    """
    # Query Metadata
    query_id: str = ""
    query_name: str = ""
    timestamp: str = ""
    sources: List[str] = field(default_factory=list)
    
    # ───────────────────────────────────────────────────────────────────────────
    # MASS PARAMETERS (all in kg)
    # ───────────────────────────────────────────────────────────────────────────
    M: Optional[float] = None               # Primary mass
    M_companion: Optional[float] = None     # Companion mass (binaries)
    M_enclosed: Optional[float] = None      # Enclosed mass (galactic dynamics)
    M_bh: Optional[float] = None            # Black hole mass
    M_halo: Optional[float] = None          # Dark matter halo mass
    M_bulge: Optional[float] = None         # Bulge mass
    M_disk: Optional[float] = None          # Disk mass
    m_eff: Optional[float] = None           # Effective/reduced mass
    
    # ───────────────────────────────────────────────────────────────────────────
    # DISTANCE/RADIUS PARAMETERS (all in m)
    # ───────────────────────────────────────────────────────────────────────────
    r: Optional[float] = None               # Primary radius/distance
    R: Optional[float] = None               # Object radius
    d: Optional[float] = None               # Distance to observer
    d_g: Optional[float] = None             # Distance to galactic center
    R_s: Optional[float] = None             # Schwarzschild radius
    R_isco: Optional[float] = None          # Innermost stable circular orbit
    R_hill: Optional[float] = None          # Hill radius
    a: Optional[float] = None               # Semi-major axis
    r_core: Optional[float] = None          # Core radius
    r_half: Optional[float] = None          # Half-light radius
    scale_height: Optional[float] = None    # Disk scale height
    
    # ───────────────────────────────────────────────────────────────────────────
    # TEMPERATURE PARAMETERS (all in K)
    # ───────────────────────────────────────────────────────────────────────────
    T: Optional[float] = None               # Temperature
    T_eff: Optional[float] = None           # Effective temperature
    T_critical: Optional[float] = None      # Critical temperature
    T_disk: Optional[float] = None          # Disk temperature
    T_corona: Optional[float] = None        # Corona temperature
    T_cmb: float = 2.725                    # CMB temperature (constant)
    
    # ───────────────────────────────────────────────────────────────────────────
    # MAGNETIC PARAMETERS (B in T, μ in J/T)
    # ───────────────────────────────────────────────────────────────────────────
    B: Optional[float] = None               # Magnetic field strength
    B_crit: Optional[float] = None          # Critical magnetic field
    B_dipole: Optional[float] = None        # Dipole magnetic field
    mu: Optional[float] = None              # Magnetic moment
    Phi_B: Optional[float] = None           # Magnetic flux (Wb)
    
    # ───────────────────────────────────────────────────────────────────────────
    # ENERGY/LUMINOSITY PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    L: Optional[float] = None               # Luminosity (W)
    L_bol: Optional[float] = None           # Bolometric luminosity
    L_X: Optional[float] = None             # X-ray luminosity
    L_Edd: Optional[float] = None           # Eddington luminosity
    E: Optional[float] = None               # Energy (J)
    E_bind: Optional[float] = None          # Binding energy
    Delta: Optional[float] = None           # Energy gap (superconducting)
    
    # ───────────────────────────────────────────────────────────────────────────
    # VELOCITY/ROTATION PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    v: Optional[float] = None               # Velocity (m/s)
    v_rot: Optional[float] = None           # Rotational velocity
    v_disp: Optional[float] = None          # Velocity dispersion σ
    v_rad: Optional[float] = None           # Radial velocity
    v_wind: Optional[float] = None          # Stellar wind velocity (m/s)
    omega: Optional[float] = None           # Angular velocity (rad/s)
    Omega_g: Optional[float] = None         # Galactic angular velocity
    P: Optional[float] = None               # Period (s)
    P_spin: Optional[float] = None          # Spin period
    P_orb: Optional[float] = None           # Orbital period
    
    # ───────────────────────────────────────────────────────────────────────────
    # DENSITY PARAMETERS (all in kg/m³)
    # ───────────────────────────────────────────────────────────────────────────
    rho: Optional[float] = None             # Density
    rho_central: Optional[float] = None     # Central density
    rho_crit: Optional[float] = None        # Critical density
    rho_wind: Optional[float] = None        # Stellar wind density (kg/m³)
    rho_fluid: Optional[float] = None       # ISM/fluid density (kg/m³)
    n_e: Optional[float] = None             # Electron number density (m⁻³)
    n_H: Optional[float] = None             # Hydrogen number density
    
    # ───────────────────────────────────────────────────────────────────────────
    # COSMOLOGICAL PARAMETERS (dimensionless unless noted)
    # ───────────────────────────────────────────────────────────────────────────
    z: Optional[float] = None               # Redshift
    H_0: float = 67.4e3 / 3.086e22          # Hubble constant (s⁻¹)
    Omega_m: float = 0.315                  # Matter density parameter
    Omega_Lambda: float = 0.685             # Dark energy density parameter
    Omega_b: float = 0.049                  # Baryon density parameter
    
    # ───────────────────────────────────────────────────────────────────────────
    # ACCRETION/FLOW PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    M_dot: Optional[float] = None           # Accretion rate (kg/s)
    M_dot_Edd: Optional[float] = None       # Eddington accretion rate
    SFR: Optional[float] = None             # Star formation rate (M_sun/yr)
    eta_acc: Optional[float] = None         # Accretion efficiency
    
    # ───────────────────────────────────────────────────────────────────────────
    # UQFF-SPECIFIC PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    psi: complex = 1.0 + 0j                 # Order parameter
    psi_magnitude: Optional[float] = None   # |ψ|²
    psi_integral: Optional[float] = None    # Wavefunction integral ∫|ψ|²
    xi: Optional[float] = None              # Coherence length
    lambda_L: Optional[float] = None        # London penetration depth
    Q: complex = 0j                         # Charge parameter
    k_coupling: float = 1.0                 # Coupling constant
    beta_coupling: float = 0.603            # Buoyancy coupling β_i
    eta_coupling: float = 1e-22             # Aether coupling η
    kappa: float = 0.0005                   # κ calibration constant
    SSq: float = 0.57                       # [SSq] constant
    
    # ───────────────────────────────────────────────────────────────────────────
    # DECAY/EVOLUTION TIMESCALES (Wolfram source14/source15/source16)
    # ───────────────────────────────────────────────────────────────────────────
    tau_B: Optional[float] = None           # Magnetic decay timescale (s)
    tau_Omega: Optional[float] = None       # Spin-down timescale (s)
    tau_acc: Optional[float] = None         # Accretion timescale (s)
    tau_SF: Optional[float] = None          # Star formation timescale (s)
    
    # ───────────────────────────────────────────────────────────────────────────
    # QUANTUM UNCERTAINTY PARAMETERS (Wolfram source14)
    # ───────────────────────────────────────────────────────────────────────────
    delta_x: Optional[float] = None         # Position uncertainty (m)
    delta_p: Optional[float] = None         # Momentum uncertainty (kg·m/s)
    
    # ───────────────────────────────────────────────────────────────────────────
    # SURFACE/VELOCITY PARAMETERS (Wolfram source14/source15)
    # ───────────────────────────────────────────────────────────────────────────
    v_surf: Optional[float] = None          # Surface velocity (m/s)
    precession_angle: Optional[float] = None  # Precession angle (rad)
    
    # ───────────────────────────────────────────────────────────────────────────
    # OBSERVATIONAL PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    mag_V: Optional[float] = None           # Visual magnitude
    mag_B: Optional[float] = None           # Blue magnitude
    color_BV: Optional[float] = None        # B-V color
    spectral_type: Optional[str] = None     # Spectral classification
    metallicity: Optional[float] = None     # [Fe/H]
    age: Optional[float] = None             # Age (s)
    parallax: Optional[float] = None        # Parallax (arcsec)
    
    # ───────────────────────────────────────────────────────────────────────────
    # MORPHOLOGICAL PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    inclination: Optional[float] = None     # Inclination angle (rad)
    eccentricity: Optional[float] = None    # Orbital eccentricity
    ellipticity: Optional[float] = None     # Shape ellipticity
    position_angle: Optional[float] = None  # Position angle (rad)
    n_sersic: Optional[float] = None        # Sersic index
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary, excluding None values."""
        return {k: v for k, v in asdict(self).items() if v is not None}
    
    def get_available_params(self) -> List[str]:
        """Return list of parameter names that have values."""
        return [k for k, v in asdict(self).items() 
                if v is not None and k not in ['query_id', 'query_name', 'timestamp', 'sources']]


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT DATA STORE
# ═══════════════════════════════════════════════════════════════════════════════

class InputDataStore:
    """
    Persistent storage for input parameters fetched from APIs.
    
    Stores all query results for:
    - Recall during calculation
    - Audit trail of data sources
    - Batch processing
    - Historical comparison
    """
    
    def __init__(self, storage_file: str = "input_data_store.json"):
        self.storage_file = storage_file
        self.queries: Dict[str, InputParameters] = {}
        self._load()
    
    def _load(self):
        """Load existing data from file."""
        if os.path.exists(self.storage_file):
            try:
                with open(self.storage_file, 'r') as f:
                    data = json.load(f)
                for qid, params in data.items():
                    self.queries[qid] = InputParameters(**params)
            except Exception as e:
                print(f"Warning: Could not load input data: {e}")
    
    def _save(self):
        """Save data to file."""
        try:
            data = {qid: p.to_dict() for qid, p in self.queries.items()}
            with open(self.storage_file, 'w') as f:
                json.dump(data, f, indent=2, default=str)
        except Exception as e:
            print(f"Warning: Could not save input data: {e}")
    
    def store(self, params: InputParameters) -> str:
        """
        Store input parameters from API fetch.
        
        Args:
            params: InputParameters dataclass with fetched data
            
        Returns:
            query_id for recall
        """
        # Generate query ID if not provided
        if not params.query_id:
            params.query_id = f"Q_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{params.query_name}"
        
        # Set timestamp
        params.timestamp = datetime.now().isoformat()
        
        # Store
        self.queries[params.query_id] = params
        self._save()
        
        return params.query_id
    
    def recall(self, query_id: str) -> Optional[InputParameters]:
        """
        Recall stored parameters by query ID.
        
        Args:
            query_id: The ID returned from store()
            
        Returns:
            InputParameters or None if not found
        """
        return self.queries.get(query_id)
    
    def recall_by_name(self, name: str) -> List[InputParameters]:
        """
        Recall all queries for a given object name.
        
        Args:
            name: Object name to search for
            
        Returns:
            List of matching InputParameters
        """
        return [p for p in self.queries.values() 
                if name.lower() in p.query_name.lower()]
    
    def get_latest(self, name: str) -> Optional[InputParameters]:
        """Get the most recent query for an object."""
        matches = self.recall_by_name(name)
        if matches:
            return max(matches, key=lambda p: p.timestamp)
        return None
    
    def list_queries(self) -> List[Dict[str, str]]:
        """List all stored queries with metadata."""
        return [
            {
                'query_id': qid,
                'name': p.query_name,
                'timestamp': p.timestamp,
                'sources': p.sources,
                'param_count': len(p.get_available_params())
            }
            for qid, p in self.queries.items()
        ]
    
    def merge_params(self, base: InputParameters, update: InputParameters) -> InputParameters:
        """
        Merge two parameter sets (update fills gaps in base).
        
        Used when multiple APIs provide complementary data.
        """
        base_dict = base.to_dict()
        update_dict = update.to_dict()
        
        # Update only where base is None
        for key, value in update_dict.items():
            if base_dict.get(key) is None and value is not None:
                base_dict[key] = value
        
        # Merge sources
        all_sources = list(set(base.sources + update.sources))
        base_dict['sources'] = all_sources
        
        return InputParameters(**base_dict)
    
    def validate(self, params: InputParameters) -> Dict[str, Any]:
        """
        Validate input parameters for physical reasonableness.
        
        Returns:
            Dict with 'valid': bool and 'warnings': list
        """
        warnings = []
        
        # Mass checks
        if params.M is not None and params.M <= 0:
            warnings.append("Mass M must be positive")
        
        # Distance checks
        if params.r is not None and params.r <= 0:
            warnings.append("Radius/distance r must be positive")
        
        # Temperature checks
        if params.T is not None and params.T < 0:
            warnings.append("Temperature T cannot be negative (0 K minimum)")
        
        # Redshift checks
        if params.z is not None and params.z < -1:
            warnings.append("Redshift z must be > -1")
        
        # Black hole mass vs radius consistency
        if params.M_bh is not None and params.R_s is not None:
            G = 6.67430e-11
            c = 299792458
            expected_R_s = 2 * G * params.M_bh / c**2
            if abs(params.R_s - expected_R_s) / expected_R_s > 0.1:
                warnings.append(f"R_s inconsistent with M_bh (expected {expected_R_s:.2e})")
        
        return {
            'valid': len(warnings) == 0,
            'warnings': warnings,
            'param_count': len(params.get_available_params())
        }
    
    def export_csv(self, filename: str = "input_data_export.csv"):
        """Export all queries to CSV."""
        if not self.queries:
            return
        
        # Get all possible keys
        all_keys = set()
        for p in self.queries.values():
            all_keys.update(p.to_dict().keys())
        
        all_keys = sorted(all_keys)
        
        with open(filename, 'w', newline='') as f:
            import csv
            writer = csv.DictWriter(f, fieldnames=['query_id'] + list(all_keys))
            writer.writeheader()
            for qid, p in self.queries.items():
                row = {'query_id': qid}
                row.update(p.to_dict())
                writer.writerow(row)
    
    def clear(self):
        """Clear all stored queries."""
        self.queries = {}
        self._save()


# ═══════════════════════════════════════════════════════════════════════════════
# MANUAL INPUT HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def create_manual_input(name: str, **kwargs) -> InputParameters:
    """
    Create InputParameters from manual user input.
    
    Example:
        params = create_manual_input(
            "My Star",
            M=2e30,      # 1 solar mass
            T=5778,      # Solar temperature
            r=1.5e11     # 1 AU
        )
    """
    params = InputParameters(
        query_name=name,
        sources=['manual_input'],
        **kwargs
    )
    return params


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

# Global input data store instance
INPUT_STORE = InputDataStore()


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def store_input(params: InputParameters) -> str:
    """Store input parameters and return query ID."""
    return INPUT_STORE.store(params)

def recall_input(query_id: str) -> Optional[InputParameters]:
    """Recall input parameters by query ID."""
    return INPUT_STORE.recall(query_id)

def get_latest_input(name: str) -> Optional[InputParameters]:
    """Get most recent input for an object name."""
    return INPUT_STORE.get_latest(name)

def list_inputs() -> List[Dict[str, str]]:
    """List all stored input queries."""
    return INPUT_STORE.list_queries()
