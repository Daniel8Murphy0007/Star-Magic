#!/usr/bin/env python3
"""
InputData.py - Input Dataset Management
========================================

Manages input datasets for CondensedPhysics.py calculations.

DATA FLOW:
    APIFetch.py → bodies_YYYYMMDD_HHMMSS.csv → InputData.py (this file)
                                                      ↓
                                         CondensedPhysics.solve(params)

PURPOSE:
    - Store fetched parameters from APIFetch.py queries
    - Validate parameter completeness and units
    - Provide parameter recall for repeat calculations
    - Track query history and timestamps

ARCHITECTURE:
    - Stateless storage (JSON/CSV files, no in-memory database)
    - One file per query session (bodies_YYYYMMDD_HHMMSS.csv)
    - Parameters stored in SI units only

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Created: February 11, 2026
"""

import json
import csv
import os
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field


# ═══════════════════════════════════════════════════════════════════════════════
# PARAMETER SCHEMA (SI UNITS ONLY)
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AstrophysicalParameters:
    """
    Standard parameter set for astrophysical objects.
    
    ALL VALUES IN SI UNITS:
        Mass: kg
        Distance: m
        Radius: m
        Temperature: K
        Luminosity: W
        Magnetic field: T
        Velocity: m/s
    """
    name: str
    source: str  # API source (SIMBAD, NED, Grok, etc.)
    timestamp: str
    
    # Core parameters
    mass: Optional[float] = None                    # kg
    distance: Optional[float] = None                # m
    radius: Optional[float] = None                  # m
    temperature: Optional[float] = None             # K
    luminosity: Optional[float] = None              # W
    
    # Electromagnetic
    magnetic_field: Optional[float] = None          # T (Tesla)
    charge: Optional[float] = None                  # C (Coulombs)
    
    # Kinematics
    velocity: Optional[float] = None                # m/s
    velocity_dispersion: Optional[float] = None     # m/s
    redshift: Optional[float] = None                # dimensionless
    
    # Composition
    metallicity: Optional[float] = None             # [Fe/H] or Z/Z_sun
    hydrogen_fraction: Optional[float] = None       # Mass fraction
    
    # Star formation
    star_formation_rate: Optional[float] = None     # M_sun/yr
    
    # Rotation
    period: Optional[float] = None                  # s
    angular_velocity: Optional[float] = None        # rad/s
    
    # Object type metadata
    object_type: Optional[str] = None              # 'galaxy', 'star', 'magnetar', 'black_hole'
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON/CSV storage."""
        return {
            'name': self.name,
            'source': self.source,
            'timestamp': self.timestamp,
            'mass': self.mass,
            'distance': self.distance,
            'radius': self.radius,
            'temperature': self.temperature,
            'luminosity': self.luminosity,
            'magnetic_field': self.magnetic_field,
            'charge': self.charge,
            'velocity': self.velocity,
            'velocity_dispersion': self.velocity_dispersion,
            'redshift': self.redshift,
            'metallicity': self.metallicity,
            'hydrogen_fraction': self.hydrogen_fraction,
            'star_formation_rate': self.star_formation_rate,
            'period': self.period,
            'angular_velocity': self.angular_velocity,
            'object_type': self.object_type
        }


# ═══════════════════════════════════════════════════════════════════════════════
# INPUT DATA MANAGER
# ═══════════════════════════════════════════════════════════════════════════════

class InputDataManager:
    """
    Manages input datasets for physics calculations.
    
    Responsibilities:
        1. Store fetched parameters from APIFetch.py
        2. Validate parameter completeness
        3. Convert to SI units if needed
        4. Provide parameter recall by query history
    """
    
    def __init__(self, data_dir: str = "."):
        """
        Initialize input data manager.
        
        Args:
            data_dir: Directory for storing bodies_*.csv files
        """
        self.data_dir = data_dir
        self.current_session_file = None
    
    def create_session(self) -> str:
        """
        Create new query session with timestamped file.
        
        Returns:
            Filename of created session file (bodies_YYYYMMDD_HHMMSS.csv)
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"bodies_{timestamp}.csv"
        filepath = os.path.join(self.data_dir, filename)
        
        # Create CSV with header
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'name', 'source', 'timestamp', 'mass', 'distance', 'radius',
                'temperature', 'luminosity', 'magnetic_field', 'charge',
                'velocity', 'velocity_dispersion', 'redshift', 'metallicity',
                'hydrogen_fraction', 'star_formation_rate', 'period',
                'angular_velocity', 'object_type'
            ])
        
        self.current_session_file = filepath
        return filename
    
    def store_parameters(self, params: AstrophysicalParameters) -> bool:
        """
        Store fetched parameters to current session file.
        
        Args:
            params: AstrophysicalParameters object with fetched data
        
        Returns:
            True if successful, False otherwise
        """
        if self.current_session_file is None:
            self.create_session()
        
        try:
            with open(self.current_session_file, 'a', newline='') as f:
                writer = csv.writer(f)
                data = params.to_dict()
                writer.writerow([
                    data['name'], data['source'], data['timestamp'],
                    data['mass'], data['distance'], data['radius'],
                    data['temperature'], data['luminosity'], data['magnetic_field'],
                    data['charge'], data['velocity'], data['velocity_dispersion'],
                    data['redshift'], data['metallicity'], data['hydrogen_fraction'],
                    data['star_formation_rate'], data['period'], data['angular_velocity'],
                    data['object_type']
                ])
            return True
        except Exception as e:
            print(f"Error storing parameters: {e}")
            return False
    
    def load_session(self, filename: str) -> List[AstrophysicalParameters]:
        """
        Load parameters from previous session file.
        
        Args:
            filename: Name of bodies_*.csv file
        
        Returns:
            List of AstrophysicalParameters objects
        """
        filepath = os.path.join(self.data_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Session file not found: {filename}")
            return []
        
        params_list = []
        
        try:
            with open(filepath, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    params = AstrophysicalParameters(
                        name=row['name'],
                        source=row['source'],
                        timestamp=row['timestamp'],
                        mass=float(row['mass']) if row['mass'] else None,
                        distance=float(row['distance']) if row['distance'] else None,
                        radius=float(row['radius']) if row['radius'] else None,
                        temperature=float(row['temperature']) if row['temperature'] else None,
                        luminosity=float(row['luminosity']) if row['luminosity'] else None,
                        magnetic_field=float(row['magnetic_field']) if row['magnetic_field'] else None,
                        charge=float(row['charge']) if row['charge'] else None,
                        velocity=float(row['velocity']) if row['velocity'] else None,
                        velocity_dispersion=float(row['velocity_dispersion']) if row['velocity_dispersion'] else None,
                        redshift=float(row['redshift']) if row['redshift'] else None,
                        metallicity=float(row['metallicity']) if row['metallicity'] else None,
                        hydrogen_fraction=float(row['hydrogen_fraction']) if row['hydrogen_fraction'] else None,
                        star_formation_rate=float(row['star_formation_rate']) if row['star_formation_rate'] else None,
                        period=float(row['period']) if row['period'] else None,
                        angular_velocity=float(row['angular_velocity']) if row['angular_velocity'] else None,
                        object_type=row['object_type'] if row['object_type'] else None
                    )
                    params_list.append(params)
        except Exception as e:
            print(f"Error loading session: {e}")
        
        return params_list
    
    def list_sessions(self) -> List[str]:
        """
        List all available session files.
        
        Returns:
            List of session filenames (bodies_YYYYMMDD_HHMMSS.csv)
        """
        files = os.listdir(self.data_dir)
        return sorted([f for f in files if f.startswith('bodies_') and f.endswith('.csv')])
    
    def validate_parameters(self, params: AstrophysicalParameters, required: List[str]) -> Tuple[bool, List[str]]:
        """
        Validate that required parameters are present.
        
        Args:
            params: AstrophysicalParameters to validate
            required: List of required parameter names (e.g., ['mass', 'radius'])
        
        Returns:
            (is_valid, missing_parameters)
        """
        missing = []
        
        for param_name in required:
            value = getattr(params, param_name, None)
            if value is None:
                missing.append(param_name)
        
        is_valid = len(missing) == 0
        return (is_valid, missing)


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("InputData.py - Input Dataset Manager")
    print("=" * 70)
    
    # Create manager
    manager = InputDataManager()
    
    # Create new session
    session_file = manager.create_session()
    print(f"Created session: {session_file}")
    
    # Store example parameters (from APIFetch.py)
    params = AstrophysicalParameters(
        name="Sagittarius A*",
        source="SIMBAD",
        timestamp=datetime.now().isoformat(),
        mass=4.1e6 * 1.989e30,      # 4.1 million solar masses
        distance=26000 * 3.086e16,  # 26 kpc
        radius=1.2e10,              # ~12 million km
        magnetic_field=1e-4,        # 0.1 mT
        object_type="black_hole"
    )
    
    manager.store_parameters(params)
    print(f"Stored parameters for: {params.name}")
    
    # Validate required parameters
    is_valid, missing = manager.validate_parameters(params, ['mass', 'distance', 'radius'])
    print(f"Validation: {is_valid}, Missing: {missing}")
    
    # List all sessions
    sessions = manager.list_sessions()
    print(f"\nAvailable sessions: {len(sessions)}")
    for s in sessions[-5:]:  # Show last 5
        print(f"  - {s}")
