"""
UQFFSystemsDatabase.py

Comprehensive database of astrophysical systems with NASA/Chandra/JWST/EHT observational
parameters for UQFF (Universal Quantum Field Framework) calculations.

Contains 40+ systems across multiple categories:
- Tidal Disruption Events (TDEs)
- Wolf-Rayet Stars
- Active Galactic Nuclei (AGN)
- X-ray Binaries and Microquasars
- Galaxy Clusters
- Intermediate-Mass Black Holes (IMBHs)
- Radio Transients
- Planetary Nebulae
- Star-Forming Regions
- Galaxies and Large-Scale Structure

Data Sources: NASA/ADS, Chandra X-ray Observatory, JWST, EHT, VLA, ALMA, Gaia DR3

Author: Daniel Murphy (Star-Magic UQFF Framework)
Created: March 2026
"""

import numpy as np
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass, field


# ================================
# Physical Constants (SI Units)
# ================================

class UQFFPhysicalConstants:
    """Physical constants for UQFF calculations."""
    
    # Fundamental constants
    c = 2.99792458e8  # Speed of light (m/s)
    G = 6.67430e-11  # Gravitational constant (m³/kg·s²)
    h = 6.62607015e-34  # Planck constant (J·s)
    hbar = 1.054571817e-34  # Reduced Planck constant (J·s)
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    
    # Astrophysical constants
    M_sun = 1.98847e30  # Solar mass (kg)
    L_sun = 3.828e26  # Solar luminosity (W)
    R_sun = 6.96e8  # Solar radius (m)
    
    # Distance units
    pc = 3.0857e16  # Parsec (m)
    kpc = 3.0857e19  # Kiloparsec (m)
    Mpc = 3.0857e22  # Megaparsec (m)
    AU = 1.496e11  # Astronomical unit (m)
    ly = 9.4607e15  # Light-year (m)
    
    # Time units
    yr = 3.15576e7  # Year (s)
    day = 86400  # Day (s)
    
    # Energy units
    eV = 1.602176634e-19  # Electron-volt (J)
    keV = 1.602176634e-16  # Kilo-electron-volt (J)
    MeV = 1.602176634e-13  # Mega-electron-volt (J)
    GeV = 1.602176634e-10  # Giga-electron-volt (J)
    
    # Magnetic field
    Gauss = 1e-4  # Gauss to Tesla conversion


# ================================
# System Parameter Dataclass
# ================================

@dataclass
class AstrophysicalSystem:
    """Dataclass for astrophysical system parameters."""
    
    name: str
    category: str
    
    # Position and distance
    ra: Optional[float] = None  # Right ascension (degrees)
    dec: Optional[float] = None  # Declination (degrees)
    distance: Optional[float] = None  # Distance (parsecs)
    redshift: Optional[float] = None  # Redshift z
    
    # Mass and size
    mass: Optional[float] = None  # Mass (solar masses)
    radius: Optional[float] = None  # Radius (solar radii or parsecs)
    
    # Luminosity and temperature
    L_bol: Optional[float] = None  # Bolometric luminosity (erg/s)
    L_X: Optional[float] = None  # X-ray luminosity (erg/s)
    T: Optional[float] = None  # Temperature (K)
    
    # Velocities
    v_exp: Optional[float] = None  # Expansion velocity (km/s)
    v_rot: Optional[float] = None  # Rotation velocity (km/s)
    v_wind: Optional[float] = None  # Wind velocity (km/s)
    
    # Magnetic field
    B: Optional[float] = None  # Magnetic field (Gauss)
    
    # Density and pressure
    n_e: Optional[float] = None  # Electron density (cm⁻³)
    n_H: Optional[float] = None  # Hydrogen density (cm⁻³)
    pressure: Optional[float] = None  # Pressure (dyn/cm²)
    
    # Star formation
    SFR: Optional[float] = None  # Star formation rate (M☉/yr)
    
    # Spectral features
    spectral_type: Optional[str] = None
    metallicity: Optional[float] = None  # [Fe/H]
    
    # Source citations
    sources: List[str] = field(default_factory=list)
    
    # Additional parameters (flexible dict for system-specific data)
    additional_params: Dict = field(default_factory=dict)


# ================================
# UQFF Systems Database
# ================================

class UQFFSystemsDatabase:
    """
    Comprehensive database of astrophysical systems for UQFF calculations.
    
    Contains observational parameters from NASA, Chandra, JWST, EHT, and other
    major astronomical facilities.
    """
    
    def __init__(self):
        """Initialize the systems database."""
        self.systems = self._initialize_database()
        self.categories = self._get_categories()
    
    def _initialize_database(self) -> Dict[str, AstrophysicalSystem]:
        """Initialize complete systems database."""
        
        systems = {}
        
        # ================================
        # TIDAL DISRUPTION EVENTS (TDEs)
        # ================================
        
        systems['AT2024tvd'] = AstrophysicalSystem(
            name='AT2024tvd',
            category='TDE',
            ra=197.5284,
            dec=-12.3471,
            distance=1.2e8,  # 120 Mpc
            redshift=0.0273,
            mass=1e6,  # SMBH mass (M☉)
            L_bol=1e44,  # Peak luminosity (erg/s)
            L_X=5e43,  # X-ray luminosity (erg/s)
            T=3e4,  # Blackbody temperature (K)
            v_exp=1e4,  # Outflow velocity (km/s)
            sources=[
                'Pasham et al. 2024, Nature',
                'NASA/Swift J0513.4-6547',
                'Chandra DDT Observation'
            ],
            additional_params={
                'peak_time': 15.2,  # days
                'decay_timescale': 33.0,  # days
                'recurrence_time': 21.0,  # days
                'disk_temp_peak': 8e4  # K
            }
        )
        
        systems['AT2019qiz'] = AstrophysicalSystem(
            name='AT2019qiz',
            category='TDE',
            ra=45.3719,
            dec=-32.9485,
            distance=2.67e8,  # 267 Mpc
            redshift=0.0151,
            mass=1e6,  # SMBH mass (M☉)
            L_bol=2.5e44,  # Peak luminosity (erg/s)
            L_X=1e43,
            T=4e4,
            v_exp=8e3,
            sources=[
                'Nicholl et al. 2020, MNRAS',
                'Swift/UVOT + XRT',
                'Las Cumbres Observatory'
            ],
            additional_params={
                'peak_time': 49.0,
                'rise_time': 32.0,
                'BH_spin': 0.4  # dimensionless
            }
        )
        
        # ================================
        # WOLF-RAYET STARS
        # ================================
        
        systems['WR124'] = AstrophysicalSystem(
            name='WR 124',
            category='Wolf-Rayet',
            ra=267.8958,
            dec=-41.3431,
            distance=5000,  # pc
            mass=20,  # Current mass (M☉)
            radius=5.5,  # R☉
            L_bol=2.5e39,  # erg/s
            T=5e4,  # K
            v_wind=2000,  # km/s
            B=100,  # Gauss
            spectral_type='WN8h',
            sources=[
                'JWST NIRCam 2022',
                'van der Hucht 2001, NewAR',
                'Chandra ACIS Observation'
            ],
            additional_params={
                'mass_loss_rate': 3e-5,  # M☉/yr
                'wind_momentum': 1e29,  # g·cm/s²
                'nebula_age': 1e4,  # years
                'nebula_mass': 3.0  # M☉
            }
        )
        
        systems['WR140'] = AstrophysicalSystem(
            name='WR 140',
            category='Wolf-Rayet',
            ra=300.9875,
            dec=43.8503,
            distance=1670,  # pc
            mass=15,  # WC7 component (M☉)
            L_bol=3e39,
            T=7e4,
            v_wind=2800,
            B=150,
            spectral_type='WC7pd + O5.5fc',
            sources=[
                'JWST MIRI 2022',
                'Williams et al. 2022, Nature',
                'Chandra HETG Observation'
            ],
            additional_params={
                'orbital_period': 2896,  # days (7.93 years)
                'eccentricity': 0.896,
                'dust_production': True,
                'pinwheel_structure': True
            }
        )
        
        # ================================
        # ACTIVE GALACTIC NUCLEI (AGN)
        # ================================
        
        systems['CentaurusA'] = AstrophysicalSystem(
            name='Centaurus A (NGC 5128)',
            category='AGN',
            ra=201.365,
            dec=-43.019,
            distance=3.8e6,  # 3.8 Mpc
            redshift=0.00183,
            mass=5.5e7,  # SMBH mass (M☉)
            L_bol=1e43,
            L_X=2e42,
            v_rot=150,  # km/s
            B=1e-4,  # Gauss (radio lobes)
            sources=[
                'EHT Collaboration 2021',
                'Chandra 1 Ms Deep Field',
                'JWST NIRCam 2023'
            ],
            additional_params={
                'jet_length': 6e5,  # pc
                'jet_power': 1e43,  # erg/s
                'inner_jet_angle': 50,  # degrees
                'dust_mass': 1e7  # M☉
            }
        )
        
        systems['M87'] = AstrophysicalSystem(
            name='M87 (Virgo A)',
            category='AGN',
            ra=187.7059,
            dec=12.3911,
            distance=1.68e7,  # 16.8 Mpc
            redshift=0.00436,
            mass=6.5e9,  # SMBH mass from EHT (M☉)
            L_bol=1e42,
            L_X=2e41,
            v_rot=550,
            B=1e-3,  # Gauss (near horizon)
            sources=[
                'EHT Collaboration 2019, ApJL',
                'Chandra HETG Observations',
                'JWST MIRI 2023'
            ],
            additional_params={
                'schwarzschild_radius': 1.95e13,  # cm
                'jet_speed': 0.99,  # c
                'accretion_rate': 9e-3,  # M☉/yr
                'BH_spin': 0.9  # dimensionless
            }
        )
        
        systems['SgrA_star'] = AstrophysicalSystem(
            name='Sagittarius A*',
            category='AGN',
            ra=266.4168,
            dec=-29.0078,
            distance=8178,  # pc (Galactic Center)
            mass=4.15e6,  # SMBH mass (M☉)
            L_bol=1e36,  # Very low luminosity
            L_X=2e33,
            B=30,  # Gauss (near horizon)
            sources=[
                'EHT Collaboration 2022, ApJL',
                'GRAVITY Collaboration 2022',
                'Chandra X-ray Observatory'
            ],
            additional_params={
                'schwarzschild_radius': 1.22e12,  # cm
                'flare_luminosity': 1e35,  # erg/s
                'accretion_rate': 1e-7,  # M☉/yr
                'orbital_period_S2': 16.05  # years
            }
        )
        
        # ================================
        # X-RAY BINARIES & MICROQUASARS
        # ================================
        
        systems['SS433'] = AstrophysicalSystem(
            name='SS 433',
            category='X-ray Binary',
            ra=287.9565,
            dec=4.9829,
            distance=5500,  # pc
            mass=10,  # BH mass (M☉)
            L_X=1e39,
            v_exp=26000,  # Relativistic jets (0.26c)
            B=1e8,  # Gauss (near compact object)
            sources=[
                'Chandra ACIS-S Observation',
                'VLA Radio Maps',
                'Gaia DR3'
            ],
            additional_params={
                'jet_precession_period': 162.5,  # days
                'orbital_period': 13.08,  # days
                'jet_opening_angle': 20,  # degrees
                'disk_inclination': 79  # degrees
            }
        )
        
        systems['CygnusX1'] = AstrophysicalSystem(
            name='Cygnus X-1',
            category='X-ray Binary',
            ra=299.5903,
            dec=35.2016,
            distance=2200,  # pc
            mass=21.2,  # BH mass (M☉)
            L_X=2e37,
            B=1e9,
            sources=[
                'Chandra HETG 2010',
                'Miller-Jones et al. 2021, Science',
                'NICER Observations'
            ],
            additional_params={
                'BH_spin': 0.998,  # Near-maximal
                'orbital_period': 5.6,  # days
                'companion_mass': 40.6,  # M☉ (O9.7Iab)
                'accretion_rate': 1e-8  # M☉/yr
            }
        )
        
        # ================================
        # GALAXY CLUSTERS
        # ================================
        
        systems['PerseusCluster'] = AstrophysicalSystem(
            name='Perseus Cluster (Abell 426)',
            category='Galaxy Cluster',
            ra=49.9508,
            dec=41.5117,
            distance=7.2e7,  # 72 Mpc
            redshift=0.0179,
            mass=6.5e14,  # Cluster mass (M☉)
            L_X=6e44,
            T=6e7,  # ICM temperature (K)
            n_e=0.04,  # cm⁻³ (central)
            B=25,  # µG (ICM magnetic field)
            sources=[
                'Chandra 1.4 Ms Deep Field',
                'Hitomi X-ray Observatory',
                'JWST NIRCam 2022'
            ],
            additional_params={
                'velocity_dispersion': 1200,  # km/s
                'cooling_flow_rate': 200,  # M☉/yr (classical value)
                'AGN_feedback_power': 1e44,  # erg/s
                'virial_radius': 2.5e6  # pc
            }
        )
        
        systems['PLCKG287'] = AstrophysicalSystem(
            name='PLCK G287.0+32.9',
            category='Galaxy Cluster',
            ra=171.25,
            dec=-28.67,
            distance=1.9e9,  # 1.9 Gpc
            redshift=0.39,
            mass=1.3e15,  # M☉
            L_X=2e45,
            T=1.2e8,  # K
            n_e=0.01,
            sources=[
                'Planck SZ Survey',
                'XMM-Newton Observation',
                'Chandra Follow-up'
            ],
            additional_params={
                'SZ_signal': 15.2,  # Compton-y parameter
                'merger_state': 'disturbed',
                'shock_velocity': 2500  # km/s
            }
        )
        
        systems['PSZ2G181'] = AstrophysicalSystem(
            name='PSZ2 G181.06+48.47',
            category='Galaxy Cluster',
            ra=132.85,
            dec=19.42,
            distance=2.1e9,  # 2.1 Gpc
            redshift=0.42,
            mass=8e14,
            T=9e7,
            sources=[
                'Planck SZ Survey',
                'SDSS Optical Confirmation'
            ]
        )
        
        # ================================
        # INTERMEDIATE-MASS BLACK HOLES
        # ================================
        
        systems['J1610'] = AstrophysicalSystem(
            name='J1610+1811',
            category='IMBH',
            ra=242.5,
            dec=18.19,
            distance=3.5e8,  # 350 Mpc
            redshift=0.078,
            mass=5e4,  # IMBH mass (M☉)
            L_X=8e41,
            sources=[
                'Chandra Source Catalog',
                'Lin et al. 2018, Nature Astronomy',
                'HST Follow-up'
            ],
            additional_params={
                'host_type': 'dwarf galaxy',
                'X_ray_variability': 0.3,
                'Eddington_ratio': 0.01
            }
        )
        
        systems['HLX1'] = AstrophysicalSystem(
            name='HLX-1 (ESO 243-49)',
            category='IMBH',
            ra=16.1833,
            dec=-49.2592,
            distance=9.5e6,  # 9.5 Mpc
            mass=2e4,  # IMBH mass (M☉)
            L_X=1e42,  # Peak
            sources=[
                'Chandra ACIS-S',
                'Farrell et al. 2009, Nature',
                'Swift XRT Monitoring'
            ],
            additional_params={
                'orbital_period': 380,  # days (candidate)
                'recurrence_time': 367,  # days
                'outburst_duration': 10  # days
            }
        )
        
        # ================================
        # RADIO TRANSIENTS
        # ================================
        
        systems['ASKAPJ1832'] = AstrophysicalSystem(
            name='ASKAP J1832-0911',
            category='Radio Transient',
            ra=278.031,
            dec=-9.185,
            distance=3000,  # pc (estimated)
            L_bol=1e32,  # Radio luminosity equivalent (erg/s)
            sources=[
                'Hurley-Walker et al. 2022, Nature',
                'ASKAP Variables and Slow Transients Survey',
                'MWA Follow-up'
            ],
            additional_params={
                'period': 1318,  # seconds
                'duty_cycle': 0.1,
                'magnetar_candidate': True,
                'radio_frequency': 154  # MHz
            }
        )
        
        systems['FRB20121102A'] = AstrophysicalSystem(
            name='FRB 20121102A',
            category='Radio Transient',
            ra=82.9975,
            dec=33.1475,
            distance=9.7e8,  # 970 Mpc
            redshift=0.193,
            L_bol=1e40,  # Peak burst luminosity
            B=1e14,  # Gauss (magnetar surface)
            sources=[
                'Chatterjee et al. 2017, Nature',
                'VLA Localization',
                'Gemini Optical Host ID'
            ],
            additional_params={
                'burst_rate': 0.02,  # per hour
                'rotation_measure': 1e5,  # rad/m²
                'burst_width': 3  # milliseconds
            }
        )
        
        # ================================
        # PLANETARY NEBULAE
        # ================================
        
        systems['HelixNebula'] = AstrophysicalSystem(
            name='Helix Nebula (NGC 7293)',
            category='Planetary Nebula',
            ra=337.4111,
            dec=-20.8378,
            distance=213,  # pc
            mass=0.6,  # Central star mass (M☉)
            radius=2.5,  # Nebula radius (pc)
            L_bol=2e34,  # erg/s
            T=1.2e5,  # Central star temperature (K)
            v_exp=20,  # km/s
            n_e=500,  # cm⁻³
            sources=[
                'JWST NIRCam/MIRI 2022',
                'Chandra ACIS-S',
                'HST WFPC2'
            ],
            additional_params={
                'age': 1e4,  # years
                'nebula_mass': 0.3,  # M☉
                'ionization_structure': 'clumpy',
                'knot_count': 3500
            }
        )
        
        systems['RAquarii'] = AstrophysicalSystem(
            name='R Aquarii',
            category='Symbiotic Binary',
            ra=355.5208,
            dec=-15.2839,
            distance=218,  # pc
            mass=0.8,  # White dwarf mass (M☉)
            L_bol=5e34,
            T=1e5,  # WD temperature
            v_exp=50,  # Jet velocity (km/s)
            sources=[
                'Chandra ACIS-S',
                'JWST MIRI 2023',
                'HST/STIS Long-term Monitoring'
            ],
            additional_params={
                'orbital_period': 44,  # years
                'jet_length': 0.15,  # pc
                'nova_eruptions': True,
                'last_outburst': 2000
            }
        )
        
        # ================================
        # STAR-FORMING REGIONS
        # ================================
        
        systems['TarantulaNebula'] = AstrophysicalSystem(
            name='Tarantula Nebula (30 Doradus)',
            category='Star-Forming Region',
            ra=84.68,
            dec=-69.1,
            distance=4.9e4,  # 49 kpc (LMC)
            mass=4.5e5,  # Total mass (M☉)
            radius=200,  # pc
            L_bol=1e40,  # erg/s
            T=1e4,  # Ionized gas temperature (K)
            n_H=100,  # cm⁻³
            SFR=0.1,  # M☉/yr
            sources=[
                'JWST NIRCam 2022',
                'Chandra ACIS-I Deep Survey',
                'ALMA 850 µm Continuum'
            ],
            additional_params={
                'age': 2e6,  # years
                'stellar_cluster': 'R136',
                'most_massive_star': 265,  # M☉ (R136a1)
                'OB_star_count': 1000
            }
        )
        
        systems['CarinaNebula'] = AstrophysicalSystem(
            name='Carina Nebula (NGC 3372)',
            category='Star-Forming Region',
            ra=161.265,
            dec=-59.875,
            distance=2300,  # pc
            mass=2e5,  # M☉
            radius=100,  # pc
            L_bol=5e39,
            T=1e4,
            n_H=300,
            SFR=0.01,
            sources=[
                'JWST NIRCam 2022',
                'Chandra Carina Complex Project',
                'HST Treasury Survey'
            ],
            additional_params={
                'age': 3e6,
                'eta_carinae_mass': 100,  # M☉
                'pillar_structures': True,
                'young_stellar_objects': 14000
            }
        )
        
        systems['OrionNebula'] = AstrophysicalSystem(
            name='Orion Nebula (M42)',
            category='Star-Forming Region',
            ra=83.8221,
            dec=-5.3911,
            distance=414,  # pc
            mass=2000,  # M☉
            radius=12,  # pc
            L_bol=1e37,
            T=1e4,
            n_H=1e4,  # High density HII region
            SFR=0.001,
            sources=[
                'JWST NIRCam 2023',
                'Chandra Orion Ultradeep Project',
                'ALMA Band 6'
            ],
            additional_params={
                'age': 3e5,  # years
                'Trapezium_cluster': True,
                'proplyds_count': 180,
                'Theta1C_mass': 40  # M☉
            }
        )
        
        # ================================
        # GALAXIES
        # ================================
        
        systems['M31'] = AstrophysicalSystem(
            name='M31 (Andromeda)',
            category='Galaxy',
            ra=10.6847,
            dec=41.2692,
            distance=7.8e5,  # 780 kpc
            mass=1.5e12,  # Total mass including dark matter (M☉)
            L_bol=3.6e43,  # erg/s
            v_rot=250,  # km/s
            SFR=0.4,  # M☉/yr
            metallicity=-0.1,  # [Fe/H]
            sources=[
                'JWST NIRCam Deep Field 2023',
                'Gaia DR3 M31 Survey',
                'Chandra M31 Survey'
            ],
            additional_params={
                'SMBH_mass': 1.4e8,  # M☉
                'bulge_mass': 2e10,  # M☉
                'disk_scale_length': 5.4e3,  # pc
                'M31_M33_merger_time': 4.5e9  # years
            }
        )
        
        systems['M87'] = systems.get('M87')  # Already defined under AGN
        
        systems['NGC1365'] = AstrophysicalSystem(
            name='NGC 1365',
            category='Galaxy',
            ra=53.4017,
            dec=-36.1403,
            distance=1.77e7,  # 17.7 Mpc
            redshift=0.00546,
            mass=2e11,  # M☉
            L_bol=2e43,
            L_X=1e41,
            v_rot=280,
            SFR=5,  # M☉/yr
            sources=[
                'JWST MIRI Barred Galaxy Survey',
                'Chandra Deep Survey',
                'ALMA 1.3mm Continuum'
            ],
            additional_params={
                'SMBH_mass': 2e6,  # M☉
                'bar_length': 5e3,  # pc
                'AGN_type': 'Seyfert 1.8',
                'nuclear_SFR': 0.5  # M☉/yr
            }
        )
        
        systems['NGC3596'] = AstrophysicalSystem(
            name='NGC 3596',
            category='Galaxy',
            ra=168.7958,
            dec=14.8522,
            distance=2.3e7,  # 23 Mpc
            redshift=0.00556,
            mass=5e10,  # M☉
            L_bol=5e42,
            v_rot=180,
            SFR=1.2,
            metallicity=-0.2,
            sources=[
                'SDSS DR17',
                'GALEX UV Survey',
                '2MASS Extended Source Catalog'
            ],
            additional_params={
                'morphology': 'SABc',
                'inclination': 65,  # degrees
                'HI_mass': 1e9  # M☉
            }
        )
        
        # ================================
        # SUPERNOVA REMNANTS
        # ================================
        
        systems['CrabNebula'] = AstrophysicalSystem(
            name='Crab Nebula (M1)',
            category='Supernova Remnant',
            ra=83.6333,
            dec=22.0145,
            distance=2000,  # pc
            mass=4.6,  # Ejecta mass (M☉)
            radius=3.5,  # pc
            L_bol=1e38,  # erg/s (including pulsar)
            T=1e4,  # Filament temperature (K)
            v_exp=1500,  # km/s
            B=3e-4,  # Gauss (average)
            sources=[
                'Chandra ACIS-S Deep Observations',
                'JWST NIRCam/MIRI 2023',
                'NICER Pulsar Timing'
            ],
            additional_params={
                'pulsar_period': 0.033,  # seconds
                'pulsar_spindown': 4.2e-13,  # s/s
                'age': 968,  # years (SN 1054)
                'pulsar_Edot': 5e38  # erg/s
            }
        )
        
        systems['CasA'] = AstrophysicalSystem(
            name='Cassiopeia A',
            category='Supernova Remnant',
            ra=350.85,
            dec=58.815,
            distance=3400,  # pc
            mass=2,  # Ejecta mass (M☉)
            radius=2.5,  # pc
            v_exp=5000,  # km/s
            B=5e-4,  # Gauss
            T=3e7,  # X-ray emitting shocked gas (K)
            sources=[
                'Chandra 1 Ms Deep Field',
                'JWST NIRCam 2023 (shocked ejecta)',
                'VLA 3 GHz Survey'
            ],
            additional_params={
                'age': 350,  # years (SN ~1680)
                'NS_candidate': True,
                'Ti44_detection': True,
                'jet_structure': True
            }
        )
        
        systems['VelaJr'] = AstrophysicalSystem(
            name='Vela Jr. (RX J0852.0-4622)',
            category='Supernova Remnant',
            ra=133.0,
            dec=-46.37,
            distance=750,  # pc
            radius=25,  # pc
            v_exp=1000,
            T=5e6,
            sources=[
                'Chandra ACIS-I Mosaic',
                'H.E.S.S. Gamma-ray Detection',
                'XMM-Newton Survey'
            ],
            additional_params={
                'age': 700,  # years
                'TeV_gamma_ray': True,
                'shock_velocity': 2000  # km/s
            }
        )
        
        # ================================
        # MOLECULAR CLOUDS
        # ================================
        
        systems['SgrB2'] = AstrophysicalSystem(
            name='Sagittarius B2',
            category='Molecular Cloud',
            ra=266.835,
            dec=-28.384,
            distance=8500,  # pc (Galactic Center)
            mass=3e6,  # M☉
            radius=15,  # pc
            T=100,  # Dust temperature (K)
            n_H=1e6,  # cm⁻³ (dense cores)
            SFR=0.01,
            sources=[
                'ALMA Band 3-10 Survey',
                'Chandra X-ray Point Sources',
                'JWST NIRCam 2023'
            ],
            additional_params={
                'molecular_species': 80,
                'HII_regions': 3,
                'masers': True,
                'complex_organics': True
            }
        )
        
        systems['W51'] = AstrophysicalSystem(
            name='W51',
            category='Molecular Cloud',
            ra=290.93,
            dec=14.52,
            distance=5400,  # pc
            mass=1e6,  # M☉
            radius=30,  # pc
            T=80,
            n_H=1e5,
            SFR=0.1,
            sources=[
                'ALMA High-resolution Survey',
                'Chandra ACIS-I',
                'VLA 6 GHz Continuum'
            ],
            additional_params={
                'UCHII_regions': 10,
                'massive_YSOs': 50,
                'outflow_activity': True
            }
        )
        
        # ================================
        # GLOBULAR CLUSTERS
        # ================================
        
        systems['OmegaCen'] = AstrophysicalSystem(
            name='Omega Centauri (NGC 5139)',
            category='Globular Cluster',
            ra=201.697,
            dec=-47.479,
            distance=5200,  # pc
            mass=4e6,  # M☉
            radius=86,  # pc (tidal radius)
            L_bol=1.2e39,  # erg/s
            v_rot=8,  # km/s (internal)
            metallicity=-1.5,  # [Fe/H] (average)
            sources=[
                'HST ACS Survey',
                'Gaia DR3 Proper Motions',
                'Chandra X-ray Binary Survey'
            ],
            additional_params={
                'age': 1.2e10,  # years
                'stellar_population': 1e7,
                'IMBH_candidate': True,
                'IMBH_mass': 4e4  # M☉ (if present)
            }
        )
        
        systems['M15'] = AstrophysicalSystem(
            name='M15 (NGC 7078)',
            category='Globular Cluster',
            ra=322.493,
            dec=12.168,
            distance=1.04e4,  # pc
            mass=5e5,  # M☉
            radius=88,  # pc
            metallicity=-2.3,  # [Fe/H]
            sources=[
                'HST WFPC2',
                'Chandra Core Survey',
                'Gaia DR3'
            ],
            additional_params={
                'age': 1.3e10,
                'core_collapsed': True,
                'X_ray_binaries': 8,
                'pulsar_count': 9
            }
        )
        
        # ================================
        # QUASARS
        # ================================
        
        systems['3C273'] = AstrophysicalSystem(
            name='3C 273',
            category='Quasar',
            ra=187.2779,
            dec=2.0524,
            distance=7.49e8,  # 749 Mpc
            redshift=0.158,
            mass=8.86e8,  # SMBH mass (M☉)
            L_bol=1e47,  # erg/s
            L_X=1e45,
            v_exp=3000,  # Broad line region (km/s)
            sources=[
                'Chandra HETG 2000-2020',
                'JWST NIRSpec 2023',
                'Fermi-LAT Gamma-ray'
            ],
            additional_params={
                'jet_length': 4e5,  # pc
                'jet_speed': 0.997,  # c
                'Eddington_ratio': 0.1,
                'variability_timescale': 1  # day
            }
        )
        
        systems['SDSSJ1030'] = AstrophysicalSystem(
            name='SDSS J1030+0524',
            category='Quasar',
            ra=157.5,
            dec=5.4,
            distance=1.8e10,  # 18 Gpc
            redshift=6.31,  # High-z quasar
            mass=2e9,  # M☉
            L_bol=1e47,
            sources=[
                'JWST NIRCam High-z Survey',
                'Keck LRIS Spectroscopy',
                'Gemini GNIRS'
            ],
            additional_params={
                'cosmic_age': 8.7e8,  # years (when observed)
                'BH_growth_timescale': 4e8,  # years
                'most_distant_at_discovery': True
            }
        )
        
        # ================================
        # NEUTRON STARS & MAGNETARS
        # ================================
        
        systems['SGR1745'] = AstrophysicalSystem(
            name='SGR 1745-2900',
            category='Magnetar',
            ra=266.415,
            dec=-29.007,
            distance=8000,  # pc (near Sgr A*)
            mass=1.4,  # NS mass (M☉)
            radius=1.2e-5,  # 12 km in R☉ units
            L_X=5e34,  # Quiescent (erg/s)
            T=5e6,  # Surface temperature (K)
            B=2e14,  # Gauss (surface dipole)
            sources=[
                'Chandra ACIS-S Monitoring',
                'NuSTAR Hard X-ray',
                'Swift BAT Burst Detection'
            ],
            additional_params={
                'spin_period': 3.76,  # seconds
                'outburst_date': '2013-04-24',
                'burst_luminosity': 1e41,  # erg/s (peak)
                'distance_to_SgrA': 2.4  # arcsec
            }
        )
        
        systems['J1745'] = systems['SGR1745']  # Alias
        
        systems['PSR_J0737'] = AstrophysicalSystem(
            name='PSR J0737-3039A/B',
            category='Binary Pulsar',
            ra=114.425,
            dec=-30.655,
            distance=1150,  # pc
            mass=1.3381,  # Pulsar A mass (M☉)
            L_bol=1e30,  # Spin-down luminosity (erg/s)
            sources=[
                'Parkes Radio Telescope',
                'GBT Timing Campaign',
                'Kramer et al. 2021, Physical Review X'
            ],
            additional_params={
                'pulsar_A_period': 0.02270,  # seconds
                'pulsar_B_period': 2.7734,  # seconds
                'orbital_period': 0.102251563,  # days (2.45 hours)
                'periastron_advance': 16.9,  # deg/year
                'GR_test_precision': 0.0002  # fractional
            }
        )

        # ================================
        # SESSIONS 58-60 SYSTEMS (grok_share_8d951e12 — March 2026)
        # ================================

        systems['SGR0501'] = AstrophysicalSystem(
            name='SGR 0501+4516',
            category='Magnetar',
            ra=75.278,
            dec=45.277,
            distance=2000,  # pc
            mass=1.4,  # M☉ (canonical neutron star)
            L_bol=1e36,  # erg/s (burst luminosity)
            B=8e13,  # T (surface magnetic field ~8e13 T)
            sources=[
                'Swift/BAT Discovery 2008',
                'Rea et al. 2009, ApJ',
                'UQFF Session 58 PAPER_226'
            ],
            additional_params={
                'spin_period': 5.76,  # seconds
                'period_derivative': 1.5e-11,  # s/s
                'burst_fluence': 1e-7,  # erg/cm²
                'GW_back_reaction': True,  # 11-term MUGE includes GW back-reaction
                'n_muge_terms': 11
            }
        )

        systems['TapestryLMC'] = AstrophysicalSystem(
            name='LMC N49 / Tapestry Star Formation',
            category='Supernova Remnant',
            ra=81.316,
            dec=-66.083,
            distance=50000,  # pc (LMC distance ~50 kpc)
            mass=20.0,  # M☉ (progenitor estimate)
            L_bol=2e38,  # erg/s
            v_wind=2000.0,  # km/s (stellar wind velocity)
            sources=[
                'Chandra X-ray Observatory',
                'Park et al. 2012, ApJ',
                'UQFF Session 58 PAPER_227'
            ],
            additional_params={
                'remnant_age': 4800,  # years
                'shock_velocity': 800,  # km/s
                'X_ray_luminosity': 2e36  # erg/s
            }
        )

        systems['Westerlund2'] = AstrophysicalSystem(
            name='Westerlund 2',
            category='Star-Forming Region',
            ra=159.775,
            dec=-57.764,
            distance=8000,  # pc (8 kpc)
            mass=1e4,  # M☉ (cluster mass)
            L_bol=1e41,  # erg/s (total OB stars)
            v_wind=2500.0,  # km/s (OB star winds)
            SFR=0.1,  # M☉/yr
            sources=[
                'HST/WFC3 Imaging',
                'Zeidler et al. 2015, A&A',
                'UQFF Session 58 PAPER_228'
            ],
            additional_params={
                'n_OB_stars': 300,
                'cluster_age': 2e6,  # years (~2 Myr)
                'wind_power': 5e38  # erg/s
            }
        )

        systems['PillarsCreation'] = AstrophysicalSystem(
            name='Pillars of Creation (M16 Eagle Nebula)',
            category='Star-Forming Region',
            ra=274.7,
            dec=-13.807,
            distance=2000,  # pc (2 kpc)
            mass=1e3,  # M☉ (pillar gas mass)
            L_bol=1e38,  # erg/s
            sources=[
                'HST/WFPC2 1995 (Hester et al.)',
                'JWST NIRCam 2022',
                'UQFF Session 58 PAPER_229'
            ],
            additional_params={
                'pillar_length': 4,  # light-years
                'photoevaporation_rate': 70,  # M☉/Myr
                'E_erosion_flag': True  # E(t) erosion model
            }
        )

        systems['NGC2525'] = AstrophysicalSystem(
            name='NGC 2525 + SN2018gv',
            category='Galaxy',
            ra=120.625,
            dec=-11.425,
            distance=2.2e7,  # pc (22 Mpc)
            mass=1e10,  # M☉ (galaxy stellar mass)
            L_bol=2e43,  # erg/s
            SFR=1.5,  # M☉/yr
            sources=[
                'HST Type Ia Cepheid programme',
                'Yang et al. 2021, ApJ',
                'UQFF Session 58 PAPER_230'
            ],
            additional_params={
                'SN2018gv_peak_mag': -19.4,  # Ia peak absolute magnitude
                'SN2018gv_date': '2018-01-15',
                'negative_g_SN': True  # NGC2525 SN negative-g model
            }
        )

        systems['HUDF_highz'] = AstrophysicalSystem(
            name='Hubble Ultra Deep Field (z=3.5)',
            category='Cosmological Field',
            ra=53.163,
            dec=-27.791,
            distance=0,  # cosmological — use redshift
            redshift=3.5,
            mass=1e11,  # M☉ (typical galaxy at z~3.5)
            L_bol=1e44,  # erg/s
            SFR=50.0,  # M☉/yr (high-z starburst)
            sources=[
                'HST ACS + WFC3 HUDF Campaign',
                'Beckwith et al. 2006, AJ',
                'UQFF Session 58 PAPER_231'
            ],
            additional_params={
                'lookback_time': 1.19e10,  # years (~11.9 Gyr)
                'n_galaxies_per_arcmin2': 10000,
                'cosmic_I_model': True  # I(t) double-modulation model
            }
        )

        systems['NGC1792'] = AstrophysicalSystem(
            name='NGC 1792',
            category='Galaxy',
            ra=76.311,
            dec=-37.981,
            distance=1.3e7,  # pc (13 Mpc)
            mass=4e10,  # M☉
            L_bol=1e44,  # erg/s
            SFR=5.0,  # M☉/yr (starburst)
            sources=[
                'Calzetti et al. 2015, ApJ',
                'GALEX + Spitzer SINGS',
                'UQFF Session 58 PAPER_232'
            ],
            additional_params={
                'nuclear_starburst': True,
                'H_alpha_luminosity': 3e40,  # erg/s
                'stellar_forge_model': True  # NGC1792 Stellar Forge F_SFR model
            }
        )

        systems['Antennae'] = AstrophysicalSystem(
            name='Antennae Galaxies (NGC 4038/4039)',
            category='Galaxy Merger',
            ra=180.474,
            dec=-18.867,
            distance=2.2e7,  # pc (22 Mpc)
            mass=2e11,  # M☉ (combined)
            L_bol=2e44,  # erg/s
            SFR=20.0,  # M☉/yr (merger-driven)
            sources=[
                'Whitmore et al. 1999, AJ',
                'Chandra X-ray Array',
                'UQFF Session 58 PAPER_235'
            ],
            additional_params={
                'merger_stage': 'ongoing',
                'n_super_star_clusters': 1000,
                'double_intensity_model': True  # double-I(t) merger model
            }
        )

        systems['RingsOfRelativity'] = AstrophysicalSystem(
            name='Rings of Relativity (Einstein Lensing Arc)',
            category='Gravitational Lens',
            ra=0.0,  # generic — depends on line of sight
            dec=0.0,
            distance=5e8,  # pc (500 Mpc typical lens)
            mass=1e13,  # M☉ (lens galaxy cluster mass)
            L_bol=1e45,  # erg/s (background source)
            sources=[
                'ESA Hubble Strong Lensing Survey',
                'UQFF Session 60 PAPER_242'
            ],
            additional_params={
                'einstein_ring_model': True,  # L_t = (GM/c²r)·(D_LS/D_S) static
                'lensing_correction': 'L_t=(GM/c²r)*(D_LS/D_S)',
                'two_mode_oscillation': True
            }
        )
        
        return systems
    
    def _get_categories(self) -> List[str]:
        """Get unique category list."""
        return sorted(list(set(system.category for system in self.systems.values())))
    
    # ================================
    # DATA ACCESS METHODS
    # ================================
    
    def get_system(self, name: str) -> Optional[AstrophysicalSystem]:
        """
        Get system by name (case-insensitive).
        
        Args:
            name: System name
            
        Returns:
            AstrophysicalSystem object or None if not found
        """
        # Try exact match first
        if name in self.systems:
            return self.systems[name]
        
        # Try case-insensitive search
        name_lower = name.lower()
        for key, system in self.systems.items():
            if key.lower() == name_lower or system.name.lower() == name_lower:
                return system
        
        return None
    
    def get_systems_by_category(self, category: str) -> List[AstrophysicalSystem]:
        """
        Get all systems in a category.
        
        Args:
            category: Category name (case-insensitive)
            
        Returns:
            List of AstrophysicalSystem objects
        """
        category_lower = category.lower()
        return [
            system for system in self.systems.values()
            if system.category.lower() == category_lower
        ]
    
    def list_systems(self) -> List[str]:
        """Get list of all system names."""
        return sorted(list(self.systems.keys()))
    
    def list_categories(self) -> List[str]:
        """Get list of all categories."""
        return self.categories
    
    def get_system_summary(self, name: str) -> str:
        """
        Get formatted summary of system parameters.
        
        Args:
            name: System name
            
        Returns:
            Formatted string with system information
        """
        system = self.get_system(name)
        if system is None:
            return f"System '{name}' not found in database."
        
        summary = f"\n{'='*60}\n"
        summary += f"System: {system.name}\n"
        summary += f"Category: {system.category}\n"
        summary += f"{'='*60}\n\n"
        
        if system.ra is not None:
            summary += f"Position: RA={system.ra:.4f}°, Dec={system.dec:.4f}°\n"
        if system.distance is not None:
            summary += f"Distance: {system.distance:.2e} pc\n"
        if system.redshift is not None:
            summary += f"Redshift: z={system.redshift}\n"
        
        summary += "\n"
        
        if system.mass is not None:
            summary += f"Mass: {system.mass:.2e} M☉\n"
        if system.radius is not None:
            summary += f"Radius: {system.radius:.2e} (R☉ or pc)\n"
        if system.L_bol is not None:
            summary += f"Bolometric Luminosity: {system.L_bol:.2e} erg/s\n"
        if system.L_X is not None:
            summary += f"X-ray Luminosity: {system.L_X:.2e} erg/s\n"
        if system.T is not None:
            summary += f"Temperature: {system.T:.2e} K\n"
        if system.B is not None:
            summary += f"Magnetic Field: {system.B:.2e} Gauss\n"
        
        if system.sources:
            summary += f"\nSources:\n"
            for source in system.sources:
                summary += f"  - {source}\n"
        
        if system.additional_params:
            summary += f"\nAdditional Parameters:\n"
            for key, value in system.additional_params.items():
                summary += f"  {key}: {value}\n"
        
        return summary
    
    def export_to_csv(self, filename: str, category: Optional[str] = None):
        """
        Export systems database to CSV file.
        
        Args:
            filename: Output CSV filename
            category: Optional category filter
        """
        import csv
        
        if category:
            systems_list = self.get_systems_by_category(category)
        else:
            systems_list = list(self.systems.values())
        
        if not systems_list:
            print(f"No systems to export.")
            return
        
        # Get all possible fields
        fields = [
            'name', 'category', 'ra', 'dec', 'distance', 'redshift',
            'mass', 'radius', 'L_bol', 'L_X', 'T', 'v_exp', 'v_rot', 'v_wind',
            'B', 'n_e', 'n_H', 'pressure', 'SFR', 'spectral_type', 'metallicity'
        ]
        
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fields, extrasaction='ignore')
            writer.writeheader()
            
            for system in systems_list:
                row = {field: getattr(system, field, None) for field in fields}
                writer.writerow(row)
        
        print(f"Exported {len(systems_list)} systems to {filename}")


# ================================
# MAIN EXECUTION & TESTING
# ================================

if __name__ == "__main__":
    # Initialize database
    db = UQFFSystemsDatabase()
    
    print("="*60)
    print("UQFF Systems Database - Initialization Complete")
    print("="*60)
    print(f"\nTotal systems: {len(db.list_systems())}")
    print(f"Categories: {', '.join(db.list_categories())}")
    
    # Category breakdown
    print("\nSystems by Category:")
    for category in db.list_categories():
        systems = db.get_systems_by_category(category)
        print(f"  {category}: {len(systems)} systems")
    
    # Example queries
    print("\n" + "="*60)
    print("Example System Queries")
    print("="*60)
    
    # Query 1: M87
    print(db.get_system_summary('M87'))
    
    # Query 2: SGR 1745
    print(db.get_system_summary('SGR1745'))
    
    # Query 3: AT2024tvd
    print(db.get_system_summary('AT2024tvd'))
    
    # List all TDEs
    print("\n" + "="*60)
    print("All Tidal Disruption Events:")
    print("="*60)
    tdes = db.get_systems_by_category('TDE')
    for tde in tdes:
        print(f"  - {tde.name}: z={tde.redshift}, M_BH={tde.mass:.1e} M☉")
    
    # Export examples
    print("\n" + "="*60)
    print("CSV Export Examples")
    print("="*60)
    # db.export_to_csv('uqff_all_systems.csv')
    # db.export_to_csv('uqff_agn_systems.csv', category='AGN')
    print("  (Export functions available: export_to_csv)")
