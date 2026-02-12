"""
QCalc Validation Data Sources
===============================

Comprehensive URL repository and fetch utilities for validating UQFF physics calculations
against observational/empirical data from major astronomical and nuclear databases.

Architecture:
    - NO HARDCODED SYSTEM DATA (only URLs and fetch utilities)
    - Parameterized queries for any astronomical object
    - Support for 26-level energy structure validation
    - Cross-validation with multiple independent sources

Data Sources:
    1. SIMBAD - Stellar/galactic parameters
    2. NED - Extragalactic database
    3. HEASARC - High-energy astrophysics
    4. Chandra - X-ray observations
    5. Fermi - Gamma-ray observations
    6. GAIA DR4 - Astrometric/photometric data
    7. LIGO - Gravitational wave events
    8. Nuclear Databases - Binding energies, shell structure
    9. Quasar Catalogs - Luminosity, redshift
    10. Solar System Ephemeris - Planetary/solar data

Usage:
    from QCalc_validation import ValidationDataFetcher
    
    fetcher = ValidationDataFetcher()
    # Fetch stellar parameters
    data = fetcher.fetch_simbad("Betelgeuse")
    # Validate 26-level energy
    nuclear_data = fetcher.fetch_nuclear_binding_energies(Z=26, A=56)
    # Check quasar luminosity
    quasar_data = fetcher.fetch_quasar_catalog(name="3C273")
"""

import requests
import json
import logging
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from datetime import datetime
import urllib.parse
import time

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════
# DATA SOURCE URLs
# ═══════════════════════════════════════════════════════════════════════════

class DataSourceURLs:
    """Central repository of all validation data source URLs."""
    
    # ═══════════════════════════════════════════════════════════════════════
    # ASTRONOMICAL DATABASES
    # ═══════════════════════════════════════════════════════════════════════
    
    # SIMBAD - Stellar/Galactic Objects
    SIMBAD_BASE = "http://simbad.u-strasbg.fr/simbad/sim-tap/sync"
    SIMBAD_API = "http://simbad.cds.unistra.fr/simbad/sim-id"
    
    # NED - NASA/IPAC Extragalactic Database
    NED_BASE = "https://ned.ipac.caltech.edu/tap/sync"
    NED_API = "https://ned.ipac.caltech.edu/srs/ObjectLookup"
    
    # HEASARC - High Energy Astrophysics Science Archive
    HEASARC_BASE = "https://heasarc.gsfc.nasa.gov/cgi-bin/vo/cone/coneGet.pl"
    HEASARC_TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap"
    
    # Specific HEASARC Catalogs
    HEASARC_XRAY = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dheasarc_xraybsc&Action=More+Options"
    HEASARC_GAMMA = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dheasarc_fermi4yr&Action=More+Options"
    HEASARC_MAGNETAR = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dmagnetar&Action=More+Options"
    
    # Chandra X-ray Observatory
    CHANDRA_DATA = "https://cda.harvard.edu/csccli/getProperties"
    CHANDRA_CATALOG = "https://cda.harvard.edu/csc2scs/cone"
    
    # Fermi Gamma-ray Space Telescope
    FERMI_LAT = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/"
    FERMI_4FGL = "https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermi4fgl&Action=More+Options"
    
    # GAIA - ESA Astrometry Mission
    GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
    GAIA_DR4 = "https://gea.esac.esa.int/data-server/data"  # DR4 anticipated 2026
    GAIA_DR3 = "https://gea.esac.esa.int/data-server/data"  # Current DR3
    
    # LIGO - Laser Interferometer Gravitational-Wave Observatory
    LIGO_GWOSC = "https://gwosc.org/eventapi/json/GWTC/"
    LIGO_GWTC4 = "https://gwosc.org/eventapi/json/GWTC-4/"  # Latest catalog
    LIGO_CATALOG = "https://gwosc.org/eventapi/html/GWTC/"
    
    # ═══════════════════════════════════════════════════════════════════════
    # NUCLEAR & PARTICLE PHYSICS DATABASES
    # ═══════════════════════════════════════════════════════════════════════
    
    # NIST Atomic Spectra Database
    NIST_ATOMIC = "https://physics.nist.gov/cgi-bin/ASD/energy1.pl"
    
    # NNDC - National Nuclear Data Center
    NNDC_BINDING = "https://www.nndc.bnl.gov/nudat3/"
    NNDC_CHART = "https://www.nndc.bnl.gov/chart/"
    
    # ENSDF - Evaluated Nuclear Structure Data File
    ENSDF = "https://www.nndc.bnl.gov/ensdf/"
    
    # Particle Data Group
    PDG_API = "https://pdglive.lbl.gov/api/v1/"
    PDG_SUMMARY = "https://pdglive.lbl.gov/DataBlock.action"
    
    # IAEA Nuclear Data Services
    IAEA_NUCLEAR = "https://www-nds.iaea.org/relnsd/vcharthtml/VChartHTML.html"
    
    # ═══════════════════════════════════════════════════════════════════════
    # QUASAR & AGN CATALOGS
    # ═══════════════════════════════════════════════════════════════════════
    
    # Sloan Digital Sky Survey (SDSS) Quasar Catalog
    SDSS_QUASAR = "https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch"
    
    # Milliquas - Million Quasars Catalog
    MILLIQUAS = "https://quasars.org/milliquas.htm"
    
    # Veron-Cetty & Veron Quasar Catalog
    VERON_QUASAR = "http://cdsarc.u-strasbg.fr/viz-bin/cat/VII/258"
    
    # ═══════════════════════════════════════════════════════════════════════
    # SOLAR SYSTEM & PLANETARY DATA
    # ═══════════════════════════════════════════════════════════════════════
    
    # JPL Horizons System (ephemeris)
    JPL_HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
    
    # NASA Planetary Fact Sheets
    NASA_PLANETS = "https://nssdc.gsfc.nasa.gov/planetary/factsheet/"
    
    # SOHO Solar Data
    SOHO_DATA = "https://sohowww.nascom.nasa.gov/data/realtime/"
    
    # ═══════════════════════════════════════════════════════════════════════
    # GALAXY ROTATION CURVES & DARK MATTER
    # ═══════════════════════════════════════════════════════════════════════
    
    # SPARC - Spitzer Photometry and Accurate Rotation Curves
    SPARC_DATABASE = "http://astroweb.cwru.edu/SPARC/"
    
    # HyperLEDA - Extragalactic Database
    HYPERLEDA = "http://leda.univ-lyon1.fr/fullsql.html"
    
    # ═══════════════════════════════════════════════════════════════════════
    # COSMOLOGICAL PARAMETERS
    # ═══════════════════════════════════════════════════════════════════════
    
    # Planck Legacy Archive
    PLANCK_ARCHIVE = "https://pla.esac.esa.int/pla/"
    
    # WMAP Legacy Archive
    WMAP_ARCHIVE = "https://lambda.gsfc.nasa.gov/product/map/"


# ═══════════════════════════════════════════════════════════════════════════
# DATA CLASSES
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class ValidationResult:
    """Result of validating UQFF calculation against observational data."""
    source: str                          # Data source (SIMBAD, HEASARC, etc.)
    object_name: str                     # Astronomical object name
    calculated_value: float              # UQFF computed value
    observed_value: Optional[float]      # Observed value from database
    unit: str                            # Unit of measurement
    percent_error: Optional[float]       # Percent difference
    match_status: str                    # "MATCH", "CLOSE", "MISMATCH", "NO_DATA"
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class NuclearData:
    """Nuclear structure data for 26-level validation."""
    Z: int                               # Atomic number
    A: int                               # Mass number
    binding_energy_MeV: Optional[float]  # Total binding energy
    shell_structure: Optional[str]       # Shell configuration
    energy_levels_J: List[float]         # Energy levels in Joules
    source: str                          # Database source


@dataclass
class AstrophysicalData:
    """Astrophysical object data from catalogs."""
    name: str
    object_type: str                     # Star, galaxy, quasar, etc.
    mass_kg: Optional[float]
    radius_m: Optional[float]
    distance_m: Optional[float]
    luminosity_W: Optional[float]
    temperature_K: Optional[float]
    magnetic_field_T: Optional[float]
    redshift: Optional[float]
    source: str
    metadata: Dict[str, Any] = field(default_factory=dict)


# ═══════════════════════════════════════════════════════════════════════════
# VALIDATION DATA FETCHER
# ═══════════════════════════════════════════════════════════════════════════

class ValidationDataFetcher:
    """
    Fetch validation data from multiple sources to verify UQFF calculations.
    
    All methods are parameterized - no hardcoded system data.
    """
    
    def __init__(self, timeout: int = 30, retry_attempts: int = 3):
        """
        Initialize fetcher.
        
        Args:
            timeout: Request timeout in seconds
            retry_attempts: Number of retry attempts on failure
        """
        self.urls = DataSourceURLs()
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'UQFF-QCalc-Validator/1.0 (Scientific Research)'
        })
    
    # ═══════════════════════════════════════════════════════════════════════
    # SIMBAD QUERIES
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_simbad(self, object_name: str) -> Optional[AstrophysicalData]:
        """
        Fetch stellar/galactic data from SIMBAD.
        
        Args:
            object_name: Name of astronomical object (e.g., "Betelgeuse", "M87")
        
        Returns:
            AstrophysicalData object or None
        """
        try:
            # SIMBAD TAP query for comprehensive data
            query = f"""
            SELECT 
                basic.main_id, basic.otype_txt,
                otypes.otype_longname,
                flux.U, flux.B, flux.V,
                mesDistance.distance_result AS distance_pc,
                mesVelocities.rv_value AS radial_velocity
            FROM basic
            LEFT JOIN flux ON basic.oid = flux.oidref
            LEFT JOIN mesDistance ON basic.oid = mesDistance.oidref
            LEFT JOIN mesVelocities ON basic.oid = mesVelocities.oidref
            LEFT JOIN otypes ON basic.otype = otypes.otype
            WHERE basic.main_id = '{object_name}'
            """
            
            params = {
                'REQUEST': 'doQuery',
                'LANG': 'ADQL',
                'FORMAT': 'json',
                'QUERY': query
            }
            
            response = self._retry_request('GET', self.urls.SIMBAD_BASE, params=params)
            
            if response and response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    obj_data = data['data'][0]
                    return AstrophysicalData(
                        name=obj_data.get('main_id', object_name),
                        object_type=obj_data.get('otype_longname', 'Unknown'),
                        mass_kg=None,  # SIMBAD doesn't directly provide mass
                        radius_m=None,
                        distance_m=self._parsec_to_meters(obj_data.get('distance_pc')),
                        luminosity_W=None,
                        temperature_K=None,
                        magnetic_field_T=None,
                        redshift=None,
                        source='SIMBAD',
                        metadata={'raw_data': obj_data}
                    )
            
            logger.warning(f"SIMBAD: No data found for {object_name}")
            return None
            
        except Exception as e:
            logger.error(f"SIMBAD fetch error for {object_name}: {e}")
            return None
    
    # ═══════════════════════════════════════════════════════════════════════
    # NED QUERIES (Extragalactic)
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_ned(self, object_name: str) -> Optional[AstrophysicalData]:
        """
        Fetch extragalactic data from NED.
        
        Args:
            object_name: Name of object (e.g., "NGC 4889", "3C273")
        
        Returns:
            AstrophysicalData object or None
        """
        try:
            params = {
                'objname': object_name,
                'of': 'json'
            }
            
            response = self._retry_request('GET', self.urls.NED_API, params=params)
            
            if response and response.status_code == 200:
                data = response.json()
                return AstrophysicalData(
                    name=data.get('Name', object_name),
                    object_type=data.get('Type', 'Unknown'),
                    mass_kg=None,
                    radius_m=None,
                    distance_m=self._redshift_to_distance(data.get('Redshift')),
                    luminosity_W=None,
                    temperature_K=None,
                    magnetic_field_T=None,
                    redshift=data.get('Redshift'),
                    source='NED',
                    metadata={'raw_data': data}
                )
            
            return None
            
        except Exception as e:
            logger.error(f"NED fetch error for {object_name}: {e}")
            return None
    
    # ═══════════════════════════════════════════════════════════════════════
    # HEASARC QUERIES (High-Energy)
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_heasarc_xray(self, ra_deg: float, dec_deg: float, radius_arcmin: float = 5.0) -> List[Dict]:
        """
        Fetch X-ray sources from HEASARC near coordinates.
        
        Args:
            ra_deg: Right Ascension (degrees)
            dec_deg: Declination (degrees)
            radius_arcmin: Search radius (arcminutes)
        
        Returns:
            List of X-ray sources
        """
        try:
            params = {
                'tablehead': 'name=heasarc_xraybsc',
                'RA': ra_deg,
                'DEC': dec_deg,
                'Radius': radius_arcmin,
                'ResultMax': 100,
                'displaymode': 'JSON'
            }
            
            response = self._retry_request('GET', self.urls.HEASARC_BASE, params=params)
            
            if response and response.status_code == 200:
                return response.json()
            
            return []
            
        except Exception as e:
            logger.error(f"HEASARC X-ray fetch error: {e}")
            return []
    
    # ═══════════════════════════════════════════════════════════════════════
    # GAIA QUERIES
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_gaia(self, object_name: str) -> Optional[Dict]:
        """
        Fetch astrometric data from GAIA DR3/DR4.
        
        Args:
            object_name: Star name or designation
        
        Returns:
            GAIA data dictionary or None
        """
        try:
            # GAIA TAP query
            query = f"""
            SELECT TOP 1
                source_id, ra, dec, parallax, parallax_error,
                pmra, pmdec, radial_velocity,
                phot_g_mean_mag, bp_rp
            FROM gaiadr3.gaia_source
            WHERE CONTAINS(
                POINT('ICRS', ra, dec),
                CIRCLE('ICRS', 
                    (SELECT ra FROM gaiadr3.dr2_neighbourhood WHERE dr2_source_id IN 
                        (SELECT source_id FROM gaiadr3.gaia_source WHERE UPPER(designation) LIKE UPPER('%{object_name}%') LIMIT 1)
                    ), 
                    (SELECT dec FROM gaiadr3.dr2_neighbourhood WHERE dr2_source_id IN 
                        (SELECT source_id FROM gaiadr3.gaia_source WHERE UPPER(designation) LIKE UPPER('%{object_name}%') LIMIT 1)
                    ),
                    0.01
                )
            )=1
            """
            
            params = {
                'REQUEST': 'doQuery',
                'LANG': 'ADQL',
                'FORMAT': 'json',
                'QUERY': query
            }
            
            response = self._retry_request('GET', self.urls.GAIA_TAP, params=params)
            
            if response and response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    return data['data'][0]
            
            return None
            
        except Exception as e:
            logger.error(f"GAIA fetch error for {object_name}: {e}")
            return None
    
    # ═══════════════════════════════════════════════════════════════════════
    # LIGO GRAVITATIONAL WAVE EVENTS
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_ligo_gwtc4(self, event_name: Optional[str] = None) -> List[Dict]:
        """
        Fetch gravitational wave events from LIGO GWTC-4.
        
        Args:
            event_name: Specific event (e.g., "GW150914") or None for all
        
        Returns:
            List of GW events
        """
        try:
            url = self.urls.LIGO_GWTC4
            if event_name:
                url += f"{event_name}/"
            
            response = self._retry_request('GET', url)
            
            if response and response.status_code == 200:
                data = response.json()
                if isinstance(data, dict) and 'events' in data:
                    return data['events']
                elif isinstance(data, list):
                    return data
                else:
                    return [data]
            
            return []
            
        except Exception as e:
            logger.error(f"LIGO GWTC-4 fetch error: {e}")
            return []
    
    # ═══════════════════════════════════════════════════════════════════════
    # NUCLEAR DATABASE QUERIES
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_nuclear_binding_energy(self, Z: int, A: int) -> Optional[NuclearData]:
        """
        Fetch nuclear binding energy from NNDC.
        
        Args:
            Z: Atomic number (protons)
            A: Mass number (protons + neutrons)
        
        Returns:
            NuclearData object or None
        """
        try:
            # Note: NNDC requires web scraping or manual download
            # This is a placeholder for the structure
            logger.info(f"Nuclear data fetch for Z={Z}, A={A}")
            logger.info(f"Manual source: {self.urls.NNDC_BINDING}")
            
            # For 26-level validation, we use theoretical values
            # Real implementation would scrape NNDC or use pre-downloaded data
            
            return NuclearData(
                Z=Z,
                A=A,
                binding_energy_MeV=None,  # Would be fetched
                shell_structure=None,
                energy_levels_J=[],
                source='NNDC (manual)'
            )
            
        except Exception as e:
            logger.error(f"Nuclear data fetch error for Z={Z}, A={A}: {e}")
            return None
    
    # ═══════════════════════════════════════════════════════════════════════
    # QUASAR CATALOG QUERIES
    # ═══════════════════════════════════════════════════════════════════════
    
    def fetch_quasar_sdss(self, quasar_name: str) -> Optional[Dict]:
        """
        Fetch quasar data from SDSS.
        
        Args:
            quasar_name: Quasar designation
        
        Returns:
            Quasar data dictionary or None
        """
        try:
            query = f"""
            SELECT TOP 1
                objID, ra, dec, z AS redshift,
                psfMag_u, psfMag_g, psfMag_r, psfMag_i, psfMag_z
            FROM SpecObj
            WHERE class = 'QSO' AND
                  bestObjID IN (SELECT objID FROM PhotoObj WHERE objID > 0)
            """
            # Note: Would need proper SDSS SQL query for specific object
            
            logger.info(f"SDSS quasar query for {quasar_name}")
            logger.info(f"Source: {self.urls.SDSS_QUASAR}")
            
            return None  # Placeholder
            
        except Exception as e:
            logger.error(f"SDSS quasar fetch error for {quasar_name}: {e}")
            return None
    
    # ═══════════════════════════════════════════════════════════════════════
    # UTILITY METHODS
    # ═══════════════════════════════════════════════════════════════════════
    
    def _retry_request(self, method: str, url: str, **kwargs) -> Optional[requests.Response]:
        """Retry HTTP request with exponential backoff."""
        for attempt in range(self.retry_attempts):
            try:
                response = self.session.request(
                    method, url, timeout=self.timeout, **kwargs
                )
                response.raise_for_status()
                return response
            except requests.exceptions.RequestException as e:
                if attempt < self.retry_attempts - 1:
                    wait_time = 2 ** attempt
                    logger.warning(f"Request failed, retrying in {wait_time}s: {e}")
                    time.sleep(wait_time)
                else:
                    logger.error(f"Request failed after {self.retry_attempts} attempts: {e}")
                    return None
        return None
    
    def _parsec_to_meters(self, parsecs: Optional[float]) -> Optional[float]:
        """Convert parsecs to meters."""
        if parsecs is None:
            return None
        return parsecs * 3.0857e16  # 1 pc = 3.0857e16 m
    
    def _redshift_to_distance(self, z: Optional[float]) -> Optional[float]:
        """Approximate comoving distance from redshift (simplified)."""
        if z is None:
            return None
        # Hubble constant H0 = 67.4 km/s/Mpc
        c = 2.998e8  # m/s
        H0 = 67.4e3 / 3.0857e22  # Convert to s^-1
        # Linear approximation (valid for z << 1)
        distance_m = (c * z) / H0
        return distance_m
    
    def validate_calculation(
        self,
        object_name: str,
        calculated_value: float,
        parameter_name: str,
        unit: str,
        sources: List[str] = ['SIMBAD', 'NED']
    ) -> List[ValidationResult]:
        """
        Validate UQFF calculation against multiple data sources.
        
        Args:
            object_name: Astronomical object name
            calculated_value: UQFF computed value
            parameter_name: What was calculated (e.g., "luminosity", "mass")
            unit: Unit of measurement
            sources: List of data sources to query
        
        Returns:
            List of ValidationResult objects
        """
        results = []
        
        for source in sources:
            if source == 'SIMBAD':
                obs_data = self.fetch_simbad(object_name)
            elif source == 'NED':
                obs_data = self.fetch_ned(object_name)
            else:
                continue
            
            if obs_data:
                # Extract relevant parameter (would need mapping)
                observed_value = None  # Would extract from obs_data
                percent_error = None
                match_status = "NO_DATA"
                
                if observed_value is not None:
                    percent_error = abs(calculated_value - observed_value) / observed_value * 100
                    if percent_error < 5:
                        match_status = "MATCH"
                    elif percent_error < 20:
                        match_status = "CLOSE"
                    else:
                        match_status = "MISMATCH"
                
                results.append(ValidationResult(
                    source=source,
                    object_name=object_name,
                    calculated_value=calculated_value,
                    observed_value=observed_value,
                    unit=unit,
                    percent_error=percent_error,
                    match_status=match_status,
                    metadata={'parameter': parameter_name}
                ))
        
        return results


# ═══════════════════════════════════════════════════════════════════════════
# VALIDATION CAMPAIGN TEMPLATES
# ═══════════════════════════════════════════════════════════════════════════

class ValidationCampaigns:
    """Pre-defined validation campaigns for systematic testing."""
    
    @staticmethod
    def validate_26_level_structure():
        """
        Validate 26-level energy structure against nuclear/cosmic observations.
        
        Test Points:
            n=1-4: Vacuum fluctuations (theoretical)
            n=5-10: Nuclear binding energies (NNDC)
            n=11-13: Molecular bonds (NIST)
            n=14-18: Stellar winds, Higgs (LHC, stellar obs)
            n=19-26: Quasar jets, cosmic rays (HEASARC, Fermi)
        """
        test_cases = [
            # Nuclear regime (n=8 → ~10^-12 J = 6 MeV)
            {'n': 8, 'expected_J': 1e-12, 'source': 'NNDC', 'system': 'Fe-56 binding'},
            
            # Higgs regime (n=18 → ~10^-2 J = 125 GeV)
            {'n': 18, 'expected_J': 1e-2, 'source': 'LHC', 'system': 'Higgs boson'},
            
            # Quasar jet regime (n=22 → ~10^2 J)
            {'n': 22, 'expected_J': 1e2, 'source': 'HEASARC', 'system': '3C273 jet'},
        ]
        return test_cases
    
    @staticmethod
    def validate_reactor_efficiency():
        """
        Validate reactor efficiency model against astrophysical luminosities.
        
        Test Points:
            - Quasar luminosities (10^39-47 W)
            - Magnetar X-ray emission
            - Planetary core heat flow
        """
        test_cases = [
            {'object': '3C273', 'expected_L_W': 1e40, 'source': 'SDSS'},
            {'object': 'SGR1745', 'expected_L_W': 1e35, 'source': 'HEASARC'},
        ]
        return test_cases
    
    @staticmethod
    def validate_galaxy_rotation():
        """
        Validate UQFF gravity against SPARC galaxy rotation curves.
        
        Test Points:
            - 175 galaxies from SPARC database
            - Compare UQFF vs. Newtonian vs. MOND
        """
        test_cases = [
            {'galaxy': 'NGC3198', 'radius_kpc': [1, 5, 10, 20], 'source': 'SPARC'},
            {'galaxy': 'UGC2885', 'radius_kpc': [10, 30, 50], 'source': 'SPARC'},
        ]
        return test_cases


# ═══════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("═" * 80)
    print("QCalc Validation Data Sources - Test Suite")
    print("═" * 80)
    
    fetcher = ValidationDataFetcher()
    
    # Test SIMBAD fetch
    print("\n[TEST 1] SIMBAD Fetch - Betelgeuse")
    print("-" * 80)
    betelgeuse = fetcher.fetch_simbad("Betelgeuse")
    if betelgeuse:
        print(f"✓ Name: {betelgeuse.name}")
        print(f"✓ Type: {betelgeuse.object_type}")
        print(f"✓ Distance: {betelgeuse.distance_m:.3e} m" if betelgeuse.distance_m else "  Distance: N/A")
    else:
        print("✗ No data returned")
    
    # Test NED fetch
    print("\n[TEST 2] NED Fetch - M87")
    print("-" * 80)
    m87 = fetcher.fetch_ned("M87")
    if m87:
        print(f"✓ Name: {m87.name}")
        print(f"✓ Type: {m87.object_type}")
        print(f"✓ Redshift: {m87.redshift}" if m87.redshift else "  Redshift: N/A")
    else:
        print("✗ No data returned")
    
    # Test LIGO fetch
    print("\n[TEST 3] LIGO GWTC-4 Fetch")
    print("-" * 80)
    gw_events = fetcher.fetch_ligo_gwtc4()
    if gw_events:
        print(f"✓ Retrieved {len(gw_events)} gravitational wave events")
        if len(gw_events) > 0:
            print(f"  Example: {gw_events[0].get('name', 'Unknown')}")
    else:
        print("✗ No events returned")
    
    # Display URL repository
    print("\n[DATA SOURCE REPOSITORY]")
    print("-" * 80)
    print(f"SIMBAD:    {fetcher.urls.SIMBAD_BASE}")
    print(f"NED:       {fetcher.urls.NED_BASE}")
    print(f"HEASARC:   {fetcher.urls.HEASARC_BASE}")
    print(f"GAIA DR3:  {fetcher.urls.GAIA_DR3}")
    print(f"LIGO:      {fetcher.urls.LIGO_GWTC4}")
    print(f"NNDC:      {fetcher.urls.NNDC_BINDING}")
    print(f"SPARC:     {fetcher.urls.SPARC_DATABASE}")
    
    print("\n" + "═" * 80)
    print("Validation infrastructure ready. Use fetcher methods to validate QCalc results.")
    print("═" * 80)
