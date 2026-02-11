#!/usr/bin/env python3
"""
APIFetch.py - UQFF API Data Fetching Layer (25 APIs)
=====================================================

Fetches astronomical parameters from external APIs for UQFF calculations.

ARCHITECTURE:
    User Query → APIFetch.py (this file) → IPData.py → QCalc.py → OPData.py
    
25 SUPPORTED APIs:
    ─────────────────────────────────────────────────────────────────────────
    ASTRONOMICAL DATABASES (Primary)
    ─────────────────────────────────────────────────────────────────────────
     1. SIMBAD       - CDS Strasbourg (~15M objects)
     2. NED          - NASA/IPAC Extragalactic Database (~300M objects)
     3. VizieR       - CDS catalog service (20,000+ catalogs)
     4. Gaia         - ESA astrometry mission (1.8B stars)
     5. MAST         - Mikulski Archive (HST, JWST, etc.)
    ─────────────────────────────────────────────────────────────────────────
    NASA/SPACE AGENCIES
    ─────────────────────────────────────────────────────────────────────────
     6. NASA Exoplanet Archive - Confirmed exoplanets
     7. NASA HEASARC - High Energy Astrophysics data
     8. NASA ADS     - Astrophysics Data System (literature)
     9. NASA JPL Horizons - Solar system ephemerides
    10. ESO Archive  - European Southern Observatory
    ─────────────────────────────────────────────────────────────────────────
    SURVEY DATABASES
    ─────────────────────────────────────────────────────────────────────────
    11. SDSS         - Sloan Digital Sky Survey
    12. 2MASS        - Two Micron All Sky Survey
    13. WISE         - Wide-field Infrared Survey Explorer
    14. Pan-STARRS   - Panoramic Survey Telescope
    15. ZTF          - Zwicky Transient Facility
    ─────────────────────────────────────────────────────────────────────────
    SPECIALIZED DATABASES
    ─────────────────────────────────────────────────────────────────────────
    16. ATNF         - Pulsar catalog
    17. McGill       - Magnetar catalog
    18. TNS          - Transient Name Server (supernovae/transients)
    19. GCN          - Gamma-ray Coordinates Network
    20. LIGO/Virgo   - Gravitational wave events
    ─────────────────────────────────────────────────────────────────────────
    COMPUTATIONAL/AI
    ─────────────────────────────────────────────────────────────────────────
    21. Wolfram Alpha - Computational knowledge engine
    22. arXiv        - Preprint server (metadata)
    23. Grok (xAI)   - AI fallback for missing data
    24. OpenAI       - Alternative AI fallback
    25. Claude       - Alternative AI fallback
    ─────────────────────────────────────────────────────────────────────────

OUTPUT: 
    All results written to IPData.py (InputParameters dataclass)
    
Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import requests
import json
import os
import csv
from datetime import datetime
from typing import Dict, List, Any, Optional, Tuple
import re

# Import IPData for storing results
from IPData import InputParameters, InputDataStore, INPUT_STORE, store_input

# ═══════════════════════════════════════════════════════════════════════════════
# UNIT CONVERSIONS (for standardizing API responses to SI units)
# ═══════════════════════════════════════════════════════════════════════════════

UNITS = {
    # Mass
    'M_sun': 1.989e30,         # Solar mass (kg)
    'M_earth': 5.972e24,       # Earth mass (kg)
    'M_jupiter': 1.898e27,     # Jupiter mass (kg)
    
    # Distance
    'pc': 3.086e16,            # Parsec (m)
    'kpc': 3.086e19,           # Kiloparsec (m)
    'Mpc': 3.086e22,           # Megaparsec (m)
    'ly': 9.461e15,            # Light-year (m)
    'AU': 1.496e11,            # Astronomical Unit (m)
    
    # Radius
    'R_sun': 6.96e8,           # Solar radius (m)
    'R_earth': 6.371e6,        # Earth radius (m)
    'R_jupiter': 6.9911e7,     # Jupiter radius (m)
    
    # Luminosity
    'L_sun': 3.828e26,         # Solar luminosity (W)
    
    # Temperature
    'K': 1.0,                  # Kelvin (base unit)
    
    # Magnetic Field
    'T': 1.0,                  # Tesla (base unit)
    'G': 1e-4,                 # Gauss to Tesla
}


# ═══════════════════════════════════════════════════════════════════════════════
# API ENDPOINTS (50 APIs Total)
# ═══════════════════════════════════════════════════════════════════════════════

ENDPOINTS = {
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP A: ASTRONOMICAL DATABASES (1-5)
    # ═══════════════════════════════════════════════════════════════════════════
    'simbad': 'https://simbad.u-strasbg.fr/simbad/sim-tap/sync',
    'ned': 'https://ned.ipac.caltech.edu/tap/sync',
    'vizier': 'https://vizier.u-strasbg.fr/viz-bin/votable',
    'gaia': 'https://gea.esac.esa.int/tap-server/tap/sync',
    'mast': 'https://mast.stsci.edu/api/v0/invoke',
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP B: NASA/SPACE AGENCIES (6-10)
    # ═══════════════════════════════════════════════════════════════════════════
    'nasa_exoplanet': 'https://exoplanetarchive.ipac.caltech.edu/TAP/sync',
    'heasarc': 'https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync',
    'ads': 'https://api.adsabs.harvard.edu/v1/search/query',
    'jpl_horizons': 'https://ssd.jpl.nasa.gov/api/horizons.api',
    'eso': 'https://archive.eso.org/tap_obs/sync',
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP C: SURVEY DATABASES (11-15)
    # ═══════════════════════════════════════════════════════════════════════════
    'sdss': 'https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch',
    '2mass': 'https://irsa.ipac.caltech.edu/TAP/sync',
    'wise': 'https://irsa.ipac.caltech.edu/TAP/sync',
    'panstarrs': 'https://catalogs.mast.stsci.edu/api/v0.1/panstarrs',
    'ztf': 'https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves',
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP D: SPECIALIZED DATABASES (16-20)
    # ═══════════════════════════════════════════════════════════════════════════
    'atnf_pulsar': 'https://www.atnf.csiro.au/research/pulsar/psrcat/proc_form.php',
    'mcgill_magnetar': 'http://www.physics.mcgill.ca/~pulsar/magnetar/TabO1.csv',
    'tns': 'https://www.wis-tns.org/api/get',
    'gcn': 'https://gcn.nasa.gov/api',
    'gwosc': 'https://www.gw-openscience.org/eventapi/json/GWTC/',
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP E: COMPUTATIONAL/AI (21-25)
    # ═══════════════════════════════════════════════════════════════════════════
    'wolfram': 'https://api.wolframalpha.com/v2/query',
    'arxiv': 'http://export.arxiv.org/api/query',
    'grok': 'https://api.x.ai/v1/chat/completions',
    'openai': 'https://api.openai.com/v1/chat/completions',
    'claude': 'https://api.anthropic.com/v1/messages',
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP F: RADIO/INFRARED SURVEYS (26-30)
    # ═══════════════════════════════════════════════════════════════════════════
    'nvss': 'https://www.cv.nrao.edu/nvss/NVSSlist.shtml',           # NRAO VLA Sky Survey
    'first': 'https://sundog.stsci.edu/cgi-bin/searchfirst',         # Faint Images of the Radio Sky
    'vlass': 'https://archive-new.nrao.edu/vlass/quicklook/',        # VLA Sky Survey
    'alma': 'https://almascience.eso.org/aq/',                       # ALMA Science Archive
    'askap': 'https://casda.csiro.au/casda_vo_tools/tap',            # ASKAP Archive
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP G: X-RAY/GAMMA-RAY (31-35)
    # ═══════════════════════════════════════════════════════════════════════════
    'chandra': 'https://cda.harvard.edu/csccli/index',               # Chandra X-ray Observatory
    'xmm': 'http://nxsa.esac.esa.int/nxsa-web/nxsa-tap/tap',         # XMM-Newton Archive
    'swift': 'https://www.swift.ac.uk/swift_live/index.php',         # Swift Gamma-Ray Burst Mission
    'fermi': 'https://fermi.gsfc.nasa.gov/ssc/data/access/',         # Fermi Gamma-ray Space Telescope
    'integral': 'https://www.isdc.unige.ch/integral/archive',        # INTEGRAL Gamma-ray Observatory
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP H: SPACE TELESCOPES (36-40)
    # ═══════════════════════════════════════════════════════════════════════════
    'hst': 'https://mast.stsci.edu/search/hst/api/v0.1/retrieve_product', # Hubble
    'jwst': 'https://mast.stsci.edu/search/jwst/api/v0.1/retrieve_product', # James Webb
    'spitzer': 'https://sha.ipac.caltech.edu/applications/Spitzer/SHA/', # Spitzer (archived)
    'kepler': 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI', # Kepler
    'tess': 'https://mast.stsci.edu/tesscut/api/v0.1',               # TESS
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP I: COSMOLOGY/CMB (41-45)
    # ═══════════════════════════════════════════════════════════════════════════
    'planck': 'https://pla.esac.esa.int/pla/aio/product-action',     # Planck CMB mission
    'wmap': 'https://lambda.gsfc.nasa.gov/product/wmap/',            # WMAP (archived)
    'des': 'https://des.ncsa.illinois.edu/desaccess/',               # Dark Energy Survey
    'desi': 'https://data.desi.lbl.gov/',                            # Dark Energy Spectroscopic Instrument
    'euclid': 'https://easotf.esac.esa.int/tap-server/tap',          # ESA Euclid mission
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GROUP J: SPECTROSCOPIC SURVEYS (46-50)
    # ═══════════════════════════════════════════════════════════════════════════
    'lamost': 'http://dr7.lamost.org/api/',                          # LAMOST Spectroscopic Survey
    'galah': 'https://www.galah-survey.org/dr3/query/',              # GALAH Survey
    'apogee': 'https://data.sdss.org/sas/dr17/apogee/',              # APOGEE Survey
    'rave': 'https://www.rave-survey.org/query/',                    # RAVE Survey
    'desi_spectra': 'https://data.desi.lbl.gov/desi/spectro/',       # DESI Spectra
}

# API Status tracking
API_STATUS = {name: {'enabled': True, 'last_call': None, 'errors': 0} for name in ENDPOINTS}


# ═══════════════════════════════════════════════════════════════════════════════
# SIMBAD FETCHER
# ═══════════════════════════════════════════════════════════════════════════════

class SIMBADFetcher:
    """
    Fetch astronomical data from SIMBAD database.
    
    SIMBAD (Set of Identifications, Measurements and Bibliography for 
    Astronomical Data) contains data for ~15 million objects.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['simbad']
        self.timeout = 30
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """
        Fetch data for an astronomical object from SIMBAD.
        
        Args:
            object_name: Name of the object (e.g., "Sagittarius A*", "Betelgeuse")
            
        Returns:
            Dictionary with available parameters or None if not found
        """
        # TAP query for basic data
        query = f"""
        SELECT TOP 1 
            main_id, ra, dec, pmra, pmdec, plx_value, rvz_radvel,
            sp_type, flux_B, flux_V, flux_J, flux_H, flux_K
        FROM basic
        WHERE main_id = '{object_name}' 
           OR ident LIKE '%{object_name}%'
        """
        
        try:
            response = requests.post(
                self.endpoint,
                data={
                    'REQUEST': 'doQuery',
                    'LANG': 'ADQL',
                    'FORMAT': 'json',
                    'QUERY': query
                },
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    return self._parse_response(data['data'][0], object_name)
            
            return None
            
        except Exception as e:
            print(f"SIMBAD fetch error: {e}")
            return None
    
    def _parse_response(self, row: List, object_name: str) -> Dict[str, Any]:
        """Parse SIMBAD response into standardized format."""
        result = {'name': object_name, 'source': 'SIMBAD'}
        
        # Parallax to distance
        plx = row[5] if len(row) > 5 else None
        if plx and plx > 0:
            # Distance in parsecs = 1000 / parallax(mas)
            distance_pc = 1000.0 / plx
            result['distance'] = distance_pc * UNITS['pc']
        
        # Radial velocity
        if len(row) > 6 and row[6]:
            result['radial_velocity'] = row[6] * 1000  # km/s to m/s
        
        # Spectral type
        if len(row) > 7 and row[7]:
            result['spectral_type'] = row[7]
            # Estimate temperature from spectral type
            result['temperature'] = self._spectral_to_temperature(row[7])
        
        return result
    
    def _spectral_to_temperature(self, spectral_type: str) -> Optional[float]:
        """Estimate effective temperature from spectral type."""
        if not spectral_type:
            return None
        
        # Simplified spectral type to temperature mapping
        spectral_temps = {
            'O': 35000, 'B': 20000, 'A': 9000, 'F': 7000,
            'G': 5500, 'K': 4500, 'M': 3200, 'L': 2000, 'T': 1200
        }
        
        first_char = spectral_type[0].upper()
        return spectral_temps.get(first_char)


# ═══════════════════════════════════════════════════════════════════════════════
# NED FETCHER
# ═══════════════════════════════════════════════════════════════════════════════

class NEDFetcher:
    """
    Fetch extragalactic data from NASA/IPAC Extragalactic Database.
    
    NED contains data on ~300 million objects including galaxies,
    quasars, black holes, and galaxy clusters.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['ned']
        self.timeout = 30
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """
        Fetch data for an extragalactic object from NED.
        
        Args:
            object_name: Name of the object (e.g., "M87", "NGC 1365")
            
        Returns:
            Dictionary with available parameters or None if not found
        """
        query = f"""
        SELECT TOP 1 
            objname, ra, dec, z, Dist_Mpc
        FROM objdir
        WHERE objname LIKE '%{object_name}%'
        """
        
        try:
            response = requests.post(
                self.endpoint,
                data={
                    'REQUEST': 'doQuery',
                    'LANG': 'ADQL',
                    'FORMAT': 'json',
                    'QUERY': query
                },
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    return self._parse_response(data['data'][0], object_name)
            
            return None
            
        except Exception as e:
            print(f"NED fetch error: {e}")
            return None
    
    def _parse_response(self, row: List, object_name: str) -> Dict[str, Any]:
        """Parse NED response into standardized format."""
        result = {'name': object_name, 'source': 'NED'}
        
        # Redshift
        if len(row) > 3 and row[3]:
            result['redshift'] = row[3]
        
        # Distance in Mpc
        if len(row) > 4 and row[4]:
            result['distance'] = row[4] * UNITS['Mpc']
        
        return result


# ═══════════════════════════════════════════════════════════════════════════════
# GROK AI FETCHER (Fallback)
# ═══════════════════════════════════════════════════════════════════════════════

class GrokFetcher:
    """
    Fallback fetcher using Grok AI to estimate missing parameters.
    
    Used when SIMBAD/NED don't have complete data.
    Requires XAI_API_KEY environment variable.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['grok']
        self.api_key = os.environ.get('XAI_API_KEY', '')
        self.timeout = 60
    
    def fetch(self, object_name: str, missing_params: List[str] = None) -> Optional[Dict[str, Any]]:
        """
        Use Grok AI to estimate astronomical parameters.
        
        Args:
            object_name: Name of the astronomical object
            missing_params: List of parameters to estimate
            
        Returns:
            Dictionary with estimated parameters or None if failed
        """
        if not self.api_key:
            print("Warning: XAI_API_KEY not set, Grok fallback unavailable")
            return None
        
        if missing_params is None:
            missing_params = ['mass', 'distance', 'radius', 'temperature', 'luminosity']
        
        prompt = f"""
        For the astronomical object "{object_name}", provide the following parameters 
        in JSON format with SI units:
        
        Parameters needed: {', '.join(missing_params)}
        
        Use these units:
        - mass: kg
        - distance: meters
        - radius: meters
        - temperature: Kelvin
        - luminosity: Watts
        - magnetic_field: Tesla
        - velocity_dispersion: m/s
        
        Respond with ONLY valid JSON, no explanation.
        Example format:
        {{"mass": 1.989e30, "distance": 2.5e20, "temperature": 5778}}
        """
        
        try:
            response = requests.post(
                self.endpoint,
                headers={
                    'Authorization': f'Bearer {self.api_key}',
                    'Content-Type': 'application/json'
                },
                json={
                    'model': 'grok-beta',
                    'messages': [
                        {'role': 'system', 'content': 'You are an astrophysics expert. Provide only JSON responses.'},
                        {'role': 'user', 'content': prompt}
                    ],
                    'temperature': 0.1
                },
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                content = data['choices'][0]['message']['content']
                # Parse JSON from response
                result = json.loads(content)
                result['name'] = object_name
                result['source'] = 'Grok'
                return result
            
            return None
            
        except Exception as e:
            print(f"Grok fetch error: {e}")
            return None


# ═══════════════════════════════════════════════════════════════════════════════
# API FETCHER PLACEHOLDERS (3-25)
# ═══════════════════════════════════════════════════════════════════════════════
# These are placeholder classes for future API implementations.
# Each follows the same interface: fetch(object_name) -> Optional[Dict[str, Any]]

class BaseFetcher:
    """Base class for all API fetchers."""
    
    def __init__(self, api_name: str):
        self.api_name = api_name
        self.endpoint = ENDPOINTS.get(api_name, '')
        self.timeout = 30
        self.enabled = True
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Override in subclass."""
        raise NotImplementedError(f"{self.api_name} fetcher not yet implemented")
    
    def _log_error(self, e: Exception):
        API_STATUS[self.api_name]['errors'] += 1
        print(f"{self.api_name} fetch error: {e}")


# ─── API 3: VizieR ──────────────────────────────────────────────────────────────

class VizieRFetcher(BaseFetcher):
    """VizieR catalog service (20,000+ astronomical catalogs)."""
    
    def __init__(self):
        super().__init__('vizier')
    
    def fetch(self, object_name: str, catalog: str = None) -> Optional[Dict[str, Any]]:
        """Fetch from VizieR catalogs. TODO: Implement."""
        # PLACEHOLDER: Implement VizieR TAP query
        return None


# ─── API 4: Gaia ────────────────────────────────────────────────────────────────

class GaiaFetcher(BaseFetcher):
    """ESA Gaia mission (1.8 billion stars with astrometry)."""
    
    def __init__(self):
        super().__init__('gaia')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from Gaia DR3. TODO: Implement."""
        # PLACEHOLDER: Implement Gaia TAP query
        # Returns: parallax, proper motion, radial velocity, photometry
        return None


# ─── API 5: MAST ────────────────────────────────────────────────────────────────

class MASTFetcher(BaseFetcher):
    """Mikulski Archive for Space Telescopes (HST, JWST, TESS, Kepler)."""
    
    def __init__(self):
        super().__init__('mast')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from MAST. TODO: Implement."""
        # PLACEHOLDER: Implement MAST API query
        return None


# ─── API 6: NASA Exoplanet Archive ──────────────────────────────────────────────

class ExoplanetFetcher(BaseFetcher):
    """NASA Exoplanet Archive (confirmed exoplanets and host stars)."""
    
    def __init__(self):
        super().__init__('nasa_exoplanet')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch exoplanet data. TODO: Implement."""
        # PLACEHOLDER: Implement Exoplanet Archive TAP query
        # Returns: planet mass, radius, orbital period, host star properties
        return None


# ─── API 7: HEASARC ─────────────────────────────────────────────────────────────

class HEASARCFetcher(BaseFetcher):
    """NASA High Energy Astrophysics Science Archive (X-ray, gamma-ray)."""
    
    def __init__(self):
        super().__init__('heasarc')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch high-energy data. TODO: Implement."""
        # PLACEHOLDER: X-ray flux, spectrum, variability
        return None


# ─── API 8: NASA ADS ────────────────────────────────────────────────────────────

class ADSFetcher(BaseFetcher):
    """NASA Astrophysics Data System (literature search)."""
    
    def __init__(self):
        super().__init__('ads')
        self.api_key = os.environ.get('ADS_API_KEY', '')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch published parameters from literature. TODO: Implement."""
        # PLACEHOLDER: Search ADS for object, extract parameters from abstracts
        return None


# ─── API 9: JPL Horizons ────────────────────────────────────────────────────────

class HorizonsFetcher(BaseFetcher):
    """JPL Horizons (solar system ephemerides)."""
    
    def __init__(self):
        super().__init__('jpl_horizons')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch solar system object ephemeris. TODO: Implement."""
        # PLACEHOLDER: Position, velocity, orbital elements
        return None


# ─── API 10: ESO Archive ────────────────────────────────────────────────────────

class ESOFetcher(BaseFetcher):
    """European Southern Observatory Archive."""
    
    def __init__(self):
        super().__init__('eso')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from ESO archive. TODO: Implement."""
        # PLACEHOLDER: Spectroscopy, imaging from VLT, ALMA
        return None


# ─── API 11: SDSS ───────────────────────────────────────────────────────────────

class SDSSFetcher(BaseFetcher):
    """Sloan Digital Sky Survey (photometry, spectroscopy)."""
    
    def __init__(self):
        super().__init__('sdss')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from SDSS. TODO: Implement."""
        # PLACEHOLDER: ugriz photometry, spectroscopic redshift
        return None


# ─── API 12: 2MASS ──────────────────────────────────────────────────────────────

class TwoMASSFetcher(BaseFetcher):
    """Two Micron All Sky Survey (near-infrared JHK photometry)."""
    
    def __init__(self):
        super().__init__('2mass')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from 2MASS. TODO: Implement."""
        # PLACEHOLDER: JHK magnitudes
        return None


# ─── API 13: WISE ───────────────────────────────────────────────────────────────

class WISEFetcher(BaseFetcher):
    """Wide-field Infrared Survey Explorer (mid-infrared)."""
    
    def __init__(self):
        super().__init__('wise')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from WISE. TODO: Implement."""
        # PLACEHOLDER: W1-W4 magnitudes
        return None


# ─── API 14: Pan-STARRS ─────────────────────────────────────────────────────────

class PanSTARRSFetcher(BaseFetcher):
    """Panoramic Survey Telescope and Rapid Response System."""
    
    def __init__(self):
        super().__init__('panstarrs')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch from Pan-STARRS. TODO: Implement."""
        # PLACEHOLDER: grizy photometry
        return None


# ─── API 15: ZTF ────────────────────────────────────────────────────────────────

class ZTFFetcher(BaseFetcher):
    """Zwicky Transient Facility (time-domain astronomy)."""
    
    def __init__(self):
        super().__init__('ztf')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch light curves from ZTF. TODO: Implement."""
        # PLACEHOLDER: Light curves, variability
        return None


# ─── API 16: ATNF Pulsar Catalog ────────────────────────────────────────────────

class ATNFPulsarFetcher(BaseFetcher):
    """Australia Telescope National Facility Pulsar Catalog."""
    
    def __init__(self):
        super().__init__('atnf_pulsar')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch pulsar data. TODO: Implement."""
        # PLACEHOLDER: Period, period derivative, DM, magnetic field
        return None


# ─── API 17: McGill Magnetar Catalog ────────────────────────────────────────────

class McGillMagnetarFetcher(BaseFetcher):
    """McGill Online Magnetar Catalog."""
    
    def __init__(self):
        super().__init__('mcgill_magnetar')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch magnetar data. TODO: Implement."""
        # PLACEHOLDER: B field, period, burst history
        return None


# ─── API 18: Transient Name Server ──────────────────────────────────────────────

class TNSFetcher(BaseFetcher):
    """IAU Transient Name Server (supernovae, transients)."""
    
    def __init__(self):
        super().__init__('tns')
        self.api_key = os.environ.get('TNS_API_KEY', '')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch transient data. TODO: Implement."""
        # PLACEHOLDER: Discovery info, classification, host galaxy
        return None


# ─── API 19: GCN ────────────────────────────────────────────────────────────────

class GCNFetcher(BaseFetcher):
    """Gamma-ray Coordinates Network (GRB alerts)."""
    
    def __init__(self):
        super().__init__('gcn')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch GRB data. TODO: Implement."""
        # PLACEHOLDER: GRB position, fluence, duration
        return None


# ─── API 20: LIGO/Virgo GWOSC ───────────────────────────────────────────────────

class GWOSCFetcher(BaseFetcher):
    """Gravitational Wave Open Science Center."""
    
    def __init__(self):
        super().__init__('gwosc')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch GW event data. TODO: Implement."""
        # PLACEHOLDER: Masses, spins, distance, SNR
        return None


# ─── API 21: Wolfram Alpha ──────────────────────────────────────────────────────

class WolframFetcher(BaseFetcher):
    """Wolfram Alpha computational knowledge engine."""
    
    def __init__(self):
        super().__init__('wolfram')
        self.api_key = os.environ.get('WOLFRAM_APP_ID', '')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Query Wolfram Alpha. TODO: Implement."""
        # PLACEHOLDER: Computed astronomical properties
        return None


# ─── API 22: arXiv ──────────────────────────────────────────────────────────────

class ArXivFetcher(BaseFetcher):
    """arXiv preprint server (metadata extraction)."""
    
    def __init__(self):
        super().__init__('arxiv')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Search arXiv for object. TODO: Implement."""
        # PLACEHOLDER: Recent papers, extract parameters from abstracts
        return None


# ─── API 24: OpenAI (Fallback) ──────────────────────────────────────────────────

class OpenAIFetcher(BaseFetcher):
    """OpenAI GPT fallback for missing data."""
    
    def __init__(self):
        super().__init__('openai')
        self.api_key = os.environ.get('OPENAI_API_KEY', '')
    
    def fetch(self, object_name: str, missing_params: List[str] = None) -> Optional[Dict[str, Any]]:
        """Use OpenAI to estimate parameters. TODO: Implement."""
        # PLACEHOLDER: Similar to GrokFetcher
        return None


# ─── API 25: Claude (Fallback) ──────────────────────────────────────────────────

class ClaudeFetcher(BaseFetcher):
    """Anthropic Claude fallback for missing data."""
    
    def __init__(self):
        super().__init__('claude')
        self.api_key = os.environ.get('ANTHROPIC_API_KEY', '')
    
    def fetch(self, object_name: str, missing_params: List[str] = None) -> Optional[Dict[str, Any]]:
        """Use Claude to estimate parameters. TODO: Implement."""
        # PLACEHOLDER: Similar to GrokFetcher
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# API FETCHER PLACEHOLDERS (26-50) - Additional APIs
# ═══════════════════════════════════════════════════════════════════════════════

# ─── GROUP F: RADIO/INFRARED SURVEYS (26-30) ────────────────────────────────────

class NVSSFetcher(BaseFetcher):
    """NRAO VLA Sky Survey (1.4 GHz radio continuum)."""
    
    def __init__(self):
        super().__init__('nvss')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch radio continuum data. TODO: Implement."""
        # PLACEHOLDER: 1.4 GHz flux density, polarization
        return None


class FIRSTFetcher(BaseFetcher):
    """Faint Images of the Radio Sky at Twenty-cm."""
    
    def __init__(self):
        super().__init__('first')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch high-resolution radio data. TODO: Implement."""
        # PLACEHOLDER: Radio flux, morphology
        return None


class VLASSFetcher(BaseFetcher):
    """VLA Sky Survey (2-4 GHz, high resolution)."""
    
    def __init__(self):
        super().__init__('vlass')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch VLASS data. TODO: Implement."""
        # PLACEHOLDER: S-band continuum
        return None


class ALMAFetcher(BaseFetcher):
    """ALMA Science Archive (mm/submm)."""
    
    def __init__(self):
        super().__init__('alma')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch ALMA observations. TODO: Implement."""
        # PLACEHOLDER: Molecular line data, continuum
        return None


class ASKAPFetcher(BaseFetcher):
    """Australian Square Kilometre Array Pathfinder."""
    
    def __init__(self):
        super().__init__('askap')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch ASKAP data. TODO: Implement."""
        # PLACEHOLDER: Radio continuum, HI surveys
        return None


# ─── GROUP G: X-RAY/GAMMA-RAY (31-35) ───────────────────────────────────────────

class ChandraFetcher(BaseFetcher):
    """Chandra X-ray Observatory Archive."""
    
    def __init__(self):
        super().__init__('chandra')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Chandra X-ray data. TODO: Implement."""
        # PLACEHOLDER: X-ray flux, spectrum, imaging
        return None


class XMMFetcher(BaseFetcher):
    """XMM-Newton Science Archive."""
    
    def __init__(self):
        super().__init__('xmm')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch XMM-Newton data. TODO: Implement."""
        # PLACEHOLDER: X-ray spectra, light curves
        return None


class SwiftFetcher(BaseFetcher):
    """Swift Gamma-Ray Burst Mission."""
    
    def __init__(self):
        super().__init__('swift')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Swift data. TODO: Implement."""
        # PLACEHOLDER: BAT, XRT, UVOT data
        return None


class FermiFetcher(BaseFetcher):
    """Fermi Gamma-ray Space Telescope."""
    
    def __init__(self):
        super().__init__('fermi')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Fermi data. TODO: Implement."""
        # PLACEHOLDER: LAT/GBM gamma-ray flux
        return None


class INTEGRALFetcher(BaseFetcher):
    """INTEGRAL Gamma-ray Observatory."""
    
    def __init__(self):
        super().__init__('integral')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch INTEGRAL data. TODO: Implement."""
        # PLACEHOLDER: Hard X-ray/gamma-ray
        return None


# ─── GROUP H: SPACE TELESCOPES (36-40) ──────────────────────────────────────────

class HSTFetcher(BaseFetcher):
    """Hubble Space Telescope Archive."""
    
    def __init__(self):
        super().__init__('hst')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch HST data. TODO: Implement."""
        # PLACEHOLDER: Optical/UV imaging and spectroscopy
        return None


class JWSTFetcher(BaseFetcher):
    """James Webb Space Telescope Archive."""
    
    def __init__(self):
        super().__init__('jwst')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch JWST data. TODO: Implement."""
        # PLACEHOLDER: Infrared imaging and spectroscopy
        return None


class SpitzerFetcher(BaseFetcher):
    """Spitzer Space Telescope Archive (legacy)."""
    
    def __init__(self):
        super().__init__('spitzer')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Spitzer data. TODO: Implement."""
        # PLACEHOLDER: Mid/far-infrared photometry
        return None


class KeplerFetcher(BaseFetcher):
    """Kepler/K2 Mission Archive."""
    
    def __init__(self):
        super().__init__('kepler')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Kepler data. TODO: Implement."""
        # PLACEHOLDER: Light curves, stellar properties
        return None


class TESSFetcher(BaseFetcher):
    """TESS (Transiting Exoplanet Survey Satellite)."""
    
    def __init__(self):
        super().__init__('tess')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch TESS data. TODO: Implement."""
        # PLACEHOLDER: Light curves, transit signals
        return None


# ─── GROUP I: COSMOLOGY/CMB (41-45) ─────────────────────────────────────────────

class PlanckFetcher(BaseFetcher):
    """Planck Cosmic Microwave Background Mission."""
    
    def __init__(self):
        super().__init__('planck')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Planck data. TODO: Implement."""
        # PLACEHOLDER: CMB, SZ effect, dust emission
        return None


class WMAPFetcher(BaseFetcher):
    """Wilkinson Microwave Anisotropy Probe (legacy)."""
    
    def __init__(self):
        super().__init__('wmap')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch WMAP data. TODO: Implement."""
        # PLACEHOLDER: CMB maps
        return None


class DESFetcher(BaseFetcher):
    """Dark Energy Survey."""
    
    def __init__(self):
        super().__init__('des')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch DES data. TODO: Implement."""
        # PLACEHOLDER: Deep optical/NIR photometry
        return None


class DESIFetcher(BaseFetcher):
    """Dark Energy Spectroscopic Instrument."""
    
    def __init__(self):
        super().__init__('desi')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch DESI data. TODO: Implement."""
        # PLACEHOLDER: Redshifts, spectroscopy
        return None


class EuclidFetcher(BaseFetcher):
    """ESA Euclid Space Mission."""
    
    def __init__(self):
        super().__init__('euclid')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Euclid data. TODO: Implement."""
        # PLACEHOLDER: Optical/NIR imaging, photometric redshifts
        return None


# ─── GROUP J: SPECTROSCOPIC SURVEYS (46-50) ─────────────────────────────────────

class LAMOSTFetcher(BaseFetcher):
    """LAMOST Spectroscopic Survey (4000+ fibers)."""
    
    def __init__(self):
        super().__init__('lamost')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch LAMOST data. TODO: Implement."""
        # PLACEHOLDER: Stellar parameters, radial velocities
        return None


class GALAHFetcher(BaseFetcher):
    """GALAH Survey (high-resolution stellar spectroscopy)."""
    
    def __init__(self):
        super().__init__('galah')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch GALAH data. TODO: Implement."""
        # PLACEHOLDER: Abundances, stellar parameters
        return None


class APOGEEFetcher(BaseFetcher):
    """APOGEE Survey (infrared stellar spectroscopy)."""
    
    def __init__(self):
        super().__init__('apogee')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch APOGEE data. TODO: Implement."""
        # PLACEHOLDER: H-band spectra, abundances
        return None


class RAVEFetcher(BaseFetcher):
    """RAVE Survey (radial velocity experiment)."""
    
    def __init__(self):
        super().__init__('rave')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch RAVE data. TODO: Implement."""
        # PLACEHOLDER: Radial velocities, stellar parameters
        return None


class DESISpectraFetcher(BaseFetcher):
    """DESI Spectral Archive."""
    
    def __init__(self):
        super().__init__('desi_spectra')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch DESI spectra. TODO: Implement."""
        # PLACEHOLDER: Galaxy/quasar spectra
        return None


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED FETCHER (All 50 APIs)
# ═══════════════════════════════════════════════════════════════════════════════

class UnifiedFetcher:
    """
    Unified API fetcher that tries all 50 sources in priority order.
    
    Priority Groups:
        A. Astronomical Databases (1-5): SIMBAD, NED, VizieR, Gaia, MAST
        B. NASA/Space Agencies (6-10): Exoplanet, HEASARC, ADS, JPL, ESO
        C. Sky Surveys (11-15): SDSS, 2MASS, WISE, Pan-STARRS, ZTF
        D. Specialized (16-20): ATNF Pulsar, McGill Magnetar, TNS, GCN, GWOSC
        E. Computational/AI (21-25): Wolfram, arXiv, Grok, OpenAI, Claude
        F. Radio/Infrared (26-30): NVSS, FIRST, VLASS, ALMA, ASKAP
        G. X-ray/Gamma-ray (31-35): Chandra, XMM, Swift, Fermi, INTEGRAL
        H. Space Telescopes (36-40): HST, JWST, Spitzer, Kepler, TESS
        I. Cosmology/CMB (41-45): Planck, WMAP, DES, DESI, Euclid
        J. Spectroscopic (46-50): LAMOST, GALAH, APOGEE, RAVE, DESI Spectra
    
    Results are stored in IPData.py for QCalc.py consumption.
    """
    
    def __init__(self):
        # Primary databases (always try first)
        self.simbad = SIMBADFetcher()
        self.ned = NEDFetcher()
        self.grok = GrokFetcher()
        
        # All 50 fetchers (will be instantiated on demand)
        self._fetchers = {
            # Group A: Astronomical Databases (1-5)
            'simbad': self.simbad,
            'ned': self.ned,
            'vizier': None,
            'gaia': None,
            'mast': None,
            # Group B: NASA/Space Agencies (6-10)
            'nasa_exoplanet': None,
            'heasarc': None,
            'ads': None,
            'jpl_horizons': None,
            'eso': None,
            # Group C: Sky Surveys (11-15)
            'sdss': None,
            '2mass': None,
            'wise': None,
            'panstarrs': None,
            'ztf': None,
            # Group D: Specialized (16-20)
            'atnf_pulsar': None,
            'mcgill_magnetar': None,
            'tns': None,
            'gcn': None,
            'gwosc': None,
            # Group E: Computational/AI (21-25)
            'wolfram': None,
            'arxiv': None,
            'grok': self.grok,
            'openai': None,
            'claude': None,
            # Group F: Radio/Infrared (26-30)
            'nvss': None,
            'first': None,
            'vlass': None,
            'alma': None,
            'askap': None,
            # Group G: X-ray/Gamma-ray (31-35)
            'chandra': None,
            'xmm': None,
            'swift': None,
            'fermi': None,
            'integral': None,
            # Group H: Space Telescopes (36-40)
            'hst': None,
            'jwst': None,
            'spitzer': None,
            'kepler': None,
            'tess': None,
            # Group I: Cosmology/CMB (41-45)
            'planck': None,
            'wmap': None,
            'des': None,
            'desi': None,
            'euclid': None,
            # Group J: Spectroscopic (46-50)
            'lamost': None,
            'galah': None,
            'apogee': None,
            'rave': None,
            'desi_spectra': None,
        }
    
    def fetch(self, object_name: str, required_params: List[str] = None) -> Dict[str, Any]:
        """
        Fetch parameters for an astronomical object from all available sources.
        
        Args:
            object_name: Name of the object to query
            required_params: List of required parameter names
            
        Returns:
            Merged dictionary with all available parameters
        """
        if required_params is None:
            required_params = ['mass', 'distance', 'temperature']
        
        result = {'name': object_name, 'sources': []}
        
        # Try SIMBAD first
        simbad_data = self.simbad.fetch(object_name)
        if simbad_data:
            result.update({k: v for k, v in simbad_data.items() if k != 'source'})
            result['sources'].append('SIMBAD')
        
        # Try NED for extragalactic objects
        ned_data = self.ned.fetch(object_name)
        if ned_data:
            for key, value in ned_data.items():
                if key not in result and key != 'source':
                    result[key] = value
            result['sources'].append('NED')
        
        # Check for missing required parameters
        missing = [p for p in required_params if p not in result]
        
        # Use Grok fallback for missing parameters
        if missing:
            grok_data = self.grok.fetch(object_name, missing)
            if grok_data:
                for key, value in grok_data.items():
                    if key not in result and key != 'source':
                        result[key] = value
                result['sources'].append('Grok')
        
        return result
    
    def fetch_to_ipdata(self, object_name: str, required_params: List[str] = None) -> Tuple[InputParameters, str]:
        """
        Fetch parameters and store in IPData.py.
        
        This is the primary method for QCalc.py integration.
        
        Args:
            object_name: Name of the object to query
            required_params: List of required parameter names
            
        Returns:
            Tuple of (InputParameters, query_id)
        """
        # Fetch from APIs
        raw_result = self.fetch(object_name, required_params)
        
        # Convert to InputParameters
        params = InputParameters(
            query_name=object_name,
            sources=raw_result.get('sources', []),
            # Map fetched values to InputParameters fields
            M=raw_result.get('mass'),
            d=raw_result.get('distance'),
            R=raw_result.get('radius'),
            T=raw_result.get('temperature'),
            L=raw_result.get('luminosity'),
            B=raw_result.get('magnetic_field'),
            z=raw_result.get('redshift'),
            v_disp=raw_result.get('velocity_dispersion'),
            v_rad=raw_result.get('radial_velocity'),
            parallax=raw_result.get('parallax'),
            spectral_type=raw_result.get('spectral_type'),
        )
        
        # Store in IPData
        query_id = store_input(params)
        
        return params, query_id
    
    def fetch_and_save(self, object_name: str, output_dir: str = ".") -> Tuple[Dict[str, Any], str]:
        """
        Fetch parameters and save to timestamped CSV file.
        Also stores in IPData.py.
        
        Args:
            object_name: Name of the object to query
            output_dir: Directory for output CSV
            
        Returns:
            Tuple of (result_dict, csv_filepath)
        """
        result = self.fetch(object_name)
        
        # Also store in IPData
        self.fetch_to_ipdata(object_name)
        
        # Generate timestamped filename
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"bodies_{timestamp}.csv"
        filepath = os.path.join(output_dir, filename)
        
        # Save to CSV
        with open(filepath, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=result.keys())
            writer.writeheader()
            writer.writerow(result)
        
        return result, filepath
    
    def list_available_apis(self) -> Dict[str, bool]:
        """List all 25 APIs and their implementation status."""
        return {
            name: (fetcher is not None and hasattr(fetcher, 'fetch'))
            for name, fetcher in self._fetchers.items()
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MANUAL INPUT HELPER
# ═══════════════════════════════════════════════════════════════════════════════

def manual_input() -> Dict[str, Any]:
    """
    Interactive manual input for astronomical parameters.
    
    Returns:
        Dictionary with user-provided parameters
    """
    print("\n" + "=" * 60)
    print("Manual Parameter Input")
    print("=" * 60)
    print("Enter values in SI units (leave blank to skip)")
    print("-" * 60)
    
    params = {}
    
    params['name'] = input("Object name: ").strip() or "manual_input"
    
    # Core parameters
    m = input("Mass (kg, or 'M_sun' for solar masses): ").strip()
    if m:
        if 'M_sun' in m:
            try:
                factor = float(m.replace('M_sun', '').replace('*', '').strip() or '1')
                params['mass'] = factor * UNITS['M_sun']
            except:
                params['mass'] = float(m)
        else:
            params['mass'] = float(m)
    
    r = input("Distance (m, or 'pc/kpc/Mpc' suffix): ").strip()
    if r:
        for unit in ['Mpc', 'kpc', 'pc', 'ly', 'AU']:
            if unit in r:
                try:
                    factor = float(r.replace(unit, '').strip())
                    params['distance'] = factor * UNITS[unit]
                    break
                except:
                    pass
        else:
            params['distance'] = float(r)
    
    T = input("Temperature (K): ").strip()
    if T:
        params['temperature'] = float(T)
    
    L = input("Luminosity (W, or 'L_sun' for solar luminosities): ").strip()
    if L:
        if 'L_sun' in L:
            try:
                factor = float(L.replace('L_sun', '').replace('*', '').strip() or '1')
                params['luminosity'] = factor * UNITS['L_sun']
            except:
                params['luminosity'] = float(L)
        else:
            params['luminosity'] = float(L)
    
    z = input("Redshift (dimensionless): ").strip()
    if z:
        params['redshift'] = float(z)
    
    B = input("Magnetic field (T): ").strip()
    if B:
        params['magnetic_field'] = float(B)
    
    params['source'] = 'manual'
    
    return params


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL FETCHER INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

FETCHER = UnifiedFetcher()


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def fetch(object_name: str) -> Dict[str, Any]:
    """Fetch parameters for an object from all sources."""
    return FETCHER.fetch(object_name)


def fetch_and_save(object_name: str, output_dir: str = ".") -> Tuple[Dict[str, Any], str]:
    """Fetch parameters and save to CSV."""
    return FETCHER.fetch_and_save(object_name, output_dir)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 80)
    print("APIFetch.py - API Data Fetching Layer Test")
    print("=" * 80)
    
    # Test with a known object
    test_object = "Betelgeuse"
    print(f"\nFetching data for: {test_object}")
    print("-" * 40)
    
    result = fetch(test_object)
    
    print("Results:")
    for key, value in result.items():
        if isinstance(value, float) and abs(value) > 1e6:
            print(f"  {key}: {value:.4e}")
        else:
            print(f"  {key}: {value}")
    
    print("\n" + "-" * 40)
    print("Sources used:", result.get('sources', []))
