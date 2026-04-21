#!/usr/bin/env python3
"""
APIFetch.py - UQFF API Data Fetching Layer (55 APIs)
=====================================================

Fetches astronomical parameters from external APIs for UQFF calculations.

ARCHITECTURE:
    User Query → APIFetch.py (this file) → IPData.py → QCalc.py → OPData.py
    
55 SUPPORTED APIs (Updated with NASA Keys):
    ─────────────────────────────────────────────────────────────────────────
    ASTRONOMICAL DATABASES (Primary)
    ─────────────────────────────────────────────────────────────────────────
     1. SIMBAD       - CDS Strasbourg (~15M objects)
     2. NED          - NASA/IPAC Extragalactic Database (~300M objects)
     3. VizieR       - CDS catalog service (20,000+ catalogs)
     4. Gaia         - ESA astrometry mission (1.8B stars)
     5. MAST         - Mikulski Archive (HST, JWST, etc.)
    ─────────────────────────────────────────────────────────────────────────
    NASA/SPACE AGENCIES (EXPANDED WITH API KEYS)
    ─────────────────────────────────────────────────────────────────────────
     6. NASA APOD    - Astronomy Picture of the Day [NEW - IMPLEMENTED]
     7. NASA NeoWs   - Near Earth Object Web Service [NEW - IMPLEMENTED]
     8. NASA Mars    - Mars Weather (InSight) [NEW]
     9. NASA EPIC    - Earth Polychromatic Imaging Camera [NEW]
    10. NASA DONKI   - Space Weather Database (CME, Solar Flares) [NEW - IMPLEMENTED]
    11. NASA Exoplanet Archive - Confirmed exoplanets
    12. NASA HEASARC - High Energy Astrophysics data
    13. NASA ADS     - Astrophysics Data System (literature)
    14. NASA JPL Horizons - Solar system ephemerides
    15. ESO Archive  - European Southern Observatory
    ─────────────────────────────────────────────────────────────────────────
    SURVEY DATABASES
    ─────────────────────────────────────────────────────────────────────────
    16. SDSS         - Sloan Digital Sky Survey
    17. 2MASS        - Two Micron All Sky Survey
    18. WISE         - Wide-field Infrared Survey Explorer
    19. Pan-STARRS   - Panoramic Survey Telescope
    20. ZTF          - Zwicky Transient Facility
    ─────────────────────────────────────────────────────────────────────────
    SPECIALIZED DATABASES
    ─────────────────────────────────────────────────────────────────────────
    21. ATNF         - Pulsar catalog
    22. McGill       - Magnetar catalog
    23. TNS          - Transient Name Server (supernovae/transients)
    24. GCN          - Gamma-ray Coordinates Network
    25. LIGO/Virgo   - Gravitational wave events
    ─────────────────────────────────────────────────────────────────────────
    COMPUTATIONAL/AI
    ─────────────────────────────────────────────────────────────────────────
    26. Wolfram Alpha - Computational knowledge engine
    27. arXiv        - Preprint server (metadata)
    28. Grok (xAI)   - AI fallback for missing data
    29. OpenAI       - Alternative AI fallback
    30. Claude       - Alternative AI fallback
    ─────────────────────────────────────────────────────────────────────────
    RADIO/INFRARED SURVEYS
    ─────────────────────────────────────────────────────────────────────────
    31. NVSS         - NRAO VLA Sky Survey
    32. FIRST        - Faint Images Radio Sky at Twenty-cm
    33. VLASS        - VLA Sky Survey
    34. ALMA         - Atacama Large Millimeter Array Archive
    35. ASKAP        - Australian Square Kilometre Array Archive
    ─────────────────────────────────────────────────────────────────────────
    X-RAY/GAMMA-RAY
    ─────────────────────────────────────────────────────────────────────────
    36. Chandra      - Chandra X-ray Observatory
    37. XMM-Newton   - XMM-Newton Archive
    38. Swift        - Swift Gamma-Ray Burst Mission
    39. Fermi        - Fermi Gamma-ray Space Telescope
    40. INTEGRAL     - INTEGRAL Gamma-ray Observatory
    ─────────────────────────────────────────────────────────────────────────
    SPACE TELESCOPES
    ─────────────────────────────────────────────────────────────────────────
    41. HST          - Hubble Space Telescope
    42. JWST         - James Webb Space Telescope
    43. Spitzer      - Spitzer Space Telescope (archived)
    44. Kepler       - Kepler Mission
    45. TESS         - Transiting Exoplanet Survey Satellite
    ─────────────────────────────────────────────────────────────────────────
    COSMOLOGY/CMB
    ─────────────────────────────────────────────────────────────────────────
    46. Planck       - Planck CMB mission
    47. WMAP         - Wilkinson Microwave Anisotropy Probe (archived)
    48. DES          - Dark Energy Survey
    49. DESI         - Dark Energy Spectroscopic Instrument
    50. Euclid       - ESA Euclid mission
    ─────────────────────────────────────────────────────────────────────────
    SPECTROSCOPIC SURVEYS
    ─────────────────────────────────────────────────────────────────────────
    51. LAMOST       - Large Sky Area Multi-Object Fiber Spec Telescope
    52. GALAH        - Galactic Archaeology with HERMES
    53. APOGEE       - Apache Point Observatory Galactic Evolution Exp
    54. RAVE         - RAdial Velocity Experiment
    55. DESI Spectra - DESI Spectroscopic Data
    ─────────────────────────────────────────────────────────────────────────

API KEYS CONFIGURED:
    - NASA_API_KEY_1: [PROMPT_FOR_NASA_API_KEY_1]
    - NASA_API_KEY_2: [PROMPT_FOR_NASA_API_KEY_2]
    - MAST_API_KEY: [PROMPT_FOR_MAST_API_KEY]
    - XAI_API_KEY (Grok): Set via environment variable XAI_API_KEY

OUTPUT: 
    All results written to IPData.py (InputParameters dataclass)
    
Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import requests
import argparse
import json
import os
import csv
import sys
from datetime import datetime, timedelta
from typing import Dict, List, Any, Optional, Tuple
import re
import xml.etree.ElementTree as ET

# Import IPData for storing results
from IPData import InputParameters, InputDataStore, INPUT_STORE, store_input

# ═══════════════════════════════════════════════════════════════════════════════
# API KEYS (Store as environment variables for security, fallback to these)
# ═══════════════════════════════════════════════════════════════════════════════

API_KEYS = {
    'NASA_API_KEY_1': os.environ.get('NASA_API_KEY_1', '[PROMPT_FOR_NASA_API_KEY_1]'),
    'NASA_API_KEY_2': os.environ.get('NASA_API_KEY_2', '[PROMPT_FOR_NASA_API_KEY_2]'),
    'MAST_API_KEY': os.environ.get('MAST_API_KEY', '[PROMPT_FOR_MAST_API_KEY]'),
    'XAI_API_KEY': os.environ.get('XAI_API_KEY'),  # Grok AI - Set environment variable
    'OPENAI_API_KEY': os.environ.get('OPENAI_API_KEY', ''),
    'ANTHROPIC_API_KEY': os.environ.get('ANTHROPIC_API_KEY', ''),
    'WOLFRAM_APP_ID': os.environ.get('WOLFRAM_APP_ID', ''),
}

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
    # GROUP B: NASA/SPACE AGENCIES (6-15)
    # ═══════════════════════════════════════════════════════════════════════════
    'nasa_apod': f'https://api.nasa.gov/planetary/apod',
    'nasa_neows': f'https://api.nasa.gov/neo/rest/v1/feed',
    'nasa_mars_weather': f'https://api.nasa.gov/insight_weather/',
    'nasa_epic': f'https://api.nasa.gov/EPIC/api/natural',
    'nasa_donki': f'https://api.nasa.gov/DONKI/CME',
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
        # Sanitize object name for ADQL (escape single quotes)
        safe_name = object_name.replace("'", "''")
        # TAP query using ident JOIN (correct SIMBAD schema)
        query = f"""
        SELECT TOP 1
            b.main_id, b.ra, b.dec, b.pmra, b.pmdec, b.plx_value, b.rvz_radvel,
            b.sp_type
        FROM basic AS b
        JOIN ident AS i ON b.oid = i.oidref
        WHERE i.id = '{safe_name}'
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
            print(f"SIMBAD fetch error: {e}", file=sys.stderr)
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
    
    # NED preferred name aliases for common short names
    _NED_ALIASES = {
        'M87': 'Messier 087', 'M31': 'Messier 031', 'M33': 'Messier 033',
        'M81': 'Messier 081', 'M82': 'Messier 082', 'M51': 'Messier 051',
        'M104': 'Messier 104', 'M101': 'Messier 101', 'M83': 'Messier 083',
    }

    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """
        Fetch data for an extragalactic object from NED.
        
        Args:
            object_name: Name of the object (e.g., "M87", "NGC 1365")
            
        Returns:
            Dictionary with available parameters or None if not found
        """
        # Sanitize for ADQL (escape single quotes)
        safe_name = object_name.replace("'", "''")
        # Try NED preferred name alias first, then original
        ned_name = self._NED_ALIASES.get(object_name, safe_name)
        query = f"""
        SELECT TOP 1
            prefname, ra, dec, z
        FROM objdir
        WHERE prefname = '{ned_name}'
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
                timeout=45
            )
            
            if response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    return self._parse_response(data['data'][0], object_name)
            
            return None
            
        except Exception as e:
            print(f"NED fetch error: {e}", file=sys.stderr)
            return None
    
    def _parse_response(self, row: List, object_name: str) -> Dict[str, Any]:
        """Parse NED response into standardized format."""
        result = {'name': object_name, 'source': 'NED'}
        
        # RA, Dec (indices 1, 2)
        if len(row) > 2 and row[1] is not None:
            result['ra'] = row[1]
            result['dec'] = row[2]
        
        # Redshift (index 3)
        if len(row) > 3 and row[3] is not None:
            result['redshift'] = row[3]
        
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
            print("Warning: XAI_API_KEY not set, Grok fallback unavailable", file=sys.stderr)
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
                    'model': 'grok-4',
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
            print(f"Grok fetch error: {e}", file=sys.stderr)
            return None


# ═══════════════════════════════════════════════════════════════════════════════
# NASA APOD AND SERVICES FETCHERS
# ═══════════════════════════════════════════════════════════════════════════════

class NASAAPODFetcher:
    """
    NASA Astronomy Picture of the Day (APOD) API.
    
    Retrieves daily astronomy images with descriptions and metadata.
    Useful for contextual information about astronomical objects.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['nasa_apod']
        self.api_key = API_KEYS['NASA_API_KEY_1']
        self.timeout = 30
    
    def fetch(self, keywords: str = None, date: str = None) -> Optional[Dict[str, Any]]:
        """
        Fetch APOD data.
        
        Args:
            keywords: Search keywords for concept tags
            date: Date in YYYY-MM-DD format (defaults to today)
            
        Returns:
            Dictionary with image URL, explanation, title, and metadata
        """
        params = {
            'api_key': self.api_key,
            'concept_tags': True
        }
        
        if keywords:
            params['keywords'] = keywords
        if date:
            params['date'] = date
        
        try:
            response = requests.get(
                self.endpoint,
                params=params,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                return {
                    'name': data.get('title', 'APOD'),
                    'source': 'NASA_APOD',
                    'description': data.get('explanation', ''),
                    'image_url': data.get('url', ''),
                    'hd_image_url': data.get('hdurl', ''),
                    'date': data.get('date', ''),
                    'media_type': data.get('media_type', ''),
                    'copyright': data.get('copyright', 'NASA')
                }
            
            return None
            
        except Exception as e:
            print(f"NASA APOD fetch error: {e}", file=sys.stderr)
            return None


class NASANeoWsFetcher:
    """
    NASA Near Earth Object Web Service (NeoWs).
    
    Provides data on near-Earth asteroids and comets.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['nasa_neows']
        self.api_key = API_KEYS['NASA_API_KEY_1']
        self.timeout = 30
    
    def fetch(self, start_date: str = None, end_date: str = None) -> Optional[Dict[str, Any]]:
        """
        Fetch Near Earth Objects.
        
        Args:
            start_date: Start date in YYYY-MM-DD format
            end_date: End date in YYYY-MM-DD format
            
        Returns:
            Dictionary with NEO data including orbital parameters
        """
        if not start_date:
            start_date = datetime.now().strftime('%Y-%m-%d')
        if not end_date:
            end_date = start_date
        
        params = {
            'start_date': start_date,
            'end_date': end_date,
            'api_key': self.api_key
        }
        
        try:
            response = requests.get(
                self.endpoint,
                params=params,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                neo_list = []
                
                for date_key in data.get('near_earth_objects', {}):
                    for neo in data['near_earth_objects'][date_key]:
                        neo_list.append({
                            'name': neo.get('name', ''),
                            'source': 'NASA_NeoWs',
                            'diameter_min': neo.get('estimated_diameter', {}).get('meters', {}).get('estimated_diameter_min'),
                            'diameter_max': neo.get('estimated_diameter', {}).get('meters', {}).get('estimated_diameter_max'),
                            'is_hazardous': neo.get('is_potentially_hazardous_asteroid', False),
                            'close_approach_date': neo.get('close_approach_data', [{}])[0].get('close_approach_date'),
                            'velocity': float(neo.get('close_approach_data', [{}])[0].get('relative_velocity', {}).get('kilometers_per_second', 0)) * 1000,  # m/s
                            'miss_distance': float(neo.get('close_approach_data', [{}])[0].get('miss_distance', {}).get('kilometers', 0)) * 1000,  # meters
                        })
                
                return {'objects': neo_list, 'count': len(neo_list)}
            
            return None
            
        except Exception as e:
            print(f"NASA NeoWs fetch error: {e}", file=sys.stderr)
            return None


class NASADONKIFetcher:
    """
    NASA Space Weather Database Of Notifications, Knowledge, Information (DONKI).
    
    Provides space weather events: CMEs, solar flares, geomagnetic storms.
    """
    
    def __init__(self):
        self.endpoint = ENDPOINTS['nasa_donki']
        self.api_key = API_KEYS['NASA_API_KEY_1']
        self.timeout = 30
    
    def fetch(self, start_date: str = None, end_date: str = None, event_type: str = 'CME') -> Optional[Dict[str, Any]]:
        """
        Fetch space weather events.
        
        Args:
            start_date: Start date YYYY-MM-DD
            end_date: End date YYYY-MM-DD
            event_type: 'CME', 'FLR', 'SEP', 'MPC', 'GST', 'RBE', 'HSS'
            
        Returns:
            Dictionary with space weather event data
        """
        if not start_date:
            start_date = (datetime.now() - timedelta(days=30)).strftime('%Y-%m-%d')
        if not end_date:
            end_date = datetime.now().strftime('%Y-%m-%d')
        
        params = {
            'startDate': start_date,
            'endDate': end_date,
            'api_key': self.api_key
        }
        
        try:
            response = requests.get(
                self.endpoint.replace('CME', event_type),
                params=params,
                timeout=self.timeout
            )
            
            if response.status_code == 200:
                data = response.json()
                return {
                    'source': 'NASA_DONKI',
                    'event_type': event_type,
                    'events': data,
                    'count': len(data) if isinstance(data, list) else 1
                }
            
            return None
            
        except Exception as e:
            print(f"NASA DONKI fetch error: {e}", file=sys.stderr)
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

    def _log_success(self):
        API_STATUS[self.api_name]['last_call'] = datetime.now().isoformat()

    def _safe_float(self, value: Any) -> Optional[float]:
        try:
            if value is None or value == '':
                return None
            return float(value)
        except (TypeError, ValueError):
            return None

    def _null_sentinel(self, value: Any, sentinels: Optional[List[float]] = None) -> Optional[float]:
        numeric = self._safe_float(value)
        if numeric is None:
            return None
        for sentinel in sentinels or [-999.0, -9999.0]:
            if abs(numeric - sentinel) < 1e-9:
                return None
        return numeric

    def _parse_votable_first_row(self, xml_text: str) -> Optional[Dict[str, Any]]:
        try:
            root = ET.fromstring(xml_text)
        except ET.ParseError:
            return None

        namespace = {'v': 'http://www.ivoa.net/xml/VOTable/v1.3'}
        fields = [field.attrib.get('name') for field in root.findall('.//v:FIELD', namespace)]
        row = root.find('.//v:TR', namespace)
        if not fields or row is None:
            return None

        values = [cell.text for cell in row.findall('v:TD', namespace)]
        if not values:
            return None

        return {
            key: values[index] if index < len(values) else None
            for index, key in enumerate(fields)
            if key
        }

    def _extract_match_count(self, payload: Any) -> int:
        if isinstance(payload, list):
            return len(payload)
        if isinstance(payload, dict):
            for key in ('totalRows', 'count', 'total', 'recordsTotal'):
                value = payload.get(key)
                numeric = self._safe_float(value)
                if numeric is not None:
                    return int(numeric)
            data = payload.get('data') or payload.get('results') or payload.get('observations')
            if isinstance(data, list):
                return len(data)
        if isinstance(payload, str):
            count_match = re.search(r'([0-9]+)\s+(?:results|observations|matches|rows)', payload, re.IGNORECASE)
            if count_match:
                return int(count_match.group(1))
        return 0

    def _build_archive_result(
        self,
        object_name: str,
        source: str,
        archive_url: str,
        description: str,
        payload: Any,
        extra: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        result = {
            'name': object_name,
            'source': source,
            'archive_url': archive_url,
            'description': description,
            'match_count': self._extract_match_count(payload),
        }
        if extra:
            result.update({key: value for key, value in extra.items() if value is not None})
        return result
    
    def _log_error(self, e: Exception):
        API_STATUS[self.api_name]['errors'] += 1
        print(f"{self.api_name} fetch error: {e}", file=sys.stderr)

    def _tap_sync_query(self, query: str, endpoint: Optional[str] = None) -> Optional[List[Any]]:
        response = requests.post(
            endpoint or self.endpoint,
            data={
                'REQUEST': 'doQuery',
                'LANG': 'ADQL',
                'FORMAT': 'json',
                'QUERY': query,
            },
            timeout=self.timeout,
        )
        response.raise_for_status()
        payload = response.json()
        rows = payload.get('data') if isinstance(payload, dict) else None
        return rows if isinstance(rows, list) else None

    def _resolve_object_coordinates(self, object_name: str) -> Optional[Tuple[float, float]]:
        simbad_data = SIMBADFetcher().fetch(object_name)
        if simbad_data and simbad_data.get('ra') is not None and simbad_data.get('dec') is not None:
            return float(simbad_data['ra']), float(simbad_data['dec'])

        ned_data = NEDFetcher().fetch(object_name)
        if ned_data and ned_data.get('ra') is not None and ned_data.get('dec') is not None:
            return float(ned_data['ra']), float(ned_data['dec'])

        return None


# ─── API 3: VizieR ──────────────────────────────────────────────────────────────

class VizieRFetcher(BaseFetcher):
    """VizieR catalog service (20,000+ astronomical catalogs)."""
    
    def __init__(self):
        super().__init__('vizier')
    
    def fetch(self, object_name: str, catalog: str = None) -> Optional[Dict[str, Any]]:
        """Fetch Pan-STARRS DR1 mean magnitudes from VizieR for the resolved target."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        params = {
            '-source': catalog or 'II/349/ps1',
            '-c': f'{ra} {dec}',
            '-c.rs': 0.02,
            '-out.max': 1,
        }

        try:
            response = requests.get(self.endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            row = self._parse_votable_first_row(response.text)
            if not row:
                return None

            result = self._build_archive_result(
                object_name,
                'VizieR',
                response.url,
                'VizieR catalog summary using the Pan-STARRS DR1 reference table.',
                response.text,
                {
                    'catalog_id': row.get('objID'),
                    'ra': self._safe_float(row.get('RAJ2000')),
                    'dec': self._safe_float(row.get('DEJ2000')),
                    'g_mag': self._null_sentinel(row.get('gmag')),
                    'r_mag': self._null_sentinel(row.get('rmag')),
                    'i_mag': self._null_sentinel(row.get('imag')),
                    'z_mag': self._null_sentinel(row.get('zmag')),
                    'y_mag': self._null_sentinel(row.get('ymag')),
                },
            )
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 4: Gaia ────────────────────────────────────────────────────────────────

class GaiaFetcher(BaseFetcher):
    """ESA Gaia mission (1.8 billion stars with astrometry)."""
    
    def __init__(self):
        super().__init__('gaia')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Gaia DR3 astrometry around a resolved target position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        query = f"""
        SELECT TOP 1
            source_id, ra, dec, parallax, pmra, pmdec, radial_velocity,
            teff_gspphot, luminosity_gspphot
        FROM gaiadr3.gaia_source
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, 0.01)
        )
        ORDER BY parallax DESC
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None

            row = rows[0]
            parallax = self._safe_float(row[3]) if len(row) > 3 else None
            radial_velocity = self._safe_float(row[6]) if len(row) > 6 else None
            luminosity = self._safe_float(row[8]) if len(row) > 8 else None
            result = self._build_archive_result(
                object_name,
                'Gaia',
                self.endpoint,
                'Gaia DR3 astrometric and stellar-parameter match.',
                rows,
                {
                    'source_id': row[0] if len(row) > 0 else None,
                    'ra': row[1] if len(row) > 1 else None,
                    'dec': row[2] if len(row) > 2 else None,
                    'parallax': parallax,
                    'proper_motion_ra': row[4] if len(row) > 4 else None,
                    'proper_motion_dec': row[5] if len(row) > 5 else None,
                    'radial_velocity': radial_velocity * 1000 if radial_velocity is not None else None,
                    'temperature': row[7] if len(row) > 7 else None,
                    'luminosity': luminosity * UNITS['L_sun'] if luminosity is not None else None,
                },
            )
            if parallax and parallax > 0:
                result['distance'] = (1000.0 / parallax) * UNITS['pc']
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 5: MAST ────────────────────────────────────────────────────────────────

class MASTFetcher(BaseFetcher):
    """Mikulski Archive for Space Telescopes (HST, JWST, TESS, Kepler)."""
    
    def __init__(self):
        super().__init__('mast')
        self.api_key = API_KEYS['MAST_API_KEY']
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Resolve a target in MAST and return archive-routing metadata."""
        request_payload = {
            'service': 'Mast.Name.Lookup',
            'params': {
                'input': object_name,
                'format': 'json',
            },
        }
        headers = {'Content-Type': 'application/x-www-form-urlencoded'}
        if self.api_key:
            headers['Authorization'] = f'token {self.api_key}'

        try:
            response = requests.post(
                self.endpoint,
                data={'request': json.dumps(request_payload)},
                headers=headers,
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
            rows = payload.get('resolvedCoordinate') or payload.get('data') or payload.get('resolved') or []
            if not rows:
                return None

            row = rows[0]
            result = self._build_archive_result(
                object_name,
                'MAST',
                response.url,
                'MAST name resolver result for HST/JWST/TESS/Kepler archive routing.',
                rows,
                {
                    'target_name': row.get('canonicalName') or row.get('objectname') or object_name,
                    'ra': self._safe_float(row.get('ra')),
                    'dec': self._safe_float(row.get('decl') if 'decl' in row else row.get('dec')),
                },
            )
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 6: NASA Exoplanet Archive ──────────────────────────────────────────────

class ExoplanetFetcher(BaseFetcher):
    """NASA Exoplanet Archive (confirmed exoplanets and host stars)."""
    
    def __init__(self):
        super().__init__('nasa_exoplanet')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch confirmed exoplanet or host-star properties from pscomppars."""
        safe_name = object_name.replace("'", "''")
        query = f"""
        SELECT TOP 1
            pl_name, hostname, pl_bmasse, pl_rade, pl_orbper,
            st_mass, st_rad, st_teff, sy_dist
        FROM pscomppars
        WHERE lower(pl_name) = lower('{safe_name}') OR lower(hostname) = lower('{safe_name}')
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None

            row = rows[0]
            star_mass = self._safe_float(row[5]) if len(row) > 5 else None
            star_radius = self._safe_float(row[6]) if len(row) > 6 else None
            distance_pc = self._safe_float(row[8]) if len(row) > 8 else None
            result = {
                'name': object_name,
                'source': 'NASA_Exoplanet',
                'planet_name': row[0] if len(row) > 0 else None,
                'host_name': row[1] if len(row) > 1 else None,
                'planet_mass': self._safe_float(row[2]) * UNITS['M_earth'] if len(row) > 2 and self._safe_float(row[2]) is not None else None,
                'planet_radius': self._safe_float(row[3]) * UNITS['R_earth'] if len(row) > 3 and self._safe_float(row[3]) is not None else None,
                'orbital_period_days': row[4] if len(row) > 4 else None,
                'mass': star_mass * UNITS['M_sun'] if star_mass is not None else None,
                'radius': star_radius * UNITS['R_sun'] if star_radius is not None else None,
                'temperature': row[7] if len(row) > 7 else None,
                'distance': distance_pc * UNITS['pc'] if distance_pc is not None else None,
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 7: HEASARC ─────────────────────────────────────────────────────────────

class HEASARCFetcher(BaseFetcher):
    """NASA High Energy Astrophysics Science Archive (X-ray, gamma-ray)."""
    
    def __init__(self):
        super().__init__('heasarc')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch HEASARC master-catalog summary around the resolved source position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        query = f"""
        SELECT TOP 1
            name, ra, dec, flux, class
        FROM heasarc_master
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, 0.05)
        )
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None
            row = rows[0]
            result = self._build_archive_result(
                object_name,
                'HEASARC',
                self.endpoint,
                'HEASARC high-energy archive summary near the resolved target.',
                rows,
                {
                    'target_name': row[0] if len(row) > 0 else object_name,
                    'ra': row[1] if len(row) > 1 else None,
                    'dec': row[2] if len(row) > 2 else None,
                    'xray_flux': row[3] if len(row) > 3 else None,
                    'high_energy_class': row[4] if len(row) > 4 else None,
                },
            )
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
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
        """Fetch recent ESO archive metadata for a target name."""
        safe_name = object_name.replace("'", "''")
        query = f"""
        SELECT TOP 1
            target_name, dp_id, instrument_name, s_ra, s_dec, t_min, em_min, em_max
        FROM ivoa.ObsCore
        WHERE lower(target_name) LIKE lower('%{safe_name}%')
        ORDER BY t_min DESC
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None
            row = rows[0]
            result = self._build_archive_result(
                object_name,
                'ESO',
                self.endpoint,
                'ESO archive observation summary for spectroscopy/imaging follow-up.',
                rows,
                {
                    'target_name': row[0] if len(row) > 0 else object_name,
                    'dataset_id': row[1] if len(row) > 1 else None,
                    'instrument': row[2] if len(row) > 2 else None,
                    'ra': row[3] if len(row) > 3 else None,
                    'dec': row[4] if len(row) > 4 else None,
                    'observation_mjd': row[5] if len(row) > 5 else None,
                    'wavelength_min_m': row[6] if len(row) > 6 else None,
                    'wavelength_max_m': row[7] if len(row) > 7 else None,
                },
            )
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 11: SDSS ───────────────────────────────────────────────────────────────

class SDSSFetcher(BaseFetcher):
    """Sloan Digital Sky Survey (photometry, spectroscopy)."""
    
    def __init__(self):
        super().__init__('sdss')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch SDSS photometry and spectroscopic redshift near the resolved target position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        query = f"""
        SELECT TOP 1
            p.objid, p.ra, p.dec, p.u, p.g, p.r, p.i, p.z,
            s.z AS specz
        FROM PhotoObjAll AS p
        LEFT JOIN SpecObjAll AS s ON p.objID = s.bestObjID
        WHERE p.ra BETWEEN {ra - 0.02} AND {ra + 0.02}
          AND p.dec BETWEEN {dec - 0.02} AND {dec + 0.02}
        ORDER BY ABS(p.ra - {ra}) + ABS(p.dec - {dec})
        """

        try:
            response = requests.get(
                self.endpoint,
                params={
                    'cmd': query,
                    'format': 'json',
                },
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
            rows = []
            if isinstance(payload, list):
                for table in payload:
                    table_rows = table.get('Rows') if isinstance(table, dict) else None
                    if table_rows:
                        rows = table_rows
                        break
            elif isinstance(payload, dict):
                rows = payload.get('Rows', [])
            if not rows:
                return None

            row = rows[0]
            result = {
                'name': object_name,
                'source': 'SDSS',
                'catalog_id': row.get('objid') or row.get('objID'),
                'ra': self._safe_float(row.get('ra')),
                'dec': self._safe_float(row.get('dec')),
                'u_mag': self._null_sentinel(row.get('u')),
                'g_mag': self._null_sentinel(row.get('g')),
                'r_mag': self._null_sentinel(row.get('r')),
                'i_mag': self._null_sentinel(row.get('i')),
                'z_mag': self._null_sentinel(row.get('z')),
                'redshift': self._safe_float(row.get('specz')),
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 12: 2MASS ──────────────────────────────────────────────────────────────

class TwoMASSFetcher(BaseFetcher):
    """Two Micron All Sky Survey (near-infrared JHK photometry)."""
    
    def __init__(self):
        super().__init__('2mass')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch near-infrared photometry from 2MASS around the target position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        query = f"""
        SELECT TOP 1 designation, ra, dec, j_m, h_m, k_m
        FROM fp_psc
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, 0.02)
        )
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None
            row = rows[0]
            result = {
                'name': object_name,
                'source': '2MASS',
                'catalog_id': row[0] if len(row) > 0 else None,
                'ra': row[1] if len(row) > 1 else None,
                'dec': row[2] if len(row) > 2 else None,
                'j_mag': row[3] if len(row) > 3 else None,
                'h_mag': row[4] if len(row) > 4 else None,
                'k_mag': row[5] if len(row) > 5 else None,
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 13: WISE ───────────────────────────────────────────────────────────────

class WISEFetcher(BaseFetcher):
    """Wide-field Infrared Survey Explorer (mid-infrared)."""
    
    def __init__(self):
        super().__init__('wise')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch mid-infrared photometry from AllWISE around the target position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        query = f"""
        SELECT TOP 1 designation, ra, dec, w1mpro, w2mpro, w3mpro, w4mpro
        FROM allwise_p3as_psd
        WHERE 1 = CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {ra}, {dec}, 0.02)
        )
        """

        try:
            rows = self._tap_sync_query(query)
            if not rows:
                return None
            row = rows[0]
            result = {
                'name': object_name,
                'source': 'WISE',
                'catalog_id': row[0] if len(row) > 0 else None,
                'ra': row[1] if len(row) > 1 else None,
                'dec': row[2] if len(row) > 2 else None,
                'w1_mag': row[3] if len(row) > 3 else None,
                'w2_mag': row[4] if len(row) > 4 else None,
                'w3_mag': row[5] if len(row) > 5 else None,
                'w4_mag': row[6] if len(row) > 6 else None,
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 14: Pan-STARRS ─────────────────────────────────────────────────────────

class PanSTARRSFetcher(BaseFetcher):
    """Panoramic Survey Telescope and Rapid Response System."""
    
    def __init__(self):
        super().__init__('panstarrs')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch Pan-STARRS DR2 mean photometry around the resolved target position."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        endpoint = f"{self.endpoint}/dr2/mean.csv"
        params = {
            'ra': ra,
            'dec': dec,
            'radius': 0.02,
            'nDetections.gte': 1,
            'pagesize': 1,
        }

        try:
            response = requests.get(endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            lines = [line for line in response.text.splitlines() if line.strip()]
            if len(lines) < 2:
                return None
            reader = csv.DictReader(lines)
            row = next(reader, None)
            if not row:
                return None

            result = {
                'name': object_name,
                'source': 'PanSTARRS',
                'catalog_id': row.get('objID') or row.get('objid'),
                'ra': self._safe_float(row.get('raMean') or row.get('ra')),
                'dec': self._safe_float(row.get('decMean') or row.get('dec')),
                'g_mag': self._null_sentinel(row.get('gMeanPSFMag') or row.get('gMeanKronMag')),
                'r_mag': self._null_sentinel(row.get('rMeanPSFMag') or row.get('rMeanKronMag')),
                'i_mag': self._null_sentinel(row.get('iMeanPSFMag') or row.get('iMeanKronMag')),
                'z_mag': self._null_sentinel(row.get('zMeanPSFMag') or row.get('zMeanKronMag')),
                'y_mag': self._null_sentinel(row.get('yMeanPSFMag') or row.get('yMeanKronMag')),
                'detection_count': self._safe_float(row.get('nDetections')),
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
            return None


# ─── API 15: ZTF ────────────────────────────────────────────────────────────────

class ZTFFetcher(BaseFetcher):
    """Zwicky Transient Facility (time-domain astronomy)."""
    
    def __init__(self):
        super().__init__('ztf')
    
    def fetch(self, object_name: str) -> Optional[Dict[str, Any]]:
        """Fetch a compact ZTF light-curve summary around the resolved target."""
        coords = self._resolve_object_coordinates(object_name)
        if not coords:
            return None

        ra, dec = coords
        params = {
            'POS': f'CIRCLE {ra} {dec} 0.001',
            'FORMAT': 'CSV',
        }

        try:
            response = requests.get(self.endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            lines = [line for line in response.text.splitlines() if line.strip()]
            if len(lines) < 2:
                return None

            reader = csv.DictReader(lines)
            rows = list(reader)
            if not rows:
                return None

            mags = [self._safe_float(row.get('mag')) for row in rows if self._safe_float(row.get('mag')) is not None]
            filters = sorted({row.get('filtercode') for row in rows if row.get('filtercode')})
            latest = rows[-1]
            result = {
                'name': object_name,
                'source': 'ZTF',
                'light_curve_points': len(rows),
                'latest_mjd': self._safe_float(latest.get('mjd')),
                'latest_mag': self._safe_float(latest.get('mag')),
                'min_mag': min(mags) if mags else None,
                'max_mag': max(mags) if mags else None,
                'mean_mag': (sum(mags) / len(mags)) if mags else None,
                'filter_bands': ','.join(filters) if filters else None,
                'ra': self._safe_float(latest.get('ra')),
                'dec': self._safe_float(latest.get('dec')),
            }
            self._log_success()
            return result
        except Exception as e:
            self._log_error(e)
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
        """Fetch ALMA observation summary for a target name."""
        params = {
            'source_name_resolver': object_name,
            'result_view': 'observation',
            'format': 'JSON',
        }

        try:
            response = requests.get(self.endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()

            payload = response.json() if 'json' in response.headers.get('Content-Type', '').lower() else response.text
            observations = []
            if isinstance(payload, dict):
                observations = payload.get('data') or payload.get('results') or payload.get('observations') or []
            elif isinstance(payload, list):
                observations = payload

            first = observations[0] if observations else {}
            self._log_success()
            return self._build_archive_result(
                object_name,
                'ALMA',
                response.url,
                'ALMA archive observation summary for millimeter/submillimeter validation.',
                payload,
                {
                    'project_code': first.get('proposal_id') or first.get('project_code'),
                    'frequency_band': first.get('band_list') or first.get('frequency_band'),
                    'target_name': first.get('target_name') or object_name,
                    'observation_date': first.get('obs_release_date') or first.get('obs_date'),
                },
            )
        except Exception as e:
            self._log_error(e)
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
        """Fetch Chandra archive search summary for a target name."""
        params = {
            'entry': object_name,
            'fields': 'all',
            'format': 'text',
        }

        try:
            response = requests.get(self.endpoint, params=params, timeout=self.timeout)
            response.raise_for_status()
            text = response.text

            observation_ids = sorted(set(re.findall(r'\b(?:obsid|ObsId)\s*[:=]?\s*([0-9]{3,6})\b', text, re.IGNORECASE)))
            energy_band_match = re.search(r'([0-9]+\.?[0-9]*)\s*[-–]\s*([0-9]+\.?[0-9]*)\s*keV', text, re.IGNORECASE)
            self._log_success()
            return self._build_archive_result(
                object_name,
                'Chandra',
                response.url,
                'Chandra archive search summary for X-ray imaging and spectroscopy.',
                text,
                {
                    'observation_ids': observation_ids[:10],
                    'energy_band_keV': energy_band_match.group(0) if energy_band_match else None,
                },
            )
        except requests.HTTPError as e:
            if getattr(e.response, 'status_code', None) == 404:
                return None
            self._log_error(e)
            return None
        except Exception as e:
            self._log_error(e)
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
        """Fetch JWST observation summary from MAST."""
        request_payload = {
            'service': 'Mast.Caom.Filtered',
            'params': {
                'obs_collection': ['JWST'],
                'target_name': object_name,
            },
            'format': 'json',
            'pagesize': 20,
            'page': 1,
        }
        headers = {'Accept': 'application/json'}
        if API_KEYS['MAST_API_KEY']:
            headers['Authorization'] = f"token {API_KEYS['MAST_API_KEY']}"

        try:
            response = requests.post(
                ENDPOINTS['mast'],
                data={'request': json.dumps(request_payload)},
                headers=headers,
                timeout=self.timeout,
            )
            response.raise_for_status()
            payload = response.json()
            observations = payload.get('data', []) if isinstance(payload, dict) else []
            first = observations[0] if observations else {}

            self._log_success()
            return self._build_archive_result(
                object_name,
                'JWST',
                response.url,
                'JWST MAST observation summary for infrared imaging and spectroscopy.',
                payload,
                {
                    'instrument_name': first.get('instrument_name'),
                    'filters': first.get('filters'),
                    'proposal_id': first.get('proposal_id'),
                    'target_name': first.get('target_name') or object_name,
                    'wavelength_region': first.get('wavelength_region'),
                },
            )
        except Exception as e:
            self._log_error(e)
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
    Unified API fetcher that tries all 55 sources in priority order.
    
    Priority Groups:
        A. Astronomical Databases (1-5): SIMBAD, NED, VizieR, Gaia, MAST
        B. NASA/Space Agencies (6-15): APOD, NeoWs, Mars, EPIC, DONKI, 
                                        Exoplanet, HEASARC, ADS, JPL, ESO
        C. Sky Surveys (16-20): SDSS, 2MASS, WISE, Pan-STARRS, ZTF
        D. Specialized (21-25): ATNF Pulsar, McGill Magnetar, TNS, GCN, GWOSC
        E. Computational/AI (26-30): Wolfram, arXiv, Grok, OpenAI, Claude
        F. Radio/Infrared (31-35): NVSS, FIRST, VLASS, ALMA, ASKAP
        G. X-ray/Gamma-ray (36-40): Chandra, XMM, Swift, Fermi, INTEGRAL
        H. Space Telescopes (41-45): HST, JWST, Spitzer, Kepler, TESS
        I. Cosmology/CMB (46-50): Planck, WMAP, DES, DESI, Euclid
        J. Spectroscopic (51-55): LAMOST, GALAH, APOGEE, RAVE, DESI Spectra
    
    Results are stored in IPData.py for QCalc.py consumption.
    """
    
    def __init__(self):
        # Primary databases (always try first)
        self.simbad = SIMBADFetcher()
        self.ned = NEDFetcher()
        self.gaia = GaiaFetcher()
        self.mast = MASTFetcher()
        self.exoplanet = ExoplanetFetcher()
        self.heasarc = HEASARCFetcher()
        self.eso = ESOFetcher()
        self.vizier = VizieRFetcher()
        self.sdss = SDSSFetcher()
        self.twomass = TwoMASSFetcher()
        self.wise = WISEFetcher()
        self.panstarrs = PanSTARRSFetcher()
        self.ztf = ZTFFetcher()
        self.grok = GrokFetcher()
        
        # NASA Services
        self.nasa_apod = NASAAPODFetcher()
        self.nasa_neows = NASANeoWsFetcher()
        self.nasa_donki = NASADONKIFetcher()
        self.alma = ALMAFetcher()
        self.chandra = ChandraFetcher()
        self.jwst = JWSTFetcher()
        self._archive_fetchers = {
            'gaia': self.gaia,
            'mast': self.mast,
            'alma': self.alma,
            'chandra': self.chandra,
            'jwst': self.jwst,
        }
        
        # All 55 fetchers (expanded from 50 with new NASA services)
        self._fetchers = {
            # Group A: Astronomical Databases (1-5)
            'simbad': self.simbad,
            'ned': self.ned,
            'vizier': self.vizier,
            'gaia': self.gaia,
            'mast': self.mast,
            # Group B: NASA/Space Agencies (6-15, expanded with new APIs)
            'nasa_apod': self.nasa_apod,           # NEW: Astronomy Picture of the Day
            'nasa_neows': self.nasa_neows,         # NEW: Near Earth Object Web Service
            'nasa_mars_weather': None,             # NEW: Mars Weather
            'nasa_epic': None,                     # NEW: Earth Polychromatic Imaging Camera
            'nasa_donki': self.nasa_donki,         # NEW: Space Weather Events
            'nasa_exoplanet': self.exoplanet,
            'heasarc': self.heasarc,
            'ads': None,
            'jpl_horizons': None,
            'eso': self.eso,
            # Group C: Sky Surveys (11-15)
            'sdss': self.sdss,
            '2mass': self.twomass,
            'wise': self.wise,
            'panstarrs': self.panstarrs,
            'ztf': self.ztf,
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
            'alma': self.alma,
            'askap': None,
            # Group G: X-ray/Gamma-ray (31-35)
            'chandra': self.chandra,
            'xmm': None,
            'swift': None,
            'fermi': None,
            'integral': None,
            # Group H: Space Telescopes (36-40)
            'hst': None,
            'jwst': self.jwst,
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

    def _merge_result_fields(self, target: Dict[str, Any], incoming: Optional[Dict[str, Any]]) -> None:
        if not incoming:
            return

        for key, value in incoming.items():
            if key in {'name', 'source', 'sources'}:
                continue
            if key not in target and value is not None and not isinstance(value, (dict, list)):
                target[key] = value

        source_name = incoming.get('source')
        if source_name and source_name not in target['sources']:
            target['sources'].append(source_name)
    
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
            self._merge_result_fields(result, simbad_data)
        
        # Try NED for extragalactic objects
        ned_data = self.ned.fetch(object_name)
        if ned_data:
            self._merge_result_fields(result, ned_data)

        archive_hits = self.fetch_archive_context(object_name)
        if archive_hits:
            result['archive_hits'] = archive_hits
            for hit in archive_hits:
                self._merge_result_fields(result, hit)

        supplemental_fetchers = [
            self.vizier,
            self.gaia,
            self.exoplanet,
            self.sdss,
            self.twomass,
            self.wise,
            self.panstarrs,
            self.ztf,
            self.mast,
        ]
        for fetcher in supplemental_fetchers:
            missing = [p for p in required_params if p not in result]
            if not missing:
                break
            self._merge_result_fields(result, fetcher.fetch(object_name))
        
        # Check for missing required parameters
        missing = [p for p in required_params if p not in result]
        
        # Use Grok fallback for missing parameters
        if missing:
            grok_data = self.grok.fetch(object_name, missing)
            if grok_data:
                self._merge_result_fields(result, grok_data)
        
        return result

    def fetch_archive_context(
        self,
        object_name: str,
        archive_names: Optional[List[str]] = None,
    ) -> List[Dict[str, Any]]:
        """Fetch archive summaries from the implemented archive bridges."""
        requested = {
            name.strip().lower()
            for name in (archive_names or list(self._archive_fetchers.keys()))
            if name and name.strip()
        }

        archive_hits: List[Dict[str, Any]] = []
        for name, fetcher in self._archive_fetchers.items():
            if name not in requested or fetcher is None:
                continue

            archive_result = fetcher.fetch(object_name)
            if archive_result:
                archive_hits.append(archive_result)

        return archive_hits
    
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
    
    def fetch_and_save(
        self,
        object_name: str,
        output_dir: str = ".",
        required_params: List[str] = None,
    ) -> Tuple[Dict[str, Any], str]:
        """
        Fetch parameters and save to timestamped CSV file.
        Also stores in IPData.py.
        
        Args:
            object_name: Name of the object to query
            output_dir: Directory for output CSV
            
        Returns:
            Tuple of (result_dict, csv_filepath)
        """
        result = self.fetch(object_name, required_params)
        
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


def fetch_and_save(
    object_name: str,
    output_dir: str = ".",
    required_params: List[str] = None,
) -> Tuple[Dict[str, Any], str]:
    """Fetch parameters and save to CSV."""
    return FETCHER.fetch_and_save(object_name, output_dir, required_params)


def _parse_csv_list(raw_value: str) -> List[str]:
    return [item.strip() for item in (raw_value or '').split(',') if item.strip()]


def run_cli(argv: Optional[List[str]] = None) -> int:
    """Command-line entry point for GUI/archive bridge execution."""
    parser = argparse.ArgumentParser(description="Fetch astronomical data and archive summaries for Star-Magic.")
    parser.add_argument('--object', dest='object_name', help='Astronomical object name to query.')
    parser.add_argument('--output-dir', default='.', help='Directory for generated CSV output when --save is used.')
    parser.add_argument('--required-params', default='', help='Comma-separated parameter names to prioritise.')
    parser.add_argument('--archives', default='gaia,mast,alma,chandra,jwst', help='Comma-separated archive bridges to query.')
    parser.add_argument('--save', action='store_true', help='Write a bodies_*.csv file in addition to JSON output.')
    parser.add_argument('--pretty', action='store_true', help='Pretty-print JSON output.')
    parser.add_argument('--manual', action='store_true', help='Use interactive manual input instead of API fetch.')
    args = parser.parse_args(argv)

    required_params = _parse_csv_list(args.required_params)
    archive_names = _parse_csv_list(args.archives)

    try:
        if args.manual:
            result = manual_input()
            csv_path = None
        else:
            object_name = (args.object_name or '').strip()
            if not object_name:
                parser.error('--object is required unless --manual is used')

            if args.save:
                result, csv_path = FETCHER.fetch_and_save(object_name, args.output_dir, required_params)
            else:
                result = FETCHER.fetch(object_name, required_params)
                csv_path = None

            if archive_names:
                result['archive_hits'] = FETCHER.fetch_archive_context(object_name, archive_names)
                source_names = {source for source in result.get('sources', []) if source}
                for hit in result['archive_hits']:
                    source = hit.get('source')
                    if source:
                        source_names.add(source)
                result['sources'] = sorted(source_names)

        payload = {
            'success': True,
            'result': result,
            'csv_path': csv_path,
        }
    except Exception as e:
        payload = {
            'success': False,
            'error': str(e),
        }

    print(json.dumps(payload, indent=2 if args.pretty else None))
    return 0 if payload.get('success') else 1


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    sys.exit(run_cli())
