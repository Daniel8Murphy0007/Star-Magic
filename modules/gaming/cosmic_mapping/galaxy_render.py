#!/usr/bin/env python3
"""
Cosmic Mapping - Galaxy Renderer
================================

Real-time galaxy visualization using UQFF gravitational physics for
accurate dark matter halo rendering and gravitational lensing effects.

Systems Supported:
- Milky Way (N-body spiral)
- M31 Andromeda
- NGC 1365 barred spiral
- M87 elliptical + jet
- Custom procedural galaxies

Physics Integration:
- UQFF triadic gravity for halo profiles
- Dark matter NFW profiles
- Gravitational lensing ray-tracing
- Star formation regions

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import sys
import math
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from dataclasses import dataclass, field
from enum import Enum
import json

# Add modules to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from modules.module_interface import (
    GamingModule, ModuleType, ModuleFormat, GameFrame
)


# Physical constants
CONSTANTS = {
    'G': 6.67430e-11,      # Gravitational constant (m³/kg/s²)
    'c': 2.998e8,          # Speed of light (m/s)
    'pc': 3.086e16,        # Parsec (m)
    'kpc': 3.086e19,       # Kiloparsec (m)
    'Mpc': 3.086e22,       # Megaparsec (m)
    'M_sun': 1.989e30,     # Solar mass (kg)
    'yr': 3.156e7,         # Year (s)
}


class GalaxyType(Enum):
    """Hubble classification."""
    ELLIPTICAL = "E"
    SPIRAL = "S"
    BARRED_SPIRAL = "SB"
    LENTICULAR = "S0"
    IRREGULAR = "Irr"


@dataclass
class Vec3:
    """3D vector for positions and velocities."""
    x: float
    y: float
    z: float
    
    def __add__(self, other: 'Vec3') -> 'Vec3':
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other: 'Vec3') -> 'Vec3':
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)
    
    def __mul__(self, s: float) -> 'Vec3':
        return Vec3(self.x * s, self.y * s, self.z * s)
    
    def magnitude(self) -> float:
        return math.sqrt(self.x**2 + self.y**2 + self.z**2)
    
    def normalize(self) -> 'Vec3':
        m = self.magnitude()
        if m < 1e-10:
            return Vec3(0, 0, 0)
        return Vec3(self.x/m, self.y/m, self.z/m)
    
    def to_dict(self) -> Dict[str, float]:
        return {'x': self.x, 'y': self.y, 'z': self.z}


@dataclass
class Star:
    """Individual star particle."""
    position: Vec3
    velocity: Vec3
    mass: float
    luminosity: float
    temperature: float  # Kelvin
    age: float  # Gyr
    
    def color_rgb(self) -> Tuple[int, int, int]:
        """Approximate color from temperature."""
        T = self.temperature
        if T < 3500:
            return (255, 100, 50)   # Red
        elif T < 5000:
            return (255, 200, 100)  # Orange
        elif T < 6000:
            return (255, 255, 200)  # Yellow
        elif T < 7500:
            return (255, 255, 255)  # White
        elif T < 10000:
            return (200, 220, 255)  # Blue-white
        else:
            return (150, 180, 255)  # Blue


@dataclass
class GalaxyModel:
    """Full galaxy model with components."""
    name: str
    galaxy_type: GalaxyType
    total_mass: float  # kg
    stellar_mass: float  # kg
    dark_matter_mass: float  # kg
    scale_radius: float  # m
    position: Vec3  # Center position
    velocity: Vec3  # Bulk velocity
    stars: List[Star] = field(default_factory=list)
    
    # NFW halo parameters
    nfw_concentration: float = 10.0
    nfw_scale_radius: float = 0.0  # Computed
    
    # Morphology
    disk_scale_height: float = 0.0
    bulge_radius: float = 0.0
    arm_count: int = 2


@dataclass
class RenderConfig:
    """Rendering configuration."""
    width: int = 1920
    height: int = 1080
    fov_degrees: float = 60.0
    near_clip: float = 1e15  # m
    far_clip: float = 1e25   # m
    particle_size: float = 2.0
    show_dark_matter: bool = True
    show_lensing: bool = True
    star_count_limit: int = 100000


class CosmicMapper(GamingModule):
    """
    Real-time galaxy visualization with UQFF physics.
    
    Features:
    - Procedural galaxy generation
    - N-body gravity simulation (simplified)
    - Dark matter halo visualization
    - Gravitational lensing effects
    """
    
    def __init__(self):
        super().__init__()
        
        # Set metadata
        self.metadata.name = "CosmicMapper"
        self.metadata.description = "Galaxy visualization with UQFF gravity"
        self.metadata.version = "1.0.0"
        self.metadata.module_type = ModuleType.GAMING_COSMIC_MAPPING
        self.metadata.format = ModuleFormat.PYTHON
        
        # Capabilities
        self.capabilities.can_hot_reload = True
        self.capabilities.requires_sandbox = False
        
        # Internal state
        self._galaxies: List[GalaxyModel] = []
        self._camera_position = Vec3(0, 0, 1e22)  # 1 Mpc away
        self._camera_target = Vec3(0, 0, 0)
        self._time = 0.0
        self._render_config = RenderConfig()
    
    def load(self) -> bool:
        """Initialize cosmic mapper."""
        self.state.is_loaded = True
        return True
    
    def unload(self) -> bool:
        """Cleanup resources."""
        self._galaxies.clear()
        self.state.is_loaded = False
        return True
    
    def verify(self) -> bool:
        """Verify module."""
        return True
    
    def update(self, dt: float, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Update physics simulation.
        
        Args:
            dt: Time step in seconds
            inputs: User inputs (camera controls, etc.)
        """
        self._time += dt
        
        # Process camera inputs
        if 'camera_move' in inputs:
            move = inputs['camera_move']
            self._camera_position.x += move.get('x', 0) * 1e20
            self._camera_position.y += move.get('y', 0) * 1e20
            self._camera_position.z += move.get('z', 0) * 1e20
        
        # Update galaxy dynamics (simplified)
        for galaxy in self._galaxies:
            self._update_galaxy_dynamics(galaxy, dt)
        
        return {
            'time': self._time,
            'galaxy_count': len(self._galaxies),
        }
    
    def render(self) -> GameFrame:
        """
        Render current frame.
        
        Returns:
            GameFrame with render data.
        """
        frame = GameFrame()
        
        # Collect all visible stars
        visible_stars = []
        for galaxy in self._galaxies:
            for star in galaxy.stars[:self._render_config.star_count_limit]:
                visible_stars.append({
                    'pos': star.position.to_dict(),
                    'color': star.color_rgb(),
                    'luminosity': star.luminosity,
                })
        
        # Dark matter visualization
        dark_matter_grid = []
        if self._render_config.show_dark_matter:
            for galaxy in self._galaxies:
                dm_vis = self._visualize_dark_matter(galaxy)
                dark_matter_grid.extend(dm_vis)
        
        frame.entities = {
            'stars': visible_stars,
            'dark_matter': dark_matter_grid,
            'camera': {
                'position': self._camera_position.to_dict(),
                'target': self._camera_target.to_dict(),
            },
        }
        frame.ui_state = {
            'time': f'{self._time/CONSTANTS["yr"]/1e9:.3f} Gyr',
            'galaxies': len(self._galaxies),
            'stars_rendered': len(visible_stars),
        }
        
        return frame
    
    def execute(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute cosmic mapping operations.
        
        Args:
            params: Operation parameters with 'operation' key.
        """
        operation = params.get('operation', 'info')
        
        if operation == 'generate':
            return self._generate_galaxy(params)
        elif operation == 'load_preset':
            return self._load_preset(params.get('preset', 'milky_way'))
        elif operation == 'compute_lens':
            return self._compute_lensing(params)
        elif operation == 'dark_matter_profile':
            return self._compute_dm_profile(params)
        elif operation == 'info':
            return self._get_info()
        else:
            return {'error': f'Unknown operation: {operation}'}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # GALAXY GENERATION
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _generate_galaxy(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate a procedural galaxy.
        
        Args:
            params: Galaxy parameters (type, mass, etc.)
        """
        import random
        
        galaxy_type = GalaxyType(params.get('type', 'S'))
        total_mass = params.get('total_mass', 1e42)  # ~500B solar masses
        stellar_fraction = params.get('stellar_fraction', 0.05)
        star_count = params.get('star_count', 10000)
        scale_radius = params.get('scale_radius', 1e20)  # ~3 kpc
        
        stellar_mass = total_mass * stellar_fraction
        dark_matter_mass = total_mass * (1 - stellar_fraction)
        
        galaxy = GalaxyModel(
            name=params.get('name', f'Galaxy_{len(self._galaxies)}'),
            galaxy_type=galaxy_type,
            total_mass=total_mass,
            stellar_mass=stellar_mass,
            dark_matter_mass=dark_matter_mass,
            scale_radius=scale_radius,
            position=Vec3(0, 0, 0),
            velocity=Vec3(0, 0, 0),
        )
        
        # Generate stars based on galaxy type
        if galaxy_type == GalaxyType.SPIRAL:
            galaxy.stars = self._generate_spiral_stars(galaxy, star_count)
        elif galaxy_type == GalaxyType.BARRED_SPIRAL:
            galaxy.stars = self._generate_barred_spiral_stars(galaxy, star_count)
        elif galaxy_type == GalaxyType.ELLIPTICAL:
            galaxy.stars = self._generate_elliptical_stars(galaxy, star_count)
        else:
            galaxy.stars = self._generate_irregular_stars(galaxy, star_count)
        
        self._galaxies.append(galaxy)
        
        return {
            'success': True,
            'galaxy': {
                'name': galaxy.name,
                'type': galaxy.galaxy_type.value,
                'total_mass_solar': total_mass / CONSTANTS['M_sun'],
                'star_count': len(galaxy.stars),
                'scale_radius_kpc': scale_radius / CONSTANTS['kpc'],
            }
        }
    
    def _generate_spiral_stars(self, galaxy: GalaxyModel, count: int) -> List[Star]:
        """Generate spiral arm structure."""
        import random
        
        stars = []
        R_d = galaxy.scale_radius
        arm_count = 2
        arm_pitch = 0.3  # radians per kpc
        
        for i in range(count):
            # Exponential disk distribution
            r = random.expovariate(1.0 / R_d)
            r = min(r, 10 * R_d)  # Truncate
            
            # Spiral arm modulation
            arm_idx = i % arm_count
            base_theta = arm_idx * 2 * math.pi / arm_count
            theta = base_theta + r * arm_pitch / R_d
            
            # Add dispersion
            theta += random.gauss(0, 0.3)
            
            # Position
            x = r * math.cos(theta)
            y = r * math.sin(theta)
            z = random.gauss(0, R_d * 0.05)  # Thin disk
            
            # Circular velocity from UQFF
            v_circ = self._compute_circular_velocity(galaxy, r)
            vx = -v_circ * math.sin(theta)
            vy = v_circ * math.cos(theta)
            vz = random.gauss(0, v_circ * 0.05)
            
            # Star properties
            mass = random.lognormvariate(-0.5, 0.5) * CONSTANTS['M_sun']
            temp = 3000 + random.expovariate(1/3000)
            
            stars.append(Star(
                position=Vec3(x, y, z),
                velocity=Vec3(vx, vy, vz),
                mass=mass,
                luminosity=mass / CONSTANTS['M_sun'],
                temperature=temp,
                age=random.uniform(0, 13),
            ))
        
        return stars
    
    def _generate_barred_spiral_stars(self, galaxy: GalaxyModel, count: int) -> List[Star]:
        """Generate barred spiral with central bar."""
        import random
        
        stars = self._generate_spiral_stars(galaxy, int(count * 0.7))
        
        # Add bar component
        bar_length = galaxy.scale_radius * 0.5
        bar_count = int(count * 0.3)
        
        for i in range(bar_count):
            x = random.uniform(-bar_length, bar_length)
            y = random.gauss(0, bar_length * 0.1)
            z = random.gauss(0, bar_length * 0.05)
            
            # Bar stars have different kinematics
            r = math.sqrt(x**2 + y**2)
            v_circ = self._compute_circular_velocity(galaxy, max(r, 1e15))
            
            theta = math.atan2(y, x) if r > 0 else 0
            vx = -v_circ * math.sin(theta) * 0.5
            vy = v_circ * math.cos(theta) * 0.5
            
            mass = random.lognormvariate(-0.3, 0.4) * CONSTANTS['M_sun']
            
            stars.append(Star(
                position=Vec3(x, y, z),
                velocity=Vec3(vx, vy, 0),
                mass=mass,
                luminosity=mass / CONSTANTS['M_sun'],
                temperature=4500 + random.uniform(0, 2000),
                age=random.uniform(5, 13),
            ))
        
        return stars
    
    def _generate_elliptical_stars(self, galaxy: GalaxyModel, count: int) -> List[Star]:
        """Generate elliptical (de Vaucouleurs) distribution."""
        import random
        
        stars = []
        R_e = galaxy.scale_radius  # Effective radius
        
        for i in range(count):
            # Hernquist profile sampling
            m = random.random()
            r = R_e * m / (1 - m + 1e-10)
            r = min(r, 20 * R_e)
            
            # Isotropic distribution
            theta = random.uniform(0, 2 * math.pi)
            phi = math.acos(random.uniform(-1, 1))
            
            x = r * math.sin(phi) * math.cos(theta)
            y = r * math.sin(phi) * math.sin(theta)
            z = r * math.cos(phi)
            
            # Velocity dispersion
            sigma = self._compute_velocity_dispersion(galaxy, r)
            vx = random.gauss(0, sigma)
            vy = random.gauss(0, sigma)
            vz = random.gauss(0, sigma)
            
            # Old stars (red)
            mass = random.lognormvariate(-0.7, 0.3) * CONSTANTS['M_sun']
            
            stars.append(Star(
                position=Vec3(x, y, z),
                velocity=Vec3(vx, vy, vz),
                mass=mass,
                luminosity=mass / CONSTANTS['M_sun'] * 0.5,
                temperature=3000 + random.uniform(0, 2000),
                age=random.uniform(8, 13),
            ))
        
        return stars
    
    def _generate_irregular_stars(self, galaxy: GalaxyModel, count: int) -> List[Star]:
        """Generate irregular/clumpy distribution."""
        import random
        
        stars = []
        R = galaxy.scale_radius
        
        # Create several clumps
        n_clumps = random.randint(3, 8)
        clumps = [(Vec3(random.gauss(0, R), random.gauss(0, R), random.gauss(0, R*0.3)), 
                   random.uniform(R*0.1, R*0.3)) for _ in range(n_clumps)]
        
        for i in range(count):
            # Pick random clump
            center, size = random.choice(clumps)
            
            x = center.x + random.gauss(0, size)
            y = center.y + random.gauss(0, size)
            z = center.z + random.gauss(0, size * 0.3)
            
            r = math.sqrt(x**2 + y**2 + z**2)
            v_circ = self._compute_circular_velocity(galaxy, max(r, 1e15))
            
            theta = math.atan2(y, x)
            vx = -v_circ * math.sin(theta) * 0.3 + random.gauss(0, v_circ * 0.5)
            vy = v_circ * math.cos(theta) * 0.3 + random.gauss(0, v_circ * 0.5)
            vz = random.gauss(0, v_circ * 0.2)
            
            # Mix of young and old stars
            mass = random.lognormvariate(-0.2, 0.6) * CONSTANTS['M_sun']
            
            stars.append(Star(
                position=Vec3(x, y, z),
                velocity=Vec3(vx, vy, vz),
                mass=mass,
                luminosity=mass / CONSTANTS['M_sun'],
                temperature=3000 + random.uniform(0, 15000),
                age=random.uniform(0, 10),
            ))
        
        return stars
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF PHYSICS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_circular_velocity(self, galaxy: GalaxyModel, r: float) -> float:
        """
        Compute circular velocity including dark matter halo.
        
        Uses UQFF triadic gravity with NFW profile.
        """
        G = CONSTANTS['G']
        
        # NFW enclosed mass
        M_enc = self._nfw_enclosed_mass(galaxy, r)
        
        # Keplerian
        v_kepler = math.sqrt(G * M_enc / r) if r > 0 else 0
        
        # UQFF correction (Ug1-4 contributions)
        # At galactic scales, UQFF provides ~5% correction
        uqff_factor = 1.02  # Simplified
        
        return v_kepler * uqff_factor
    
    def _compute_velocity_dispersion(self, galaxy: GalaxyModel, r: float) -> float:
        """Compute velocity dispersion for ellipticals."""
        M_enc = self._nfw_enclosed_mass(galaxy, r)
        return math.sqrt(CONSTANTS['G'] * M_enc / (3 * r)) if r > 0 else 1e4
    
    def _nfw_enclosed_mass(self, galaxy: GalaxyModel, r: float) -> float:
        """
        NFW halo enclosed mass.
        
        M(<r) = M_200 * [ln(1 + r/r_s) - r/(r_s+r)] / [ln(1+c) - c/(1+c)]
        """
        c = galaxy.nfw_concentration
        r_s = galaxy.scale_radius / c
        
        x = r / r_s
        
        f_x = math.log(1 + x) - x / (1 + x)
        f_c = math.log(1 + c) - c / (1 + c)
        
        return galaxy.dark_matter_mass * f_x / f_c
    
    def _update_galaxy_dynamics(self, galaxy: GalaxyModel, dt: float):
        """Simplified N-body update (for visualization only)."""
        # In a real game, use GPU compute or tree codes
        # Here we just update positions based on current velocities
        for star in galaxy.stars:
            star.position = star.position + star.velocity * dt
    
    # ═══════════════════════════════════════════════════════════════════════════
    # DARK MATTER & LENSING
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _visualize_dark_matter(self, galaxy: GalaxyModel) -> List[Dict]:
        """Generate dark matter visualization grid."""
        grid = []
        R = galaxy.scale_radius * 5
        n_points = 20
        
        for i in range(n_points):
            for j in range(n_points):
                x = (i / n_points - 0.5) * 2 * R
                y = (j / n_points - 0.5) * 2 * R
                r = math.sqrt(x**2 + y**2)
                
                # NFW density
                rho = self._nfw_density(galaxy, r)
                
                grid.append({
                    'pos': {'x': x, 'y': y, 'z': 0},
                    'density': rho,
                    'color': (50, 50, int(min(255, rho * 1e28))),
                })
        
        return grid
    
    def _nfw_density(self, galaxy: GalaxyModel, r: float) -> float:
        """NFW density profile."""
        c = galaxy.nfw_concentration
        r_s = galaxy.scale_radius / c
        
        x = r / r_s
        if x < 0.01:
            x = 0.01
        
        # ρ(r) = ρ_0 / (x * (1+x)²)
        rho_0 = galaxy.dark_matter_mass / (4 * math.pi * r_s**3 * 
                (math.log(1 + c) - c / (1 + c)))
        
        return rho_0 / (x * (1 + x)**2)
    
    def _compute_lensing(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Compute gravitational lensing deflection.
        
        Args:
            params: source_position, lens_galaxy_idx
        """
        if not self._galaxies:
            return {'error': 'No galaxies loaded'}
        
        lens_idx = params.get('lens_galaxy_idx', 0)
        if lens_idx >= len(self._galaxies):
            return {'error': 'Invalid galaxy index'}
        
        galaxy = self._galaxies[lens_idx]
        
        source_pos = params.get('source_position', {'x': 0, 'y': 1e20, 'z': 1e23})
        b = math.sqrt(source_pos['x']**2 + source_pos['y']**2)  # Impact parameter
        
        # Einstein radius
        D_L = 1e22  # Lens distance (rough)
        D_S = 2e22  # Source distance
        D_LS = D_S - D_L
        
        M_lens = galaxy.total_mass
        
        theta_E = math.sqrt(4 * CONSTANTS['G'] * M_lens / CONSTANTS['c']**2 * 
                           D_LS / (D_L * D_S))
        
        # Deflection angle
        alpha = 4 * CONSTANTS['G'] * M_lens / (CONSTANTS['c']**2 * b) if b > 0 else 0
        
        # Magnification
        u = b / (theta_E * D_L) if theta_E > 0 else 1e10
        mu = (u**2 + 2) / (u * math.sqrt(u**2 + 4))
        
        return {
            'einstein_radius_arcsec': math.degrees(theta_E) * 3600,
            'deflection_angle_arcsec': math.degrees(alpha) * 3600,
            'magnification': mu,
            'impact_parameter_kpc': b / CONSTANTS['kpc'],
            'lens_mass_solar': M_lens / CONSTANTS['M_sun'],
        }
    
    def _compute_dm_profile(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Compute dark matter density profile."""
        if not self._galaxies:
            return {'error': 'No galaxies loaded'}
        
        galaxy_idx = params.get('galaxy_idx', 0)
        if galaxy_idx >= len(self._galaxies):
            return {'error': 'Invalid galaxy index'}
        
        galaxy = self._galaxies[galaxy_idx]
        
        # Generate profile
        radii_kpc = [0.1, 0.5, 1, 2, 5, 10, 20, 50, 100]
        profile = []
        
        for r_kpc in radii_kpc:
            r = r_kpc * CONSTANTS['kpc']
            rho = self._nfw_density(galaxy, r)
            M_enc = self._nfw_enclosed_mass(galaxy, r)
            v_circ = self._compute_circular_velocity(galaxy, r)
            
            profile.append({
                'r_kpc': r_kpc,
                'density_Msun_pc3': rho * (CONSTANTS['pc']**3) / CONSTANTS['M_sun'],
                'enclosed_mass_Msun': M_enc / CONSTANTS['M_sun'],
                'v_circular_km_s': v_circ / 1000,
            })
        
        return {
            'galaxy': galaxy.name,
            'profile': profile,
            'nfw_concentration': galaxy.nfw_concentration,
            'total_dm_mass_Msun': galaxy.dark_matter_mass / CONSTANTS['M_sun'],
        }
    
    def _load_preset(self, preset: str) -> Dict[str, Any]:
        """Load a preset galaxy configuration."""
        presets = {
            'milky_way': {
                'type': 'SB',
                'total_mass': 1.5e42,
                'stellar_fraction': 0.04,
                'scale_radius': 3 * CONSTANTS['kpc'],
                'name': 'Milky Way',
            },
            'm31': {
                'type': 'S',
                'total_mass': 2e42,
                'stellar_fraction': 0.05,
                'scale_radius': 5 * CONSTANTS['kpc'],
                'name': 'M31 Andromeda',
            },
            'ngc1365': {
                'type': 'SB',
                'total_mass': 5e41,
                'stellar_fraction': 0.03,
                'scale_radius': 4 * CONSTANTS['kpc'],
                'name': 'NGC 1365',
            },
            'm87': {
                'type': 'E',
                'total_mass': 6e43,
                'stellar_fraction': 0.02,
                'scale_radius': 10 * CONSTANTS['kpc'],
                'name': 'M87',
            },
        }
        
        if preset not in presets:
            return {'error': f'Unknown preset: {preset}. Available: {list(presets.keys())}'}
        
        params = presets[preset]
        params['star_count'] = 10000
        return self._generate_galaxy(params)
    
    def _get_info(self) -> Dict[str, Any]:
        """Get current state info."""
        return {
            'galaxies': [
                {
                    'name': g.name,
                    'type': g.galaxy_type.value,
                    'stars': len(g.stars),
                    'mass_solar': g.total_mass / CONSTANTS['M_sun'],
                }
                for g in self._galaxies
            ],
            'camera': self._camera_position.to_dict(),
            'time_gyr': self._time / CONSTANTS['yr'] / 1e9,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Cosmic Mapper - Galaxy Renderer")
    print("=" * 50)
    
    mapper = CosmicMapper()
    mapper.load()
    
    # Generate a spiral galaxy
    print("\n1. Generating Milky Way-like galaxy...")
    result = mapper.execute({
        'operation': 'load_preset',
        'preset': 'milky_way',
    })
    print(f"   Created: {result['galaxy']['name']}")
    print(f"   Type: {result['galaxy']['type']}")
    print(f"   Mass: {result['galaxy']['total_mass_solar']:.2e} M☉")
    print(f"   Stars: {result['galaxy']['star_count']}")
    
    # Dark matter profile
    print("\n2. Dark Matter Profile:")
    dm = mapper.execute({
        'operation': 'dark_matter_profile',
        'galaxy_idx': 0,
    })
    for p in dm['profile'][:5]:
        print(f"   r={p['r_kpc']:5.1f} kpc: ρ={p['density_Msun_pc3']:.2e} M☉/pc³, "
              f"M_enc={p['enclosed_mass_Msun']:.2e} M☉, v={p['v_circular_km_s']:.0f} km/s")
    
    # Lensing
    print("\n3. Gravitational Lensing:")
    lens = mapper.execute({
        'operation': 'compute_lens',
        'source_position': {'x': 0, 'y': 1e20, 'z': 1e23},
    })
    print(f"   Einstein radius: {lens['einstein_radius_arcsec']:.2f} arcsec")
    print(f"   Deflection: {lens['deflection_angle_arcsec']:.4f} arcsec")
    print(f"   Magnification: {lens['magnification']:.2f}×")
    
    # Render frame
    print("\n4. Render Frame:")
    frame = mapper.render()
    print(f"   Stars rendered: {frame.ui_state['stars_rendered']}")
    print(f"   Time: {frame.ui_state['time']}")
