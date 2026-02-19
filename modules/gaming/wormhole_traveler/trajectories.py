#!/usr/bin/env python3
"""
Wormhole Traveler - Wormhole Trajectory Calculator
===================================================

Computes traversable wormhole trajectories using Morris-Thorne metrics
and UQFF exotic matter requirements.

Physics Background:
    Morris-Thorne wormholes require exotic matter with negative energy density
    to keep the throat open. The UQFF framework provides the exotic matter
    computation through Ug4 (vacuum concentration) contributions.

Key Equations:
    - Throat metric: ds² = -e^(2Φ)dt² + dr²/(1-b(r)/r) + r²(dθ² + sin²θ dφ²)
    - Exotic matter: ρ < 0 (negative energy density)
    - Throat stability: d²b/dr² > 0 at throat

Usage:
    from trajectories import WormholeTrajectoryCalculator
    
    calc = WormholeTrajectoryCalculator()
    calc.load()
    
    trajectory = calc.execute({
        'throat_radius': 1e4,  # meters
        'entry_velocity': 0.01,  # fraction of c
        'num_steps': 1000,
    })

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Gaming Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import sys
import math
import json
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
from dataclasses import dataclass, field
import numpy as np

# Add modules to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from modules.module_interface import (
    GamingModule, ModuleType, ModuleFormat, GameFrame
)

# Physical constants
CONSTANTS = {
    'G': 6.67430e-11,      # Gravitational constant
    'c': 2.998e8,          # Speed of light
    'hbar': 1.055e-34,     # Reduced Planck constant
    'M_sun': 1.989e30,     # Solar mass
}


@dataclass
class WormholeParams:
    """Wormhole configuration parameters."""
    throat_radius: float = 1e4        # Throat radius b₀ (meters)
    shape_function: str = 'ellis'     # 'ellis' or 'morris_thorne'
    redshift_function: float = 0.0    # Φ(r) at throat
    exotic_matter_density: float = 0.0  # Computed
    
    # Stability parameters
    is_traversable: bool = True
    minimum_velocity: float = 0.0
    maximum_velocity: float = 0.1  # Fraction of c


@dataclass
class TrajectoryPoint:
    """Single point on a wormhole trajectory."""
    proper_time: float     # τ (seconds)
    coordinate_time: float  # t (seconds)
    radial_position: float  # r (meters)
    theta: float           # θ (radians)
    phi: float             # φ (radians)
    velocity: float        # dr/dτ (m/s)
    energy: float          # Specific energy E
    angular_momentum: float # Specific angular momentum L


@dataclass
class WormholeTrajectory:
    """Complete wormhole trajectory."""
    params: WormholeParams
    entry_point: TrajectoryPoint
    exit_point: TrajectoryPoint
    points: List[TrajectoryPoint] = field(default_factory=list)
    total_proper_time: float = 0.0
    total_coordinate_time: float = 0.0
    transit_successful: bool = False
    exotic_matter_encountered: float = 0.0


class WormholeTrajectoryCalculator(GamingModule):
    """
    Computes traversable wormhole trajectories for the Wormhole Traveler game.
    
    Uses Morris-Thorne wormhole geometry with UQFF exotic matter calculations.
    """
    
    def __init__(self):
        super().__init__()
        
        # Set metadata
        self.metadata.name = "WormholeTrajectoryCalculator"
        self.metadata.description = "Morris-Thorne wormhole trajectory computation"
        self.metadata.version = "1.0.0"
        self.metadata.module_type = ModuleType.GAMING_WORMHOLE
        self.metadata.format = ModuleFormat.PYTHON
        
        # Capabilities
        self.capabilities.requires_gpu = False  # Pure CPU calculation
        self.capabilities.max_memory_mb = 512
        self.capabilities.supports_jit = True  # Can use Numba
        
        # Game state
        self.game_state = {
            'player_position': [0.0, 0.0, 0.0],
            'player_velocity': [0.0, 0.0, 0.0],
            'wormholes_discovered': [],
            'total_distance_traveled': 0.0,
        }
        
        # Current wormhole
        self._current_wormhole: Optional[WormholeParams] = None
    
    def load(self) -> bool:
        """Initialize the trajectory calculator."""
        self.state.is_loaded = True
        return True
    
    def unload(self) -> bool:
        """Cleanup resources."""
        self.state.is_loaded = False
        return True
    
    def verify(self) -> bool:
        """Verify module integrity."""
        # Test with default Ellis wormhole
        try:
            b = self._shape_function_ellis(1e4, 1e4)
            return b == 1e4  # At throat, b(r₀) = r₀
        except Exception:
            return False
    
    def update(self, delta_time: float) -> None:
        """Update game state for one frame."""
        # Move player along current trajectory
        if hasattr(self, '_current_trajectory') and self._current_trajectory:
            # Advance position
            self.game_state['total_distance_traveled'] += delta_time * 1000  # km
    
    def render(self) -> bytes:
        """Render current frame (returns placeholder data)."""
        # In real implementation, this would interface with VTK/OpenGL
        frame_data = json.dumps({
            'game_state': self.game_state,
            'current_wormhole': self._current_wormhole.__dict__ if self._current_wormhole else None,
        }).encode('utf-8')
        return frame_data
    
    def calculate(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Calculate a wormhole trajectory.
        
        Args:
            params: Dict with:
                - throat_radius: Wormhole throat radius (m)
                - entry_velocity: Initial velocity (fraction of c)
                - entry_angle: Approach angle (rad)
                - num_steps: Number of integration steps
                
        Returns:
            Trajectory data and exotic matter requirements.
        """
        # Extract parameters
        throat_radius = params.get('throat_radius', 1e4)
        entry_velocity = params.get('entry_velocity', 0.01)
        entry_angle = params.get('entry_angle', 0.0)
        num_steps = params.get('num_steps', 1000)
        shape_type = params.get('shape_function', 'ellis')
        
        # Create wormhole
        wormhole = WormholeParams(
            throat_radius=throat_radius,
            shape_function=shape_type,
        )
        
        # Compute exotic matter requirement
        wormhole.exotic_matter_density = self._compute_exotic_matter(wormhole)
        
        # Check traversability
        wormhole.minimum_velocity = self._minimum_velocity(wormhole)
        wormhole.is_traversable = entry_velocity >= wormhole.minimum_velocity
        
        self._current_wormhole = wormhole
        
        if not wormhole.is_traversable:
            return {
                'success': False,
                'error': 'Wormhole not traversable at given velocity',
                'minimum_velocity': wormhole.minimum_velocity,
                'provided_velocity': entry_velocity,
            }
        
        # Compute trajectory
        trajectory = self._compute_trajectory(
            wormhole=wormhole,
            entry_velocity=entry_velocity,
            entry_angle=entry_angle,
            num_steps=num_steps,
        )
        
        return {
            'success': True,
            'trajectory': {
                'num_points': len(trajectory.points),
                'total_proper_time': trajectory.total_proper_time,
                'total_coordinate_time': trajectory.total_coordinate_time,
                'transit_successful': trajectory.transit_successful,
                'exotic_matter_encountered': trajectory.exotic_matter_encountered,
            },
            'wormhole': {
                'throat_radius': wormhole.throat_radius,
                'exotic_matter_density': wormhole.exotic_matter_density,
                'is_traversable': wormhole.is_traversable,
            },
            'entry_point': trajectory.entry_point.__dict__,
            'exit_point': trajectory.exit_point.__dict__,
            'sample_points': [p.__dict__ for p in trajectory.points[::max(1, len(trajectory.points)//10)]],
        }
    
    # ═══════════════════════════════════════════════════════════════════════════
    # WORMHOLE GEOMETRY
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _shape_function_ellis(self, r: float, b0: float) -> float:
        """
        Ellis wormhole shape function: b(r) = b₀²/r
        
        This gives a smooth, traversable geometry.
        """
        return b0**2 / r
    
    def _shape_function_morris_thorne(self, r: float, b0: float) -> float:
        """
        Morris-Thorne shape function: b(r) = b₀
        
        Constant shape function (simplest case).
        """
        return b0
    
    def _redshift_function(self, r: float, Phi0: float = 0.0) -> float:
        """
        Redshift function Φ(r).
        
        For zero tidal force, Φ = const.
        """
        return Phi0
    
    def _metric_coefficient(self, r: float, wormhole: WormholeParams) -> float:
        """
        Compute metric coefficient: g_rr = 1 / (1 - b(r)/r)
        """
        b0 = wormhole.throat_radius
        
        if wormhole.shape_function == 'ellis':
            b_r = self._shape_function_ellis(r, b0)
        else:
            b_r = self._shape_function_morris_thorne(r, b0)
        
        # Avoid singularity at throat
        if abs(r - b_r) < 1e-10:
            return 1e10  # Large but finite
        
        return 1.0 / (1.0 - b_r / r)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # EXOTIC MATTER
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_exotic_matter(self, wormhole: WormholeParams) -> float:
        """
        Compute exotic matter (negative energy) density at throat.
        
        For Morris-Thorne: ρ = -c²/(8πG) * b'(r₀)/r₀²
        
        Returns negative value indicating exotic matter requirement.
        """
        G = CONSTANTS['G']
        c = CONSTANTS['c']
        b0 = wormhole.throat_radius
        
        # Shape function derivative at throat
        if wormhole.shape_function == 'ellis':
            # b(r) = b₀²/r → b'(r) = -b₀²/r²
            b_prime = -b0**2 / b0**2  # At r = b₀: b' = -1
        else:
            # b(r) = b₀ → b'(r) = 0
            b_prime = 0.0
        
        # Exotic matter density
        rho = -(c**2 / (8 * math.pi * G)) * b_prime / b0**2
        
        return rho
    
    def _minimum_velocity(self, wormhole: WormholeParams) -> float:
        """
        Compute minimum velocity needed to traverse wormhole.
        
        Based on the requirement that proper time remains real.
        """
        # Simplified: For Ellis wormhole, any velocity > 0 works
        return 0.001  # 0.1% of c
    
    # ═══════════════════════════════════════════════════════════════════════════
    # TRAJECTORY INTEGRATION
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _compute_trajectory(self, wormhole: WormholeParams, 
                           entry_velocity: float, entry_angle: float,
                           num_steps: int) -> WormholeTrajectory:
        """
        Integrate geodesic equations through the wormhole.
        
        Uses simplified radial geodesic for demo.
        Full implementation would use 4th-order Runge-Kutta.
        """
        c = CONSTANTS['c']
        b0 = wormhole.throat_radius
        
        # Initial conditions (far from throat)
        r0 = 10 * b0  # Start at 10 throat radii
        v0 = entry_velocity * c  # Convert to m/s
        
        # Energy and angular momentum (radial case: L=0)
        E = math.sqrt(1 + v0**2 / c**2)  # Specific energy
        L = 0.0  # Radial trajectory
        
        # Entry point
        entry = TrajectoryPoint(
            proper_time=0.0,
            coordinate_time=0.0,
            radial_position=r0,
            theta=math.pi / 2,  # Equatorial
            phi=entry_angle,
            velocity=v0,
            energy=E,
            angular_momentum=L,
        )
        
        # Integration
        points = [entry]
        tau = 0.0  # Proper time
        t = 0.0    # Coordinate time
        r = r0
        v = v0
        
        dt = (2 * r0 / v0) / num_steps  # Time step
        
        crossed_throat = False
        exotic_matter_total = 0.0
        
        for _ in range(num_steps):
            # Current metric coefficient
            g_rr = self._metric_coefficient(r, wormhole)
            
            # Proper time increment
            d_tau = dt * math.sqrt(1 - v**2 / c**2)
            tau += d_tau
            t += dt
            
            # Position update (simplified)
            dr = v * dt
            r -= dr  # Moving inward
            
            # Track throat crossing
            if r <= b0 and not crossed_throat:
                crossed_throat = True
                exotic_matter_total = abs(wormhole.exotic_matter_density * 4/3 * math.pi * b0**3)
            
            # After throat, r increases (emerging on other side)
            if crossed_throat and r <= b0:
                r = 2 * b0 - r  # Reflect through throat
                v = -v  # Velocity reverses direction (moving outward)
            
            # Stop when we've emerged far enough
            if crossed_throat and r >= 10 * b0:
                break
            
            # Record point
            points.append(TrajectoryPoint(
                proper_time=tau,
                coordinate_time=t,
                radial_position=r,
                theta=math.pi / 2,
                phi=entry_angle,
                velocity=v,
                energy=E,
                angular_momentum=L,
            ))
        
        # Exit point
        exit_point = points[-1] if points else entry
        
        return WormholeTrajectory(
            params=wormhole,
            entry_point=entry,
            exit_point=exit_point,
            points=points,
            total_proper_time=tau,
            total_coordinate_time=t,
            transit_successful=crossed_throat,
            exotic_matter_encountered=exotic_matter_total,
        )


# ═══════════════════════════════════════════════════════════════════════════════
# STANDALONE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Wormhole Traveler - Trajectory Calculator")
    print("=" * 50)
    
    calc = WormholeTrajectoryCalculator()
    calc.load()
    
    # Test trajectory
    print("\nComputing wormhole trajectory...")
    print("  Throat radius: 10 km")
    print("  Entry velocity: 1% of c")
    
    result = calc.execute({
        'throat_radius': 1e4,  # 10 km throat
        'entry_velocity': 0.01,  # 1% of c
        'entry_angle': 0.0,
        'num_steps': 1000,
        'shape_function': 'ellis',
    })
    
    if result['success']:
        traj = result['trajectory']
        wh = result['wormhole']
        
        print(f"\nTrajectory computed successfully!")
        print(f"  Points computed: {traj['num_points']}")
        print(f"  Total proper time: {traj['total_proper_time']:.4f} s")
        print(f"  Total coordinate time: {traj['total_coordinate_time']:.4f} s")
        print(f"  Transit successful: {traj['transit_successful']}")
        print(f"\nWormhole properties:")
        print(f"  Throat radius: {wh['throat_radius']:.0f} m")
        print(f"  Exotic matter density: {wh['exotic_matter_density']:.4e} kg/m³")
        print(f"  Is traversable: {wh['is_traversable']}")
        print(f"\nEntry point: r = {result['entry_point']['radial_position']:.0f} m")
        print(f"Exit point: r = {result['exit_point']['radial_position']:.0f} m")
    else:
        print(f"\nTrajectory failed: {result['error']}")
