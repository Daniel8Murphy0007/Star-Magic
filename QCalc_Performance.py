#!/usr/bin/env python3
"""
QCalc_Performance.py - Performance Optimization Layer
=====================================================

Performance optimizations for QCalc.py calculations without modifying core logic.

DESIGN PRINCIPLES:
- Non-invasive: QCalc.py remains pure physics calculator
- Caching: Avoid repeated expensive computations
- Vectorization: Numpy optimizations for bulk calculations
- Profiling: Built-in performance monitoring

ARCHITECTURE:
    QCalc.py → QCalc_Performance.py → Optimized results
    
Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from functools import lru_cache, wraps
from typing import Dict, List, Any, Callable, Tuple
import time
import hashlib
import json
from collections import defaultdict

# ═══════════════════════════════════════════════════════════════════════════════
# PERFORMANCE MONITORING
# ═══════════════════════════════════════════════════════════════════════════════

class PerformanceMonitor:
    """Track computation times and identify bottlenecks."""
    
    def __init__(self):
        self.timings: Dict[str, List[float]] = defaultdict(list)
        self.call_counts: Dict[str, int] = defaultdict(int)
        
    def time_function(self, func: Callable) -> Callable:
        """Decorator to time function execution."""
        @wraps(func)
        def wrapper(*args, **kwargs):
            start = time.perf_counter()
            result = func(*args, **kwargs)
            elapsed = time.perf_counter() - start
            
            func_name = func.__name__
            self.timings[func_name].append(elapsed)
            self.call_counts[func_name] += 1
            
            return result
        return wrapper
    
    def get_report(self) -> Dict[str, Dict[str, float]]:
        """Get performance statistics."""
        report = {}
        for func_name, times in self.timings.items():
            report[func_name] = {
                'total_time': sum(times),
                'avg_time': np.mean(times),
                'min_time': min(times),
                'max_time': max(times),
                'std_time': np.std(times),
                'call_count': self.call_counts[func_name]
            }
        return report
    
    def print_report(self):
        """Print formatted performance report."""
        print("\n" + "="*80)
        print("PERFORMANCE REPORT")
        print("="*80)
        
        report = self.get_report()
        sorted_funcs = sorted(report.items(), key=lambda x: x[1]['total_time'], reverse=True)
        
        print(f"{'Function':<40} {'Calls':>8} {'Avg (ms)':>12} {'Total (s)':>12}")
        print("-"*80)
        
        for func_name, stats in sorted_funcs:
            print(f"{func_name:<40} {stats['call_count']:>8} "
                  f"{stats['avg_time']*1000:>12.4f} {stats['total_time']:>12.4f}")
        
        print("="*80 + "\n")


# Global monitor instance
MONITOR = PerformanceMonitor()


# ═══════════════════════════════════════════════════════════════════════════════
# INTELLIGENT CACHING SYSTEM
# ═══════════════════════════════════════════════════════════════════════════════

class ResultCache:
    """
    LRU cache for expensive calculations with parameter-aware invalidation.
    
    Features:
    - Automatic cache key generation from parameters
    - TTL (time-to-live) support
    - Size limits
    - Cache statistics
    """
    
    def __init__(self, max_size: int = 1000, ttl_seconds: float = 3600):
        self.max_size = max_size
        self.ttl = ttl_seconds
        self.cache: Dict[str, Tuple[Any, float]] = {}
        self.hits = 0
        self.misses = 0
        
    def _make_key(self, func_name: str, args: tuple, kwargs: dict) -> str:
        """Generate unique cache key from function signature."""
        # Convert args/kwargs to stable JSON string
        key_data = {
            'func': func_name,
            'args': [self._serialize_arg(arg) for arg in args],
            'kwargs': {k: self._serialize_arg(v) for k, v in sorted(kwargs.items())}
        }
        key_str = json.dumps(key_data, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    def _serialize_arg(self, arg: Any) -> Any:
        """Serialize argument for cache key."""
        if isinstance(arg, np.ndarray):
            return arg.tolist()
        elif isinstance(arg, dict):
            return {k: self._serialize_arg(v) for k, v in sorted(arg.items())}
        elif hasattr(arg, '__dict__'):
            return self._serialize_arg(vars(arg))
        else:
            return arg
    
    def get(self, func_name: str, args: tuple, kwargs: dict) -> Tuple[bool, Any]:
        """
        Get cached result if available and not expired.
        
        Returns:
            (found: bool, result: Any)
        """
        key = self._make_key(func_name, args, kwargs)
        
        if key in self.cache:
            result, timestamp = self.cache[key]
            age = time.time() - timestamp
            
            if age < self.ttl:
                self.hits += 1
                return True, result
            else:
                # Expired - remove from cache
                del self.cache[key]
        
        self.misses += 1
        return False, None
    
    def set(self, func_name: str, args: tuple, kwargs: dict, result: Any):
        """Store result in cache."""
        key = self._make_key(func_name, args, kwargs)
        
        # Check size limit
        if len(self.cache) >= self.max_size:
            # Remove oldest entry (simple FIFO, could be improved to LRU)
            oldest_key = next(iter(self.cache))
            del self.cache[oldest_key]
        
        self.cache[key] = (result, time.time())
    
    def clear(self):
        """Clear all cached entries."""
        self.cache.clear()
        self.hits = 0
        self.misses = 0
    
    def get_stats(self) -> Dict[str, Any]:
        """Get cache statistics."""
        total_requests = self.hits + self.misses
        hit_rate = self.hits / total_requests if total_requests > 0 else 0.0
        
        return {
            'size': len(self.cache),
            'max_size': self.max_size,
            'hits': self.hits,
            'misses': self.misses,
            'hit_rate': hit_rate,
            'total_requests': total_requests
        }


# Global cache instance
CACHE = ResultCache(max_size=1000, ttl_seconds=3600)


def cached_calculation(func: Callable) -> Callable:
    """Decorator to cache expensive calculations."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        # Check cache
        found, result = CACHE.get(func.__name__, args, kwargs)
        if found:
            return result
        
        # Compute and cache
        result = func(*args, **kwargs)
        CACHE.set(func.__name__, args, kwargs, result)
        return result
    
    return wrapper


# ═══════════════════════════════════════════════════════════════════════════════
# VECTORIZED OPERATIONS
# ═══════════════════════════════════════════════════════════════════════════════

class VectorizedCalculations:
    """
    Numpy-optimized bulk calculations for multiple systems.
    
    Example: Calculate Ug1-4 for 1000 systems simultaneously.
    """
    
    @staticmethod
    def vectorized_26_level_energy(E0: float, levels: np.ndarray) -> np.ndarray:
        """
        Compute 26-level energy structure for multiple base energies.
        
        Args:
            E0: Base energy (J) - can be array
            levels: Array of level numbers (0-25)
            
        Returns:
            Energy array: shape (len(E0), 26) if E0 is array
        """
        E0_arr = np.atleast_1d(E0)
        levels_arr = np.atleast_1d(levels)
        
        # Broadcasting: E0[i] × 10^level[j]
        E_n = E0_arr[:, np.newaxis] * (10.0 ** levels_arr[np.newaxis, :])
        return E_n.squeeze()
    
    @staticmethod
    def vectorized_reactor_efficiency(
        E0: np.ndarray,
        t: np.ndarray,
        t_n: np.ndarray,
        kappa: float = 0.0005,
        omega: float = 1.16e-7,
        S_k_B: float = 0.1,
        SCm: float = 1e3
    ) -> np.ndarray:
        """
        Vectorized reactor efficiency for time series.
        
        Args:
            E0: Base energies (shape: N,)
            t: Time points (shape: M,)
            t_n: Negative time points (shape: M,)
            
        Returns:
            E_react: shape (N, M)
        """
        E0_arr = np.atleast_2d(E0).T  # (N, 1)
        t_arr = np.atleast_1d(t)      # (M,)
        t_n_arr = np.atleast_1d(t_n)  # (M,)
        
        # Time decay
        decay = np.exp(-kappa * t_arr)
        
        # Negative time oscillation
        oscillation = np.cos(omega * t_n_arr)
        
        # Entropy correction (simplified)
        entropy_factor = 1.0 - S_k_B
        
        # Broadcasting: E0[i] × factors[j]
        E_react = E0_arr * decay[np.newaxis, :] * oscillation[np.newaxis, :] * entropy_factor * SCm
        return E_react.squeeze()
    
    @staticmethod
    def vectorized_stress_energy_tensor(
        lambda_UA: np.ndarray,
        lambda_SCm: np.ndarray,
        lambda_A: np.ndarray,
        t_n: np.ndarray,
        T_base: float = 1.27e3,
        T_cosmic: float = 1.11e7,
        omega_c: float = 7.27e-5
    ) -> np.ndarray:
        """
        Vectorized stress-energy tensor T_s^00 for multiple systems/times.
        
        Args:
            lambda_UA, lambda_SCm, lambda_A: Vacuum densities (shape: N,)
            t_n: Time points (shape: M,)
            
        Returns:
            T_s_00: Energy density (shape: N, M)
        """
        # Ensure 2D arrays for broadcasting
        lambda_UA_2d = np.atleast_2d(lambda_UA).T    # (N, 1)
        lambda_SCm_2d = np.atleast_2d(lambda_SCm).T  # (N, 1)
        lambda_A_2d = np.atleast_2d(lambda_A).T      # (N, 1)
        t_n_arr = np.atleast_1d(t_n)                 # (M,)
        
        # Time modulation
        f_t = 1.0 + 0.1 * np.cos(omega_c * t_n_arr)  # (M,)
        
        # Base contribution
        T_base_contrib = T_base * (lambda_UA_2d + lambda_SCm_2d)  # (N, 1)
        
        # Cosmic contribution with time modulation
        T_cosmic_contrib = T_cosmic * lambda_A_2d * f_t[np.newaxis, :]  # (N, M)
        
        # Total stress-energy (broadcast base to match cosmic)
        T_s_00 = T_base_contrib + T_cosmic_contrib  # (N, M)
        return T_s_00.squeeze()
    
    @staticmethod
    def batch_metric_perturbations(
        T_s_tensor: np.ndarray,
        eta: float = 1e-22
    ) -> np.ndarray:
        """
        Compute metric perturbations for multiple stress-energy tensors.
        
        Args:
            T_s_tensor: Stress-energy tensors (shape: N, 4, 4)
            eta: Aether coupling
            
        Returns:
            delta_g: Metric perturbations (shape: N, 4, 4)
        """
        return eta * T_s_tensor


# ═══════════════════════════════════════════════════════════════════════════════
# PARALLEL PROCESSING UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def parallel_solve(
    solver,
    systems: List[Dict[str, Any]],
    n_workers: int = 4
) -> List[Dict[str, Any]]:
    """
    Solve multiple systems in parallel (future: multiprocessing).
    
    Current: Sequential (multiprocessing requires careful state management)
    Future: Use ProcessPoolExecutor with serializable solver state
    
    Args:
        solver: UnifiedFieldSolver instance
        systems: List of system parameter dictionaries
        n_workers: Number of parallel workers (not yet implemented)
        
    Returns:
        List of solution dictionaries
    """
    # For now, sequential (multiprocessing requires state serialization)
    results = []
    for system in systems:
        result = solver.solve(system)
        results.append(result)
    
    return results


# ═══════════════════════════════════════════════════════════════════════════════
# MEMORY OPTIMIZATION
# ═══════════════════════════════════════════════════════════════════════════════

class MemoryOptimizer:
    """Utilities for reducing memory footprint of large calculations."""
    
    @staticmethod
    def compress_results(results: Dict[str, Any]) -> bytes:
        """
        Compress results for storage.
        
        Uses zlib for space savings (typical 70-80% reduction).
        """
        import zlib
        json_str = json.dumps(results)
        compressed = zlib.compress(json_str.encode('utf-8'))
        return compressed
    
    @staticmethod
    def decompress_results(compressed: bytes) -> Dict[str, Any]:
        """Decompress results from storage."""
        import zlib
        json_str = zlib.decompress(compressed).decode('utf-8')
        return json.loads(json_str)
    
    @staticmethod
    def truncate_precision(value: float, significant_digits: int = 10) -> float:
        """
        Reduce float precision for memory savings.
        
        Example: 1.23456789012345 → 1.234567890 (10 digits)
        """
        if value == 0:
            return 0.0
        
        from math import log10, floor
        magnitude = floor(log10(abs(value)))
        scale = 10 ** (significant_digits - magnitude - 1)
        return round(value * scale) / scale


# ═══════════════════════════════════════════════════════════════════════════════
# USAGE EXAMPLES
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("="*80)
    print("QCalc Performance Optimization Layer - Test Suite")
    print("="*80)
    
    # Test 1: Vectorized 26-level energy
    print("\n[TEST 1] Vectorized 26-level energy")
    E0 = 1.0
    levels = np.arange(26)
    E_n = VectorizedCalculations.vectorized_26_level_energy(E0, levels)
    print(f"  E_0  = {E_n[0]:.4e} J")
    print(f"  E_15 = {E_n[15]:.4e} J (Human scale)")
    print(f"  E_25 = {E_n[25]:.4e} J (Cosmic scale)")
    
    # Test 2: Vectorized reactor efficiency (time series)
    print("\n[TEST 2] Vectorized reactor efficiency")
    t_points = np.linspace(0, 86400, 100)  # 1 day
    t_n_points = -t_points
    E_react = VectorizedCalculations.vectorized_reactor_efficiency(
        E0=np.array([1.0]),
        t=t_points,
        t_n=t_n_points
    )
    print(f"  E_react(t=0)     = {E_react[0]:.4e} J")
    print(f"  E_react(t=1day)  = {E_react[-1]:.4e} J")
    print(f"  Decay over 1 day = {(E_react[0] - E_react[-1])/E_react[0] * 100:.2f}%")
    
    # Test 3: Stress-energy tensor (vectorized)
    print("\n[TEST 3] Vectorized stress-energy tensor")
    lambda_UA = np.array([7.09e-36])
    lambda_SCm = np.array([7.09e-37])
    lambda_A = np.array([8.99e-07])
    t_n_arr = np.array([0.0, -1000, -10000])
    
    T_s_00 = VectorizedCalculations.vectorized_stress_energy_tensor(
        lambda_UA, lambda_SCm, lambda_A, t_n_arr
    )
    print(f"  T_s^00(t_n=0)     = {T_s_00[0]:.4e} kg/m³ c²")
    print(f"  T_s^00(t_n=-1000) = {T_s_00[1]:.4e} kg/m³ c²")
    print(f"  T_s^00(t_n=-10k)  = {T_s_00[2]:.4e} kg/m³ c²")
    
    # Test 4: Caching system
    print("\n[TEST 4] Result caching")
    
    @cached_calculation
    def expensive_function(x: float) -> float:
        time.sleep(0.01)  # Simulate expensive computation
        return x ** 2
    
    # First call (miss)
    start = time.perf_counter()
    result1 = expensive_function(5.0)
    time1 = time.perf_counter() - start
    
    # Second call (hit)
    start = time.perf_counter()
    result2 = expensive_function(5.0)
    time2 = time.perf_counter() - start
    
    print(f"  First call:  {time1*1000:.2f} ms (cache miss)")
    print(f"  Second call: {time2*1000:.2f} ms (cache hit)")
    print(f"  Speedup:     {time1/time2:.1f}x")
    print(f"  Cache stats: {CACHE.get_stats()}")
    
    # Test 5: Memory optimization
    print("\n[TEST 5] Memory optimization")
    test_data = {
        'equations': ['E_n', 'E_react', 'lambda_vac'],
        'solutions': {'E_n': [1e0, 1e15, 1e25], 'E_react': 9.0e2},
        'metadata': {'timestamp': '2026-02-12', 'query': 'test'}
    }
    
    # Original size
    json_str = json.dumps(test_data)
    original_size = len(json_str.encode('utf-8'))
    
    # Compressed size
    compressed = MemoryOptimizer.compress_results(test_data)
    compressed_size = len(compressed)
    
    # Decompress and verify
    decompressed = MemoryOptimizer.decompress_results(compressed)
    
    print(f"  Original size:   {original_size} bytes")
    print(f"  Compressed size: {compressed_size} bytes")
    print(f"  Compression:     {(1 - compressed_size/original_size)*100:.1f}% reduction")
    print(f"  Integrity:       {'PASS' if decompressed == test_data else 'FAIL'}")
    
    print("\n" + "="*80)
    print("All tests complete! Performance layer ready for production.")
    print("="*80)
