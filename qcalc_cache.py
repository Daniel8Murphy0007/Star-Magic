#!/usr/bin/env python3
"""
qcalc_cache.py
==============
Phase 3: LRU caching for repeated UQFF queries

Features:
- LRU (Least Recently Used) cache with configurable size
- TTL (Time-To-Live) for cache entries
- Cache hit/miss tracking
- Query normalization for cache key generation
- Persistent cache (optional JSON file)

Author: Phase 3 Polish
Date: March 3, 2026
"""

import json
import time
import hashlib
from typing import Dict, Any, Optional, Tuple
from collections import OrderedDict
from dataclasses import dataclass, asdict
from datetime import datetime, timedelta

#══════════════════════════════════════════════════════════════════════════════
# CACHE CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CacheConfig:
    """Cache configuration"""
    max_size: int = 100                    # Maximum entries
    ttl_seconds: float = 3600.0            # 1 hour default
    enable_persistence: bool = False       # Save to disk
    persistence_file: str = "qcalc_cache.json"


# ═══════════════════════════════════════════════════════════════════════════════
# CACHE ENTRY
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CacheEntry:
    """Single cache entry with metadata"""
    key: str
    result: Dict[str, Any]
    created_at: float
    last_accessed: float
    access_count: int
    ttl: float
    
    def is_expired(self) -> bool:
        """Check if entry is expired"""
        return time.time() - self.created_at > self.ttl
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'key': self.key,
            'result': self.result,
            'created_at': self.created_at,
            'last_accessed': self.last_accessed,
            'access_count': self.access_count,
            'ttl': self.ttl
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'CacheEntry':
        """Create from dictionary"""
        return cls(**data)


# ═══════════════════════════════════════════════════════════════════════════════
# LRU CACHE
# ═══════════════════════════════════════════════════════════════════════════════

class LRUCache:
    """
    LRU cache with TTL support for UQFF calculations.
    
    Features:
    - Automatic eviction of least recently used entries
    - TTL expiration
    - Hit/miss statistics
    - Optional disk persistence
    """
    
    def __init__(self, config: Optional[CacheConfig] = None):
        """
        Initialize LRU cache.
        
        Args:
            config: Cache configuration
        """
        self.config = config or CacheConfig()
        self.cache: OrderedDict[str, CacheEntry] = OrderedDict()
        
        # Statistics
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self.expirations = 0
        
        # Load from disk if persistence enabled
        if self.config.enable_persistence:
            self._load_from_disk()
    
    def _generate_key(self, query_data: Dict[str, Any]) -> str:
        """
        Generate normalized cache key from query parameters.
        
        Args:
            query_data: Query parameters
            
        Returns:
            Cache key (MD5 hash of normalized params)
        """
        # Extract relevant parameters for caching
        cache_params = {
            'object_name': query_data.get('object_name', ''),
            'M': query_data.get('M', 0),
            'r': query_data.get('r', 0),
            'z': query_data.get('z', 0),
            'B': query_data.get('B', 0),
            'T': query_data.get('T', 0),
            'SFR': query_data.get('SFR', 0),
        }
        
        # Sort keys for consistency
        sorted_params = json.dumps(cache_params, sort_keys=True)
        
        # Generate MD5 hash
        return hashlib.md5(sorted_params.encode()).hexdigest()
    
    def get(self, query_data: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Retrieve cached result if available.
        
        Args:
            query_data: Query parameters
            
        Returns:
            Cached result or None if not found/expired
        """
        key = self._generate_key(query_data)
        
        # Check if key exists
        if key not in self.cache:
            self.misses += 1
            return None
        
        entry = self.cache[key]
        
        # Check if expired
        if entry.is_expired():
            del self.cache[key]
            self.expirations += 1
            self.misses += 1
            return None
        
        # Update access metadata
        entry.last_accessed = time.time()
        entry.access_count += 1
        
        # Move to end (most recently used)
        self.cache.move_to_end(key)
        
        self.hits += 1
        
        # Add cache metadata to result
        result = entry.result.copy()
        result['__cache_hit__'] = True
        result['__cache_age_seconds__'] = time.time() - entry.created_at
        result['__cache_access_count__'] = entry.access_count
        
        return result
    
    def put(self, query_data: Dict[str, Any], result: Dict[str, Any], 
            ttl: Optional[float] = None):
        """
        Store result in cache.
        
        Args:
            query_data: Query parameters
            result: Calculation result
            ttl: Optional custom TTL (uses config default if None)
        """
        key = self._generate_key(query_data)
        
        # Use custom TTL or default
        entry_ttl = ttl if ttl is not None else self.config.ttl_seconds
        
        # Create entry
        entry = CacheEntry(
            key=key,
            result=result.copy(),
            created_at=time.time(),
            last_accessed=time.time(),
            access_count=0,
            ttl=entry_ttl
        )
        
        # Check if we need to evict
        if len(self.cache) >= self.config.max_size:
            # Remove least recently used (first item)
            self.cache.popitem(last=False)
            self.evictions += 1
        
        # Add to cache (end = most recently used)
        self.cache[key] = entry
        
        # Persist if enabled
        if self.config.enable_persistence:
            self._save_to_disk()
    
    def invalidate(self, query_data: Dict[str, Any]) -> bool:
        """
        Invalidate specific cache entry.
        
        Args:
            query_data: Query parameters
            
        Returns:
            True if entry was found and removed
        """
        key = self._generate_key(query_data)
        
        if key in self.cache:
            del self.cache[key]
            return True
        
        return False
    
    def clear(self):
        """Clear entire cache"""
        self.cache.clear()
        self.hits = 0
        self.misses = 0
        self.evictions = 0
        self.expirations = 0
        
        if self.config.enable_persistence:
            self._save_to_disk()
    
    def get_stats(self) -> Dict[str, Any]:
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache metrics
        """
        total_requests = self.hits + self.misses
        hit_rate = (self.hits / total_requests * 100) if total_requests > 0 else 0.0
        
        # Calculate average age of entries
        if self.cache:
            avg_age = sum(time.time() - e.created_at for e in self.cache.values()) / len(self.cache)
        else:
            avg_age = 0.0
        
        return {
            'size': len(self.cache),
            'max_size': self.config.max_size,
            'hits': self.hits,
            'misses': self.misses,
            'hit_rate_percent': round(hit_rate, 2),
            'evictions': self.evictions,
            'expirations': self.expirations,
            'total_requests': total_requests,
            'avg_entry_age_seconds': round(avg_age, 1),
            'ttl_seconds': self.config.ttl_seconds
        }
    
    def _save_to_disk(self):
        """Save cache to disk (JSON)"""
        try:
            cache_data = {
                'config': asdict(self.config),
                'stats': {
                    'hits': self.hits,
                    'misses': self.misses,
                    'evictions': self.evictions,
                    'expirations': self.expirations
                },
                'entries': [entry.to_dict() for entry in self.cache.values()]
            }
            
            with open(self.config.persistence_file, 'w') as f:
                json.dump(cache_data, f, indent=2)
                
        except Exception as e:
            print(f"[Cache] Failed to save to disk: {str(e)}", file=sys.stderr)
    
    def _load_from_disk(self):
        """Load cache from disk (JSON)"""
        try:
            with open(self.config.persistence_file, 'r') as f:
                cache_data = json.load(f)
            
            # Restore stats
            stats = cache_data.get('stats', {})
            self.hits = stats.get('hits', 0)
            self.misses = stats.get('misses', 0)
            self.evictions = stats.get('evictions', 0)
            self.expirations = stats.get('expirations', 0)
            
            # Restore entries (skip expired)
            for entry_data in cache_data.get('entries', []):
                entry = CacheEntry.from_dict(entry_data)
                if not entry.is_expired():
                    self.cache[entry.key] = entry
                else:
                    self.expirations += 1
            
            print(f"[Cache] Loaded {len(self.cache)} entries from disk", file=sys.stderr)
            
        except FileNotFoundError:
            print(f"[Cache] No cache file found, starting fresh", file=sys.stderr)
        except Exception as e:
            print(f"[Cache] Failed to load from disk: {str(e)}", file=sys.stderr)


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL CACHE INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

# Global cache with default configuration
GLOBAL_CACHE = LRUCache(CacheConfig(
    max_size=100,              # Cache up to 100 queries
    ttl_seconds=3600.0,        # 1 hour TTL
    enable_persistence=False   # Disabled by default (enable for production)
))


if __name__ == '__main__':
    import sys
    
    print("Testing LRU cache...")
    
    # Create test cache
    cache = LRUCache(CacheConfig(max_size=3, ttl_seconds=2.0))
    
    # Test data
    query1 = {'object_name': 'Sun', 'M': 1.989e30, 'r': 1.496e11}
    result1 = {'F_U': 5.93e-3, 'success': True}
    
    query2 = {'object_name': 'Sgr A*', 'M': 4.15e6 * 1.989e30, 'r': 8.5e3 * 3.086e16}
    result2 = {'F_U': 1.23e-5, 'success': True}
    
    # Test cache miss
    print("\n1. Cache miss test:")
    cached = cache.get(query1)
    print(f"   Result: {cached} (expected: None)")
    
    # Test cache put and hit
    print("\n2. Cache put and hit test:")
    cache.put(query1, result1)
    cached = cache.get(query1)
    print(f"   Result: {cached is not None} (expected: True)")
    print(f"   Cache hit: {cached.get('__cache_hit__')} (expected: True)")
    
    # Test LRU eviction
    print("\n3. LRU eviction test (max_size=3):")
    cache.put(query2, result2)
    cache.put({'object_name': 'Query3', 'M': 1, 'r': 1}, {'F_U': 1})
    cache.put({'object_name': 'Query4', 'M': 1, 'r': 1}, {'F_U': 1})  # Should evict query1
    print(f"   Cache size: {len(cache.cache)} (expected: 3)")
    print(f"   Query1 still cached: {cache.get(query1) is not None} (expected: False)")
    
    # Test TTL expiration
    print("\n4. TTL expiration test (ttl=2s):")
    cache_ttl = LRUCache(CacheConfig(max_size=10, ttl_seconds=1.0))
    cache_ttl.put(query1, result1)
    print("   Sleeping 1.5s...")
    time.sleep(1.5)
    expired = cache_ttl.get(query1)
    print(f"   Result after TTL: {expired} (expected: None)")
    
    # Test statistics
    print("\n5. Cache statistics:")
    stats = cache.get_stats()
    print(f"   Hit rate: {stats['hit_rate_percent']}%")
    print(f"   Total requests: {stats['total_requests']}")
    print(f"   Evictions: {stats['evictions']}")
    
    print("\n✓ Cache module ready for Phase 3 integration!")
