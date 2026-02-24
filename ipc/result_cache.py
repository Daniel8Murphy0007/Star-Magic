"""
result_cache.py - LRU Result Cache for UQFF Physics Calculations
=================================================================

Caches physics calculation results to avoid redundant computation.
Uses content-based hashing to detect identical queries.

Storage:
- In-memory LRU cache (fast)
- SQLite persistent cache (durable)
- File-based JSON cache (debuggable)

Usage:
    from ipc.result_cache import cache_physics_result, get_cached_result
    
    # Check cache before computing
    cached = get_cached_result(params)
    if cached:
        return cached
    
    # Compute and cache
    result = heavy_physics_calculation(params)
    cache_physics_result(params, result)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

import os
import json
import hashlib
import sqlite3
import logging
import threading
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
from datetime import datetime, timedelta
from collections import OrderedDict
from dataclasses import dataclass, asdict

logger = logging.getLogger("UQFF_ResultCache")


def _hash_params(params: Dict[str, Any]) -> str:
    """Create deterministic hash of parameters"""
    # Sort keys for consistency
    canonical = json.dumps(params, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


@dataclass
class CacheEntry:
    """Cached result entry"""
    key: str
    params_hash: str
    result: Dict[str, Any]
    created_at: datetime
    access_count: int = 0
    last_accessed: datetime = None
    
    def __post_init__(self):
        if self.last_accessed is None:
            self.last_accessed = self.created_at


class ResultCache:
    """
    Multi-tier result cache with LRU eviction
    
    Tiers:
    1. In-memory LRU (fastest, volatile)
    2. SQLite (fast, persistent)
    3. JSON files (debuggable, persistent)
    """
    
    def __init__(self,
                 cache_dir: str = ".cache",
                 max_memory_entries: int = 1000,
                 max_db_entries: int = 10000,
                 ttl_hours: int = 24):
        """
        Initialize ResultCache
        
        Args:
            cache_dir: Directory for persistent cache
            max_memory_entries: Max entries in memory LRU
            max_db_entries: Max entries in SQLite database
            ttl_hours: Time-to-live in hours
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        
        self.max_memory = max_memory_entries
        self.max_db = max_db_entries
        self.ttl = timedelta(hours=ttl_hours)
        
        # In-memory LRU cache
        self._memory_cache: OrderedDict[str, CacheEntry] = OrderedDict()
        self._lock = threading.RLock()
        
        # SQLite database
        self._db_path = self.cache_dir / "result_cache.db"
        self._init_database()
        
        # Statistics
        self.hits = 0
        self.misses = 0
    
    def _init_database(self) -> None:
        """Initialize SQLite database"""
        conn = sqlite3.connect(str(self._db_path))
        conn.execute("""
            CREATE TABLE IF NOT EXISTS cache (
                key TEXT PRIMARY KEY,
                params_hash TEXT NOT NULL,
                result_json TEXT NOT NULL,
                created_at TEXT NOT NULL,
                last_accessed TEXT NOT NULL,
                access_count INTEGER DEFAULT 0
            )
        """)
        conn.execute("CREATE INDEX IF NOT EXISTS idx_hash ON cache(params_hash)")
        conn.execute("CREATE INDEX IF NOT EXISTS idx_accessed ON cache(last_accessed)")
        conn.commit()
        conn.close()
    
    def _evict_memory_if_needed(self) -> None:
        """Evict oldest entries if memory limit exceeded"""
        while len(self._memory_cache) >= self.max_memory:
            # Remove oldest (LRU)
            self._memory_cache.popitem(last=False)
    
    def _evict_db_if_needed(self, conn: sqlite3.Connection) -> None:
        """Evict oldest entries from database if limit exceeded"""
        cursor = conn.execute("SELECT COUNT(*) FROM cache")
        count = cursor.fetchone()[0]
        
        if count >= self.max_db:
            # Delete oldest 10%
            delete_count = max(100, self.max_db // 10)
            conn.execute(f"""
                DELETE FROM cache WHERE key IN (
                    SELECT key FROM cache 
                    ORDER BY last_accessed ASC 
                    LIMIT {delete_count}
                )
            """)
            conn.commit()
    
    def get(self, params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Get cached result for parameters
        
        Args:
            params: Physics calculation parameters
            
        Returns:
            Cached result or None if not found
        """
        params_hash = _hash_params(params)
        
        with self._lock:
            # Check memory cache first
            if params_hash in self._memory_cache:
                entry = self._memory_cache[params_hash]
                
                # Check TTL
                if datetime.now() - entry.created_at > self.ttl:
                    del self._memory_cache[params_hash]
                else:
                    # Update access info
                    entry.access_count += 1
                    entry.last_accessed = datetime.now()
                    # Move to end (most recently used)
                    self._memory_cache.move_to_end(params_hash)
                    self.hits += 1
                    return entry.result
        
        # Check SQLite
        try:
            conn = sqlite3.connect(str(self._db_path))
            cursor = conn.execute(
                "SELECT result_json, created_at, access_count FROM cache WHERE params_hash = ?",
                (params_hash,)
            )
            row = cursor.fetchone()
            
            if row:
                created_at = datetime.fromisoformat(row[1])
                
                # Check TTL
                if datetime.now() - created_at > self.ttl:
                    conn.execute("DELETE FROM cache WHERE params_hash = ?", (params_hash,))
                    conn.commit()
                else:
                    # Update access info
                    conn.execute("""
                        UPDATE cache SET 
                            access_count = access_count + 1,
                            last_accessed = ?
                        WHERE params_hash = ?
                    """, (datetime.now().isoformat(), params_hash))
                    conn.commit()
                    
                    result = json.loads(row[0])
                    
                    # Promote to memory cache
                    with self._lock:
                        self._evict_memory_if_needed()
                        self._memory_cache[params_hash] = CacheEntry(
                            key=params_hash,
                            params_hash=params_hash,
                            result=result,
                            created_at=created_at,
                            access_count=row[2] + 1
                        )
                    
                    self.hits += 1
                    conn.close()
                    return result
            
            conn.close()
        except Exception as e:
            logger.error(f"SQLite cache error: {e}")
        
        self.misses += 1
        return None
    
    def put(self, params: Dict[str, Any], result: Dict[str, Any]) -> None:
        """
        Store result in cache
        
        Args:
            params: Physics calculation parameters
            result: Calculation result to cache
        """
        params_hash = _hash_params(params)
        now = datetime.now()
        
        # Store in memory
        with self._lock:
            self._evict_memory_if_needed()
            self._memory_cache[params_hash] = CacheEntry(
                key=params_hash,
                params_hash=params_hash,
                result=result,
                created_at=now
            )
        
        # Store in SQLite
        try:
            conn = sqlite3.connect(str(self._db_path))
            self._evict_db_if_needed(conn)
            
            conn.execute("""
                INSERT OR REPLACE INTO cache 
                (key, params_hash, result_json, created_at, last_accessed, access_count)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                params_hash,
                params_hash,
                json.dumps(result, default=str),
                now.isoformat(),
                now.isoformat(),
                0
            ))
            conn.commit()
            conn.close()
        except Exception as e:
            logger.error(f"SQLite cache write error: {e}")
    
    def invalidate(self, params: Dict[str, Any] = None) -> None:
        """
        Invalidate cache entries
        
        Args:
            params: If provided, invalidate specific entry; else clear all
        """
        if params:
            params_hash = _hash_params(params)
            
            with self._lock:
                if params_hash in self._memory_cache:
                    del self._memory_cache[params_hash]
            
            try:
                conn = sqlite3.connect(str(self._db_path))
                conn.execute("DELETE FROM cache WHERE params_hash = ?", (params_hash,))
                conn.commit()
                conn.close()
            except Exception as e:
                logger.error(f"Cache invalidation error: {e}")
        else:
            # Clear all
            with self._lock:
                self._memory_cache.clear()
            
            try:
                conn = sqlite3.connect(str(self._db_path))
                conn.execute("DELETE FROM cache")
                conn.commit()
                conn.close()
            except Exception as e:
                logger.error(f"Cache clear error: {e}")
    
    def stats(self) -> Dict[str, Any]:
        """Get cache statistics"""
        total = self.hits + self.misses
        hit_rate = self.hits / total if total > 0 else 0
        
        try:
            conn = sqlite3.connect(str(self._db_path))
            cursor = conn.execute("SELECT COUNT(*) FROM cache")
            db_count = cursor.fetchone()[0]
            conn.close()
        except:
            db_count = 0
        
        return {
            'memory_entries': len(self._memory_cache),
            'db_entries': db_count,
            'hits': self.hits,
            'misses': self.misses,
            'hit_rate': hit_rate,
            'max_memory': self.max_memory,
            'max_db': self.max_db,
            'ttl_hours': self.ttl.total_seconds() / 3600
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE-LEVEL CACHE INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════
_result_cache: Optional[ResultCache] = None


def get_result_cache() -> ResultCache:
    """Get or create the module-level result cache"""
    global _result_cache
    if _result_cache is None:
        _result_cache = ResultCache()
    return _result_cache


def cache_physics_result(params: Dict[str, Any], result: Dict[str, Any]) -> None:
    """Cache a physics calculation result"""
    get_result_cache().put(params, result)


def get_cached_result(params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Get a cached physics result"""
    return get_result_cache().get(params)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("ResultCache Test")
    print("=" * 40)
    
    cache = ResultCache(cache_dir=".test_cache")
    
    # Test caching
    params = {'body': 'NGC3596', 'M': 1e11, 'r': 1e20}
    result = {'F_U_Bi_i': 1.83e71, 'timestamp': datetime.now().isoformat()}
    
    # Miss
    cached = cache.get(params)
    print(f"First get (miss): {cached}")
    
    # Store
    cache.put(params, result)
    print(f"Stored result")
    
    # Hit
    cached = cache.get(params)
    print(f"Second get (hit): {cached is not None}")
    
    # Stats
    print(f"\nStats: {cache.stats()}")
    
    print("\nTest complete")
