"""
state_sync.py - Multi-Program State Synchronization
====================================================

Maintains consistent state across the 5 Principal Programs using
version vectors and atomic file operations.

State includes:
- Current input parameters (from APIFetch.py)
- Active calculations in progress
- Last results from each calculator
- Configuration/settings

Usage:
    from ipc.state_sync import sync_state, get_state_synchronizer
    
    # Get synchronized state
    state = sync_state()
    
    # Update local state
    state.update('QCalc', {'last_result': result})
    
    # Get state from other programs
    coangqi_state = state.get_program_state('MAIN_1_CoAnQi')

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

import os
import sys
import json
import logging
import threading
import time
from pathlib import Path
from typing import Dict, Any, Optional, List
from datetime import datetime
from dataclasses import dataclass, field, asdict
from copy import deepcopy

logger = logging.getLogger("UQFF_StateSync")


@dataclass
class ProgramState:
    """State for a single program"""
    name: str
    status: str = "unknown"  # 'running', 'idle', 'error', 'unknown'
    last_heartbeat: str = ""
    last_result: Dict[str, Any] = field(default_factory=dict)
    current_params: Dict[str, Any] = field(default_factory=dict)
    version: int = 0
    
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'ProgramState':
        return cls(**data)


@dataclass
class StateSnapshot:
    """Complete state snapshot across all programs"""
    timestamp: str = ""
    version: int = 0
    programs: Dict[str, ProgramState] = field(default_factory=dict)
    shared: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        if not self.timestamp:
            self.timestamp = datetime.now().isoformat()
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            'timestamp': self.timestamp,
            'version': self.version,
            'programs': {k: v.to_dict() for k, v in self.programs.items()},
            'shared': self.shared
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'StateSnapshot':
        programs = {}
        for name, pdata in data.get('programs', {}).items():
            programs[name] = ProgramState.from_dict(pdata)
        
        return cls(
            timestamp=data.get('timestamp', ''),
            version=data.get('version', 0),
            programs=programs,
            shared=data.get('shared', {})
        )


class StateSynchronizer:
    """
    Coordinates state across multiple programs using file-based sync
    
    Uses atomic write-then-rename pattern for corruption safety.
    Version vectors detect conflicts and enable merge.
    """
    
    PROGRAM_NAMES = [
        'source2',           # Principal GUI
        'QCalc',            # Pure Physics Solver
        'CondensedPhysics', # Calculator
        'MAIN_1_CoAnQi',    # Physics Engine
        'vr_runtime'        # VR/VM Simulator
    ]
    
    def __init__(self,
                 state_dir: str = ".state",
                 program_name: str = "Python",
                 sync_interval: float = 5.0):
        """
        Initialize StateSynchronizer
        
        Args:
            state_dir: Directory for state files
            program_name: This program's identifier
            sync_interval: Seconds between sync checks
        """
        self.state_dir = Path(state_dir)
        self.state_dir.mkdir(exist_ok=True)
        
        self.program_name = program_name
        self.sync_interval = sync_interval
        
        self._state_file = self.state_dir / "state_snapshot.json"
        self._lock_file = self.state_dir / "state.lock"
        
        self._local_state: Optional[StateSnapshot] = None
        self._lock = threading.RLock()
        
        self._running = False
        self._sync_thread: Optional[threading.Thread] = None
        
        # Callbacks
        self.on_state_changed: Optional[callable] = None
        
        # Load or create initial state
        self._load_state()
    
    def _load_state(self) -> None:
        """Load state from file or create new"""
        try:
            if self._state_file.exists():
                with open(self._state_file, 'r') as f:
                    data = json.load(f)
                self._local_state = StateSnapshot.from_dict(data)
                logger.info(f"Loaded state v{self._local_state.version}")
            else:
                self._local_state = StateSnapshot()
                self._save_state()
                logger.info("Created new state file")
        except Exception as e:
            logger.error(f"State load error: {e}")
            self._local_state = StateSnapshot()
    
    def _save_state(self) -> None:
        """Save state with atomic write"""
        try:
            # Write to temp file first
            tmp_file = self._state_file.with_suffix('.tmp')
            with open(tmp_file, 'w') as f:
                json.dump(self._local_state.to_dict(), f, indent=2, default=str)
            
            # Atomic rename
            tmp_file.replace(self._state_file)
            logger.debug(f"Saved state v{self._local_state.version}")
        except Exception as e:
            logger.error(f"State save error: {e}")
    
    def _acquire_lock(self, timeout: float = 5.0) -> bool:
        """Acquire file lock for atomic operations"""
        start = time.time()
        while time.time() - start < timeout:
            try:
                # Create lock file with exclusive access
                fd = os.open(str(self._lock_file), 
                            os.O_CREAT | os.O_EXCL | os.O_WRONLY)
                os.write(fd, f"{os.getpid()}\n".encode())
                os.close(fd)
                return True
            except FileExistsError:
                time.sleep(0.1)
        return False
    
    def _release_lock(self) -> None:
        """Release file lock"""
        try:
            self._lock_file.unlink()
        except:
            pass
    
    def get_state(self) -> StateSnapshot:
        """Get current state snapshot"""
        with self._lock:
            self._load_state()
            return deepcopy(self._local_state)
    
    def get_program_state(self, program: str) -> Optional[ProgramState]:
        """Get state for a specific program"""
        with self._lock:
            return self._local_state.programs.get(program)
    
    def update(self, key: str, value: Any, scope: str = 'program') -> bool:
        """
        Update state atomically
        
        Args:
            key: State key to update
            value: New value
            scope: 'program' (this program) or 'shared' (global)
        
        Returns:
            True if update succeeded
        """
        if not self._acquire_lock():
            logger.warning("Failed to acquire state lock")
            return False
        
        try:
            with self._lock:
                # Reload to get latest
                self._load_state()
                
                if scope == 'program':
                    # Update this program's state
                    if self.program_name not in self._local_state.programs:
                        self._local_state.programs[self.program_name] = ProgramState(
                            name=self.program_name
                        )
                    
                    prog_state = self._local_state.programs[self.program_name]
                    
                    if key == 'status':
                        prog_state.status = value
                    elif key == 'last_result':
                        prog_state.last_result = value
                    elif key == 'current_params':
                        prog_state.current_params = value
                    else:
                        # Generic update
                        setattr(prog_state, key, value)
                    
                    prog_state.version += 1
                    prog_state.last_heartbeat = datetime.now().isoformat()
                
                elif scope == 'shared':
                    # Update shared state
                    self._local_state.shared[key] = value
                
                # Increment global version
                self._local_state.version += 1
                self._local_state.timestamp = datetime.now().isoformat()
                
                self._save_state()
                return True
                
        finally:
            self._release_lock()
    
    def heartbeat(self) -> bool:
        """Send heartbeat to indicate program is alive"""
        return self.update('status', 'running')
    
    def set_idle(self) -> bool:
        """Mark program as idle"""
        return self.update('status', 'idle')
    
    def set_error(self, error: str = None) -> bool:
        """Mark program as having error"""
        self.update('last_error', error, scope='program')
        return self.update('status', 'error')
    
    def _sync_loop(self) -> None:
        """Background sync loop"""
        while self._running:
            try:
                old_version = self._local_state.version if self._local_state else 0
                self._load_state()
                
                if self._local_state and self._local_state.version != old_version:
                    # State changed externally
                    if self.on_state_changed:
                        self.on_state_changed(self._local_state)
                
                # Send heartbeat
                self.heartbeat()
                
            except Exception as e:
                logger.error(f"Sync error: {e}")
            
            time.sleep(self.sync_interval)
    
    def start(self) -> None:
        """Start background sync"""
        if self._running:
            return
        
        self._running = True
        self._sync_thread = threading.Thread(target=self._sync_loop, daemon=True)
        self._sync_thread.start()
        logger.info(f"StateSynchronizer started: {self.program_name}")
    
    def stop(self) -> None:
        """Stop background sync"""
        self._running = False
        if self._sync_thread:
            self._sync_thread.join(timeout=2.0)
        
        # Mark as not running
        self.update('status', 'stopped')
        logger.info("StateSynchronizer stopped")
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, *args):
        self.stop()


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE-LEVEL INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════
_synchronizer: Optional[StateSynchronizer] = None


def get_state_synchronizer(program_name: str = "Python") -> StateSynchronizer:
    """Get or create the module-level state synchronizer"""
    global _synchronizer
    if _synchronizer is None:
        _synchronizer = StateSynchronizer(program_name=program_name)
    return _synchronizer


def sync_state(program_name: str = "Python") -> StateSnapshot:
    """Convenience function to get current synchronized state"""
    return get_state_synchronizer(program_name).get_state()


# ═══════════════════════════════════════════════════════════════════════════════
# TEST
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("StateSynchronizer Test")
    print("=" * 40)
    
    sync = StateSynchronizer(
        state_dir=".test_state",
        program_name="TestPython"
    )
    
    # Update state
    sync.update('status', 'running')
    sync.update('last_result', {'F_U_Bi_i': 1.83e71})
    sync.update('active_query', 'NGC3596', scope='shared')
    
    # Get state
    state = sync.get_state()
    print(f"Version: {state.version}")
    print(f"Programs: {list(state.programs.keys())}")
    print(f"Shared: {state.shared}")
    
    print("\nTest complete")
