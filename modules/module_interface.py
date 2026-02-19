#!/usr/bin/env python3
"""
Module Interface - Abstract Base Class for Plug/Play Modules
=============================================================

All modules (AI Clones, Gaming, etc.) must inherit from ModuleInterface.

Architecture:
    ModuleInterface (ABC)
    ├── AICloneModule          # Encryption, Astrological, Quantum calculators
    ├── GamingModule           # Wormhole Traveler, Cosmic Mapping, Physics Paradigm
    └── DebugModule            # Health monitoring, hot-reload, diagnostics

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Callable
from enum import Enum, auto
import hashlib
import json
import uuid
from datetime import datetime


class ModuleType(Enum):
    """Module category enumeration."""
    AI_CLONE_ENCRYPTION = auto()
    AI_CLONE_ASTROLOGICAL = auto()
    AI_CLONE_QUANTUM = auto()
    GAMING_WORMHOLE = auto()
    GAMING_COSMIC_MAP = auto()
    GAMING_COSMIC_MAPPING = auto()  # Galaxy renderer
    GAMING_PHYSICS = auto()
    GAMING_PHYSICS_PARADIGM = auto()  # UQFF vs GR comparison
    DEBUG_HEALTH = auto()
    DEBUG_HOTRELOAD = auto()
    DEBUG_DIAGNOSTIC = auto()
    DEBUG_SANDBOX = auto()


class ModuleFormat(Enum):
    """Supported module file formats."""
    CPP = "cpp"
    JSON = "json"
    PYTHON = "py"
    QT = "qt"
    DATABASE = "db"
    JIT = "jit"


@dataclass
class ModuleCapabilities:
    """Declares what a module can do."""
    can_hot_reload: bool = False
    requires_sandbox: bool = False
    supports_jit: bool = False
    requires_gpu: bool = False
    max_memory_mb: int = 512
    max_cpu_percent: float = 25.0
    allowed_syscalls: List[str] = field(default_factory=list)
    dependencies: List[str] = field(default_factory=list)


@dataclass
class ModuleMetadata:
    """Module identification and version info."""
    uuid: str = field(default_factory=lambda: str(uuid.uuid4()))
    name: str = ""
    version: str = "1.0.0"
    author: str = "Daniel T. Murphy"
    description: str = ""
    module_type: ModuleType = ModuleType.AI_CLONE_ENCRYPTION
    format: ModuleFormat = ModuleFormat.PYTHON
    created_at: str = field(default_factory=lambda: datetime.now().isoformat())
    checksum: str = ""
    signature: str = ""


@dataclass
class ModuleState:
    """Runtime state of a module."""
    is_loaded: bool = False
    is_running: bool = False
    last_execution: Optional[str] = None
    execution_count: int = 0
    error_count: int = 0
    last_error: Optional[str] = None
    data: Dict[str, Any] = field(default_factory=dict)


class ModuleInterface(ABC):
    """
    Abstract base class for all plug/play modules.
    
    Every module must implement:
    - load(): Initialize the module
    - unload(): Cleanup resources
    - execute(): Run the module's main logic
    - verify(): Validate module integrity
    
    Optional:
    - on_hot_reload(): Handle live code updates
    - export_state(): Serialize state for persistence
    - import_state(): Restore state from persistence
    """
    
    def __init__(self):
        self.metadata = ModuleMetadata()
        self.capabilities = ModuleCapabilities()
        self.state = ModuleState()
        self._event_handlers: Dict[str, List[Callable]] = {}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # ABSTRACT METHODS (MUST IMPLEMENT)
    # ═══════════════════════════════════════════════════════════════════════════
    
    @abstractmethod
    def load(self) -> bool:
        """
        Initialize the module. Called once when module is first loaded.
        
        Returns:
            bool: True if successfully loaded, False otherwise.
        """
        pass
    
    @abstractmethod
    def unload(self) -> bool:
        """
        Cleanup resources. Called when module is being removed.
        
        Returns:
            bool: True if successfully unloaded, False otherwise.
        """
        pass
    
    @abstractmethod
    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run the module's main computation.
        
        Args:
            inputs: Dictionary of input parameters.
            
        Returns:
            Dictionary of output results.
        """
        pass
    
    @abstractmethod
    def verify(self) -> bool:
        """
        Validate module integrity (checksum, signature, dependencies).
        
        Returns:
            bool: True if module passes all integrity checks.
        """
        pass
    
    # ═══════════════════════════════════════════════════════════════════════════
    # OPTIONAL METHODS (OVERRIDE IF NEEDED)
    # ═══════════════════════════════════════════════════════════════════════════
    
    def on_hot_reload(self, new_code: bytes) -> bool:
        """
        Handle live code updates without restart.
        
        Args:
            new_code: The updated module code.
            
        Returns:
            bool: True if hot-reload succeeded.
        """
        if not self.capabilities.can_hot_reload:
            return False
        # Default: no-op, subclasses override
        return True
    
    def export_state(self) -> bytes:
        """
        Serialize module state for persistence.
        
        Returns:
            Serialized state as bytes.
        """
        return json.dumps({
            'metadata': {
                'uuid': self.metadata.uuid,
                'name': self.metadata.name,
                'version': self.metadata.version,
            },
            'state': {
                'execution_count': self.state.execution_count,
                'error_count': self.state.error_count,
                'data': self.state.data,
            }
        }).encode('utf-8')
    
    def import_state(self, data: bytes) -> bool:
        """
        Restore module state from persistence.
        
        Args:
            data: Serialized state bytes.
            
        Returns:
            bool: True if state was successfully restored.
        """
        try:
            state_dict = json.loads(data.decode('utf-8'))
            self.state.execution_count = state_dict.get('state', {}).get('execution_count', 0)
            self.state.error_count = state_dict.get('state', {}).get('error_count', 0)
            self.state.data = state_dict.get('state', {}).get('data', {})
            return True
        except Exception:
            return False
    
    # ═══════════════════════════════════════════════════════════════════════════
    # EVENT SYSTEM
    # ═══════════════════════════════════════════════════════════════════════════
    
    def on(self, event: str, handler: Callable) -> None:
        """Register an event handler."""
        if event not in self._event_handlers:
            self._event_handlers[event] = []
        self._event_handlers[event].append(handler)
    
    def emit(self, event: str, data: Any = None) -> None:
        """Emit an event to all registered handlers."""
        for handler in self._event_handlers.get(event, []):
            try:
                handler(data)
            except Exception as e:
                self.state.error_count += 1
                self.state.last_error = str(e)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UTILITY METHODS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def compute_checksum(self, data: bytes) -> str:
        """Compute SHA-256 checksum of data."""
        return hashlib.sha256(data).hexdigest()
    
    def get_info(self) -> Dict[str, Any]:
        """Get module information as a dictionary."""
        return {
            'uuid': self.metadata.uuid,
            'name': self.metadata.name,
            'version': self.metadata.version,
            'type': self.metadata.module_type.name,
            'format': self.metadata.format.value,
            'is_loaded': self.state.is_loaded,
            'is_running': self.state.is_running,
            'execution_count': self.state.execution_count,
            'error_count': self.state.error_count,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# SPECIALIZED BASE CLASSES
# ═══════════════════════════════════════════════════════════════════════════════

class AICloneModule(ModuleInterface):
    """Base class for AI Clone calculators (encryption, astrological, quantum)."""
    
    def __init__(self):
        super().__init__()
        self.capabilities.can_hot_reload = True
        self.capabilities.requires_sandbox = True
    
    @abstractmethod
    def calculate(self, params: Dict[str, Any]) -> Dict[str, Any]:
        """Perform the AI clone's primary calculation."""
        pass
    
    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Execute wraps calculate with state tracking."""
        self.state.is_running = True
        self.state.execution_count += 1
        self.state.last_execution = datetime.now().isoformat()
        
        try:
            result = self.calculate(inputs)
            result['success'] = True
            return result
        except Exception as e:
            self.state.error_count += 1
            self.state.last_error = str(e)
            return {'success': False, 'error': str(e)}
        finally:
            self.state.is_running = False


class GamingModule(ModuleInterface):
    """Base class for gaming modules (Wormhole Traveler, Cosmic Mapping, etc.)."""
    
    def __init__(self):
        super().__init__()
        self.capabilities.requires_gpu = True
        self.capabilities.max_memory_mb = 2048
        self.game_state: Dict[str, Any] = {}
    
    @abstractmethod
    def update(self, delta_time: float) -> None:
        """Update game state for one frame."""
        pass
    
    @abstractmethod
    def render(self) -> bytes:
        """Render current frame."""
        pass
    
    def execute(self, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """Execute game loop iteration."""
        delta_time = inputs.get('delta_time', 1.0/60.0)
        self.update(delta_time)
        frame_data = self.render()
        return {
            'frame': frame_data,
            'game_state': self.game_state,
            'success': True,
        }


class DebugModule(ModuleInterface):
    """Base class for debug/maintenance modules."""
    
    def __init__(self):
        super().__init__()
        self.capabilities.can_hot_reload = True
        self.capabilities.requires_sandbox = False  # Debug modules are trusted
    
    @abstractmethod
    def diagnose(self, target: Any) -> Dict[str, Any]:
        """Run diagnostics on target module/system."""
        pass


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE RESULT TYPES
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class CalculationResult:
    """Result from an AI Clone calculation."""
    success: bool
    value: Any
    equation_latex: str = ""
    computation_time_ms: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class GameFrame:
    """Single frame from a gaming module."""
    frame_number: int
    timestamp: float
    render_data: bytes
    game_state: Dict[str, Any]


@dataclass
class DiagnosticReport:
    """Report from a debug module."""
    target_uuid: str
    timestamp: str
    health_score: float  # 0.0 to 1.0
    issues: List[str]
    recommendations: List[str]
    metrics: Dict[str, Any]


if __name__ == "__main__":
    # Example: Show module type info
    print("Star-Magic Module Interface v3.0")
    print("=" * 50)
    print("\nModule Types:")
    for mt in ModuleType:
        print(f"  - {mt.name}")
    print("\nModule Formats:")
    for mf in ModuleFormat:
        print(f"  - {mf.value}")
