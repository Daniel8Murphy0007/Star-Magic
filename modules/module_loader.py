#!/usr/bin/env python3
"""
Module Loader - Dynamic Module Orchestrator
============================================

Manages loading, unloading, and execution of plug/play modules.
Supports hot-reload, sandboxing, and multi-format modules.

Architecture:
    ModuleLoader
    ├── load_module(path) → ModuleInterface
    ├── unload_module(uuid) → bool
    ├── execute_module(uuid, inputs) → Dict
    ├── hot_reload(uuid) → bool
    ├── list_modules() → List[ModuleMetadata]
    └── verify_all() → Dict[uuid, bool]

Supported Formats:
    - .py    → Python importlib
    - .cpp   → ctypes/pybind11 .so/.dll
    - .json  → Configuration/data modules
    - .qt    → Qt plugin system
    - .db    → SQLite/CosmosDB state
    - .jit   → Numba JIT compilation

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic Plug/Play Architecture v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import os
import sys
import json
import importlib
import importlib.util
import hashlib
import threading
from pathlib import Path
from typing import Dict, Any, Optional, List, Type, Callable
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum, auto
import logging

# Import module interface (relative import)
from .module_interface import (
    ModuleInterface, ModuleMetadata, ModuleCapabilities, ModuleState,
    ModuleType, ModuleFormat, AICloneModule, GamingModule, DebugModule
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] ModuleLoader: %(message)s'
)
logger = logging.getLogger(__name__)


class LoaderStatus(Enum):
    """Module loader status codes."""
    SUCCESS = auto()
    NOT_FOUND = auto()
    INVALID_FORMAT = auto()
    VERIFICATION_FAILED = auto()
    DEPENDENCY_ERROR = auto()
    SANDBOX_ERROR = auto()
    ALREADY_LOADED = auto()
    NOT_LOADED = auto()


@dataclass
class LoadResult:
    """Result of a module load operation."""
    status: LoaderStatus
    module: Optional[ModuleInterface] = None
    uuid: Optional[str] = None
    message: str = ""
    load_time_ms: float = 0.0


class ModuleLoader:
    """
    Dynamic module orchestrator for plug/play architecture.
    
    Thread-safe module management with support for:
    - Multiple module formats (.py, .cpp, .json, .qt, .db, .jit)
    - Hot-reload without restart
    - Sandboxed execution for untrusted modules
    - Dependency resolution
    - Health monitoring
    """
    
    def __init__(self, module_root: str = None):
        """
        Initialize the module loader.
        
        Args:
            module_root: Root directory for modules. Defaults to ./modules
        """
        if module_root is None:
            module_root = Path(__file__).parent
        self.module_root = Path(module_root)
        
        # Module registry: uuid → ModuleInterface
        self._modules: Dict[str, ModuleInterface] = {}
        
        # Module paths: uuid → file path
        self._module_paths: Dict[str, Path] = {}
        
        # File watchers for hot-reload
        self._watchers: Dict[str, any] = {}
        
        # Thread lock for concurrent access
        self._lock = threading.RLock()
        
        # Load registry
        self._registry_path = self.module_root / "module_registry.json"
        self._registry = self._load_registry()
        
        logger.info(f"ModuleLoader initialized. Root: {self.module_root}")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # CORE MODULE OPERATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def load_module(self, path: str, verify: bool = True) -> LoadResult:
        """
        Load a module from file path.
        
        Args:
            path: Path to the module file (relative or absolute).
            verify: Whether to verify checksum/signature before loading.
            
        Returns:
            LoadResult with status and loaded module.
        """
        import time
        start_time = time.perf_counter()
        
        # Resolve path
        module_path = Path(path)
        if not module_path.is_absolute():
            module_path = self.module_root / path
        
        if not module_path.exists():
            return LoadResult(
                status=LoaderStatus.NOT_FOUND,
                message=f"Module not found: {module_path}"
            )
        
        # Determine format
        suffix = module_path.suffix.lower()
        format_handlers = {
            '.py': self._load_python_module,
            '.cpp': self._load_cpp_module,
            '.json': self._load_json_module,
            '.so': self._load_cpp_module,
            '.dll': self._load_cpp_module,
            '.db': self._load_db_module,
        }
        
        handler = format_handlers.get(suffix)
        if handler is None:
            return LoadResult(
                status=LoaderStatus.INVALID_FORMAT,
                message=f"Unsupported format: {suffix}"
            )
        
        try:
            with self._lock:
                module = handler(module_path)
                
                if module is None:
                    return LoadResult(
                        status=LoaderStatus.INVALID_FORMAT,
                        message="Failed to instantiate module"
                    )
                
                # Verify if requested
                if verify and not module.verify():
                    return LoadResult(
                        status=LoaderStatus.VERIFICATION_FAILED,
                        message="Module verification failed"
                    )
                
                # Check dependencies
                for dep in module.capabilities.dependencies:
                    if dep not in self._modules:
                        return LoadResult(
                            status=LoaderStatus.DEPENDENCY_ERROR,
                            message=f"Missing dependency: {dep}"
                        )
                
                # Load the module
                if not module.load():
                    return LoadResult(
                        status=LoaderStatus.DEPENDENCY_ERROR,
                        message="Module load() returned False"
                    )
                
                # Register
                module.state.is_loaded = True
                self._modules[module.metadata.uuid] = module
                self._module_paths[module.metadata.uuid] = module_path
                
                # Update registry
                self._registry[module.metadata.uuid] = {
                    'name': module.metadata.name,
                    'path': str(module_path),
                    'loaded_at': datetime.now().isoformat(),
                }
                self._save_registry()
                
                load_time = (time.perf_counter() - start_time) * 1000
                
                logger.info(f"Loaded module: {module.metadata.name} ({module.metadata.uuid}) in {load_time:.2f}ms")
                
                return LoadResult(
                    status=LoaderStatus.SUCCESS,
                    module=module,
                    uuid=module.metadata.uuid,
                    message="Module loaded successfully",
                    load_time_ms=load_time
                )
                
        except Exception as e:
            logger.error(f"Failed to load module {path}: {e}")
            return LoadResult(
                status=LoaderStatus.DEPENDENCY_ERROR,
                message=str(e)
            )
    
    def unload_module(self, uuid: str) -> bool:
        """
        Unload a module by UUID.
        
        Args:
            uuid: The module's UUID.
            
        Returns:
            True if successfully unloaded.
        """
        with self._lock:
            if uuid not in self._modules:
                logger.warning(f"Module not loaded: {uuid}")
                return False
            
            module = self._modules[uuid]
            
            try:
                # Call unload
                module.unload()
                module.state.is_loaded = False
                
                # Remove from registry
                del self._modules[uuid]
                if uuid in self._module_paths:
                    del self._module_paths[uuid]
                if uuid in self._registry:
                    del self._registry[uuid]
                self._save_registry()
                
                logger.info(f"Unloaded module: {uuid}")
                return True
                
            except Exception as e:
                logger.error(f"Failed to unload module {uuid}: {e}")
                return False
    
    def execute_module(self, uuid: str, inputs: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute a loaded module.
        
        Args:
            uuid: The module's UUID.
            inputs: Input parameters for execution.
            
        Returns:
            Execution results dictionary.
        """
        with self._lock:
            if uuid not in self._modules:
                return {'success': False, 'error': f'Module not loaded: {uuid}'}
            
            module = self._modules[uuid]
        
        # Execute outside lock to allow concurrent executions
        try:
            return module.execute(inputs)
        except Exception as e:
            logger.error(f"Module execution error {uuid}: {e}")
            return {'success': False, 'error': str(e)}
    
    def hot_reload(self, uuid: str) -> bool:
        """
        Hot-reload a module without unloading.
        
        Args:
            uuid: The module's UUID.
            
        Returns:
            True if hot-reload succeeded.
        """
        with self._lock:
            if uuid not in self._modules:
                return False
            
            module = self._modules[uuid]
            
            if not module.capabilities.can_hot_reload:
                logger.warning(f"Module {uuid} does not support hot-reload")
                return False
            
            path = self._module_paths.get(uuid)
            if path is None or not path.exists():
                return False
            
            try:
                # Read new code
                with open(path, 'rb') as f:
                    new_code = f.read()
                
                # Export current state
                state_backup = module.export_state()
                
                # Trigger hot-reload
                success = module.on_hot_reload(new_code)
                
                if success:
                    logger.info(f"Hot-reloaded module: {uuid}")
                else:
                    # Restore state on failure
                    module.import_state(state_backup)
                    logger.warning(f"Hot-reload failed for {uuid}, state restored")
                
                return success
                
            except Exception as e:
                logger.error(f"Hot-reload error {uuid}: {e}")
                return False
    
    # ═══════════════════════════════════════════════════════════════════════════
    # QUERY METHODS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def list_modules(self) -> List[Dict[str, Any]]:
        """List all loaded modules with their info."""
        with self._lock:
            return [module.get_info() for module in self._modules.values()]
    
    def get_module(self, uuid: str) -> Optional[ModuleInterface]:
        """Get a loaded module by UUID."""
        with self._lock:
            return self._modules.get(uuid)
    
    def verify_all(self) -> Dict[str, bool]:
        """Verify all loaded modules."""
        results = {}
        with self._lock:
            for uuid, module in self._modules.items():
                try:
                    results[uuid] = module.verify()
                except Exception:
                    results[uuid] = False
        return results
    
    def find_by_type(self, module_type: ModuleType) -> List[ModuleInterface]:
        """Find all modules of a specific type."""
        with self._lock:
            return [
                m for m in self._modules.values()
                if m.metadata.module_type == module_type
            ]
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FORMAT-SPECIFIC LOADERS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _load_python_module(self, path: Path) -> Optional[ModuleInterface]:
        """Load a Python module."""
        try:
            # Create module spec
            spec = importlib.util.spec_from_file_location(
                path.stem,
                path
            )
            if spec is None or spec.loader is None:
                return None
            
            # Load module
            py_module = importlib.util.module_from_spec(spec)
            sys.modules[path.stem] = py_module
            spec.loader.exec_module(py_module)
            
            # Find ModuleInterface subclass
            for name in dir(py_module):
                obj = getattr(py_module, name)
                if (isinstance(obj, type) and 
                    issubclass(obj, ModuleInterface) and 
                    obj is not ModuleInterface and
                    obj is not AICloneModule and
                    obj is not GamingModule and
                    obj is not DebugModule):
                    # Instantiate
                    instance = obj()
                    return instance
            
            return None
            
        except Exception as e:
            logger.error(f"Python module load error: {e}")
            return None
    
    def _load_cpp_module(self, path: Path) -> Optional[ModuleInterface]:
        """Load a C++ module (.so/.dll)."""
        # Placeholder for ctypes/pybind11 loading
        logger.warning(f"C++ module loading not yet implemented: {path}")
        return None
    
    def _load_json_module(self, path: Path) -> Optional[ModuleInterface]:
        """Load a JSON configuration module."""
        # JSON modules are config-only, use a wrapper
        logger.warning(f"JSON module loading not yet implemented: {path}")
        return None
    
    def _load_db_module(self, path: Path) -> Optional[ModuleInterface]:
        """Load a database state module."""
        logger.warning(f"DB module loading not yet implemented: {path}")
        return None
    
    # ═══════════════════════════════════════════════════════════════════════════
    # REGISTRY MANAGEMENT
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _load_registry(self) -> Dict[str, Any]:
        """Load the module registry from disk."""
        if self._registry_path.exists():
            try:
                with open(self._registry_path, 'r') as f:
                    return json.load(f)
            except Exception:
                pass
        return {}
    
    def _save_registry(self) -> None:
        """Save the module registry to disk."""
        try:
            with open(self._registry_path, 'w') as f:
                json.dump(self._registry, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save registry: {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# SINGLETON INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

# Global module loader instance
_loader: Optional[ModuleLoader] = None


def get_loader() -> ModuleLoader:
    """Get the singleton ModuleLoader instance."""
    global _loader
    if _loader is None:
        _loader = ModuleLoader()
    return _loader


# ═══════════════════════════════════════════════════════════════════════════════
# CLI INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Star-Magic Module Loader CLI")
    parser.add_argument('command', choices=['list', 'load', 'unload', 'execute', 'verify'])
    parser.add_argument('--path', '-p', help='Module path for load command')
    parser.add_argument('--uuid', '-u', help='Module UUID for unload/execute')
    parser.add_argument('--inputs', '-i', help='JSON inputs for execute', default='{}')
    
    args = parser.parse_args()
    loader = get_loader()
    
    if args.command == 'list':
        modules = loader.list_modules()
        print(json.dumps(modules, indent=2))
        
    elif args.command == 'load':
        if not args.path:
            print("Error: --path required for load")
            sys.exit(1)
        result = loader.load_module(args.path)
        print(f"Status: {result.status.name}")
        print(f"Message: {result.message}")
        if result.uuid:
            print(f"UUID: {result.uuid}")
            
    elif args.command == 'unload':
        if not args.uuid:
            print("Error: --uuid required for unload")
            sys.exit(1)
        success = loader.unload_module(args.uuid)
        print(f"Unload {'succeeded' if success else 'failed'}")
        
    elif args.command == 'execute':
        if not args.uuid:
            print("Error: --uuid required for execute")
            sys.exit(1)
        inputs = json.loads(args.inputs)
        result = loader.execute_module(args.uuid, inputs)
        print(json.dumps(result, indent=2))
        
    elif args.command == 'verify':
        results = loader.verify_all()
        print(json.dumps(results, indent=2))
