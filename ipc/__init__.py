r"""
ipc - UQFF Inter-Process Communication Package
===============================================

Provides IPC channels for the 5 Principal Programs in Star-Magic:
1. source2.cpp (Principal GUI)
2. QCalc.py (Pure Physics Solver)
3. CondensedPhysics.py (Calculator)
4. MAIN_1_CoAnQi.cpp (Physics Engine)
5. source2(HEAD PROGRAM).cpp (VR/VM Simulator)

Channels:
---------
- NamedPipeChannel: Windows named pipes (\\.\pipe\StarMagic_UQFF)
- SharedMemoryChannel: Cross-platform shared memory (mmap)
- GrpcChannel: Cross-platform gRPC (protobuf)

Utilities:
----------
- FileWatcher: Watch for bodies_*.csv changes
- ResultCache: LRU cache for physics results
- StateSynchronizer: Multi-program state sync
- FileLock: Lock-free file access patterns

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

# Define __all__ early for wildcard imports
__all__ = [
    # Core IPC
    'MessageType',
    'MessageHeader',
    'IPCChannel',
    'NamedPipeChannel',
    'SharedMemoryChannel',
    'UQFFIPCClient',
    'get_ipc_client',
    'ipc_connected',
    
    # File Watcher
    'FileWatcher',
    'BodiesFileWatcher',
    'watch_bodies_directory',
    
    # Result Cache
    'ResultCache',
    'get_result_cache',
    'cache_physics_result',
    'get_cached_result',
    
    # State Sync
    'StateSynchronizer',
    'StateSnapshot',
    'get_state_synchronizer',
    'sync_state',
    
    # File Lock
    'FileLock',
    'atomic_write',
    'safe_read',
    
    # gRPC (conditional)
    'GRPC_AVAILABLE',
]

__version__ = '3.0.0'

# Core IPC classes
from .uqff_ipc import (
    MessageType,
    MessageHeader,
    IPCChannel,
    NamedPipeChannel,
    SharedMemoryChannel,
    UQFFIPCClient,
    get_ipc_client,
    ipc_connected,
)

# File system utilities
from .file_watcher import (
    FileWatcher,
    BodiesFileWatcher,
    watch_bodies_directory,
)

from .result_cache import (
    ResultCache,
    get_result_cache,
    cache_physics_result,
    get_cached_result,
)

from .state_sync import (
    StateSynchronizer,
    StateSnapshot,
    get_state_synchronizer,
    sync_state,
)

from .file_lock import (
    FileLock,
    atomic_write,
    safe_read,
)

# gRPC channel (optional - requires grpcio)
try:
    from .grpc_channel import (
        GrpcChannel,
        GrpcServer,
        start_grpc_server,
    )
    GRPC_AVAILABLE = True
    __all__.extend(['GrpcChannel', 'GrpcServer', 'start_grpc_server'])
except ImportError:
    GRPC_AVAILABLE = False

__version__ = '3.0.0'
