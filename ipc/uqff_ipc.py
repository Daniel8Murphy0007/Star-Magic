#!/usr/bin/env python3
"""
uqff_ipc.py - Python IPC Bridge for UQFF Simultaneous Joint Operation Pipeline
===============================================================================

Provides Python access to the shared memory and named pipe IPC channels
used by the 5 Principal Programs in the Star-Magic pipeline:

    1. source2.cpp (Principal GUI) - USER STARTS HERE
    2. QCalc.py (Pure Physics Solver) - Uses this module
    3. CondensedPhysics.py (Calculator) - Uses this module
    4. MAIN_1_CoAnQi.cpp (Physics Engine)
    5. source2(HEAD PROGRAM).cpp (VR/VM Simulator)

Data Flow:
    USER → source2.cpp → APIFetch.py → bodies_*.csv → IPData.py
                │
    ┌───────────┼───────────┬───────────┐
    ▼           ▼           ▼           ▼
  QCalc.py   CondensedPhys  MAIN_1     source2(HEAD
              ics.py        CoAnQi.cpp  PROGRAM).cpp
    │           │           │           │
    └───────────┴───────────┴───────────┘
                     │
                OPData.py → RECALL

Named Pipe: \\\\.\\pipe\\StarMagic_UQFF
Ports: 990 (FTPS), 21 (FTP+TLS), 3141 (HTTP), 8443 (HTTPS)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
Phase: 1 - IPC Pipeline Connection
"""

import os
import sys
import struct
import json
import mmap
import time
import logging
from typing import Dict, Any, Optional, Tuple
from enum import IntEnum
from dataclasses import dataclass, field, asdict
from datetime import datetime

# Setup logging
logger = logging.getLogger("UQFF_IPC")

# ═══════════════════════════════════════════════════════════════════════════════
# MESSAGE TYPES (Must match uqff_ipc.h MessageType enum)
# ═══════════════════════════════════════════════════════════════════════════════
class MessageType(IntEnum):
    """IPC Message types - mirror C++ enum from ipc/uqff_ipc.h"""
    PHYSICS_REQUEST = 0
    PHYSICS_RESPONSE = 1
    STATE_SYNC = 2
    CONTROL = 3
    HEARTBEAT = 4
    ERROR = 5
    JSON_DATA = 6

# ═══════════════════════════════════════════════════════════════════════════════
# MESSAGE HEADER (Must match uqff_ipc.h MessageHeader struct)
# ═══════════════════════════════════════════════════════════════════════════════
@dataclass
class MessageHeader:
    """IPC Message header - mirror C++ struct"""
    msg_type: MessageType = MessageType.JSON_DATA
    payload_size: int = 0
    sequence_id: int = 0
    sender_id: int = 0
    timestamp: float = 0.0
    
    HEADER_SIZE = 28  # 4 + 4 + 4 + 4 + 8 + padding = 28 bytes
    HEADER_FORMAT = "=IIIII d"  # unsigned ints + double
    
    def pack(self) -> bytes:
        """Pack header to bytes"""
        return struct.pack(
            self.HEADER_FORMAT,
            int(self.msg_type),
            self.payload_size,
            self.sequence_id,
            self.sender_id,
            0,  # reserved
            self.timestamp
        )
    
    @classmethod
    def unpack(cls, data: bytes) -> 'MessageHeader':
        """Unpack bytes to header"""
        if len(data) < cls.HEADER_SIZE:
            raise ValueError(f"Header data too short: {len(data)} < {cls.HEADER_SIZE}")
        values = struct.unpack(cls.HEADER_FORMAT, data[:cls.HEADER_SIZE])
        return cls(
            msg_type=MessageType(values[0]),
            payload_size=values[1],
            sequence_id=values[2],
            sender_id=values[3],
            timestamp=values[5]
        )

# ═══════════════════════════════════════════════════════════════════════════════
# IPC CHANNEL BASE
# ═══════════════════════════════════════════════════════════════════════════════
class IPCChannel:
    """Base IPC channel interface"""
    
    def __init__(self, name: str = "Python"):
        self.name = name
        self.connected = False
        self.sequence_counter = 0
        self.sender_id = hash(name) & 0xFFFFFFFF
    
    def connect(self) -> bool:
        """Connect to the IPC channel"""
        raise NotImplementedError
    
    def disconnect(self) -> None:
        """Disconnect from the IPC channel"""
        raise NotImplementedError
    
    def send(self, msg_type: MessageType, payload: bytes) -> bool:
        """Send a message"""
        raise NotImplementedError
    
    def receive(self, timeout_ms: int = 1000) -> Tuple[Optional[MessageHeader], Optional[bytes]]:
        """Receive a message"""
        raise NotImplementedError
    
    def send_json(self, msg_type: MessageType, data: Dict[str, Any]) -> bool:
        """Send JSON data"""
        payload = json.dumps(data).encode('utf-8')
        return self.send(msg_type, payload)
    
    def receive_json(self, timeout_ms: int = 1000) -> Tuple[Optional[MessageHeader], Optional[Dict[str, Any]]]:
        """Receive JSON data"""
        header, payload = self.receive(timeout_ms)
        if header and payload:
            try:
                data = json.loads(payload.decode('utf-8'))
                return header, data
            except json.JSONDecodeError as e:
                logger.error(f"JSON decode error: {e}")
                return header, None
        return None, None

# ═══════════════════════════════════════════════════════════════════════════════
# NAMED PIPE CHANNEL (Windows)
# ═══════════════════════════════════════════════════════════════════════════════
class NamedPipeChannel(IPCChannel):
    """Named pipe IPC channel for Windows"""
    
    PIPE_NAME = r"\\.\pipe\StarMagic_UQFF"
    
    def __init__(self, name: str = "Python"):
        super().__init__(name)
        self.pipe_handle = None
    
    def connect(self) -> bool:
        """Connect to the named pipe"""
        if sys.platform != 'win32':
            logger.warning("Named pipes only supported on Windows")
            return False
        
        try:
            import win32pipe
            import win32file
            import pywintypes
            
            # Try to connect to existing pipe (client mode)
            try:
                self.pipe_handle = win32file.CreateFile(
                    self.PIPE_NAME,
                    win32file.GENERIC_READ | win32file.GENERIC_WRITE,
                    0, None,
                    win32file.OPEN_EXISTING,
                    0, None
                )
                self.connected = True
                logger.info(f"IPC: {self.name} connected to named pipe")
                return True
            except pywintypes.error as e:
                if e.winerror == 2:  # Pipe does not exist
                    logger.debug(f"Named pipe not available (standalone mode)")
                else:
                    logger.warning(f"Pipe connection error: {e}")
                return False
                
        except ImportError:
            logger.debug("pywin32 not installed - named pipe unavailable")
            return False
        except Exception as e:
            logger.warning(f"Named pipe connect error: {e}")
            return False
    
    def disconnect(self) -> None:
        """Disconnect from the named pipe"""
        if self.pipe_handle:
            try:
                import win32file
                win32file.CloseHandle(self.pipe_handle)
            except Exception:
                pass
        self.pipe_handle = None
        self.connected = False
    
    def send(self, msg_type: MessageType, payload: bytes) -> bool:
        """Send a message through the pipe"""
        if not self.connected or not self.pipe_handle:
            return False
        
        try:
            import win32file
            
            header = MessageHeader(
                msg_type=msg_type,
                payload_size=len(payload),
                sequence_id=self.sequence_counter,
                sender_id=self.sender_id,
                timestamp=time.time()
            )
            self.sequence_counter += 1
            
            message = header.pack() + payload
            win32file.WriteFile(self.pipe_handle, message)
            return True
        except Exception as e:
            logger.error(f"Pipe send error: {e}")
            return False
    
    def receive(self, timeout_ms: int = 1000) -> Tuple[Optional[MessageHeader], Optional[bytes]]:
        """Receive a message from the pipe"""
        if not self.connected or not self.pipe_handle:
            return None, None
        
        try:
            import win32file
            
            # Read header first
            _, header_data = win32file.ReadFile(self.pipe_handle, MessageHeader.HEADER_SIZE)
            header = MessageHeader.unpack(header_data)
            
            # Read payload
            if header.payload_size > 0:
                _, payload = win32file.ReadFile(self.pipe_handle, header.payload_size)
            else:
                payload = b''
            
            return header, payload
        except Exception as e:
            logger.error(f"Pipe receive error: {e}")
            return None, None

# ═══════════════════════════════════════════════════════════════════════════════
# SHARED MEMORY CHANNEL
# ═══════════════════════════════════════════════════════════════════════════════
class SharedMemoryChannel(IPCChannel):
    """Shared memory IPC channel (cross-platform)"""
    
    SHM_NAME = "UQFF_SharedMemory"
    SHM_SIZE = 1024 * 1024  # 1 MB
    
    def __init__(self, name: str = "Python"):
        super().__init__(name)
        self.shm = None
        self.mm = None
    
    def connect(self) -> bool:
        """Connect to shared memory"""
        try:
            if sys.platform == 'win32':
                return self._connect_windows()
            else:
                return self._connect_posix()
        except Exception as e:
            logger.debug(f"Shared memory unavailable: {e}")
            return False
    
    def _connect_windows(self) -> bool:
        """Connect using Windows shared memory"""
        try:
            import win32api
            import win32file
            import win32security
            
            # Try to open existing shared memory
            try:
                self.shm = mmap.mmap(-1, self.SHM_SIZE, tagname=self.SHM_NAME, access=mmap.ACCESS_WRITE)
                self.connected = True
                logger.info(f"IPC: {self.name} connected to shared memory")
                return True
            except Exception:
                logger.debug("Shared memory not available (standalone mode)")
                return False
                
        except ImportError:
            # Fallback without win32api
            try:
                self.shm = mmap.mmap(-1, self.SHM_SIZE, tagname=self.SHM_NAME, access=mmap.ACCESS_WRITE)
                self.connected = True
                return True
            except Exception:
                return False
    
    def _connect_posix(self) -> bool:
        """Connect using POSIX shared memory"""
        try:
            # POSIX shared memory
            fd = os.open(f"/dev/shm/{self.SHM_NAME}", os.O_RDWR | os.O_CREAT, 0o666)
            os.ftruncate(fd, self.SHM_SIZE)
            self.mm = mmap.mmap(fd, self.SHM_SIZE)
            os.close(fd)
            self.connected = True
            return True
        except Exception:
            return False
    
    def disconnect(self) -> None:
        """Disconnect from shared memory"""
        if self.mm:
            self.mm.close()
        if self.shm:
            self.shm.close()
        self.mm = None
        self.shm = None
        self.connected = False
    
    def send(self, msg_type: MessageType, payload: bytes) -> bool:
        """Send a message through shared memory"""
        mem = self.mm or self.shm
        if not self.connected or not mem:
            return False
        
        try:
            header = MessageHeader(
                msg_type=msg_type,
                payload_size=len(payload),
                sequence_id=self.sequence_counter,
                sender_id=self.sender_id,
                timestamp=time.time()
            )
            self.sequence_counter += 1
            
            message = header.pack() + payload
            
            # Write to shared memory with lock
            mem.seek(0)
            mem.write(struct.pack("=I", len(message)))  # Length prefix
            mem.write(message)
            mem.flush()
            
            return True
        except Exception as e:
            logger.error(f"Shared memory send error: {e}")
            return False
    
    def receive(self, timeout_ms: int = 1000) -> Tuple[Optional[MessageHeader], Optional[bytes]]:
        """Receive a message from shared memory"""
        mem = self.mm or self.shm
        if not self.connected or not mem:
            return None, None
        
        try:
            mem.seek(0)
            length_data = mem.read(4)
            length = struct.unpack("=I", length_data)[0]
            
            if length == 0:
                return None, None
            
            header_data = mem.read(MessageHeader.HEADER_SIZE)
            header = MessageHeader.unpack(header_data)
            
            if header.payload_size > 0:
                payload = mem.read(header.payload_size)
            else:
                payload = b''
            
            return header, payload
        except Exception as e:
            logger.error(f"Shared memory receive error: {e}")
            return None, None

# ═══════════════════════════════════════════════════════════════════════════════
# UQFF IPC CLIENT (Main Interface)
# ═══════════════════════════════════════════════════════════════════════════════
class UQFFIPCClient:
    """
    Main IPC client for Python calculators.
    
    Usage:
        from ipc.uqff_ipc import UQFFIPCClient
        
        ipc = UQFFIPCClient("QCalc")
        if ipc.connect():
            ipc.send_physics_result({...})
    """
    
    def __init__(self, name: str = "Python"):
        """Initialize IPC client"""
        self.name = name
        self.connected = False
        self.channel: Optional[IPCChannel] = None
    
    def connect(self) -> bool:
        """Connect to IPC pipeline (tries named pipe first, then shared memory)"""
        # Try named pipe first (Windows)
        if sys.platform == 'win32':
            pipe = NamedPipeChannel(self.name)
            if pipe.connect():
                self.channel = pipe
                self.connected = True
                return True
        
        # Fall back to shared memory
        shm = SharedMemoryChannel(self.name)
        if shm.connect():
            self.channel = shm
            self.connected = True
            return True
        
        # Standalone mode (no IPC)
        logger.info(f"IPC: {self.name} running in standalone mode")
        self.connected = False
        return False
    
    def disconnect(self) -> None:
        """Disconnect from IPC pipeline"""
        if self.channel:
            self.channel.disconnect()
        self.channel = None
        self.connected = False
    
    def send_physics_result(self, result: Dict[str, Any]) -> bool:
        """Send physics calculation result to pipeline"""
        if not self.connected or not self.channel:
            return False
        
        return self.channel.send_json(MessageType.PHYSICS_RESPONSE, result)
    
    def request_physics(self, params: Dict[str, Any]) -> bool:
        """Request physics calculation from another program"""
        if not self.connected or not self.channel:
            return False
        
        return self.channel.send_json(MessageType.PHYSICS_REQUEST, params)
    
    def send_state_sync(self, state: Dict[str, Any]) -> bool:
        """Send state synchronization to pipeline"""
        if not self.connected or not self.channel:
            return False
        
        return self.channel.send_json(MessageType.STATE_SYNC, state)
    
    def receive(self, timeout_ms: int = 1000) -> Tuple[Optional[MessageType], Optional[Dict[str, Any]]]:
        """Receive message from pipeline"""
        if not self.connected or not self.channel:
            return None, None
        
        header, data = self.channel.receive_json(timeout_ms)
        if header and data:
            return header.msg_type, data
        return None, None
    
    def heartbeat(self) -> bool:
        """Send heartbeat to pipeline"""
        if not self.connected or not self.channel:
            return False
        
        return self.channel.send_json(MessageType.HEARTBEAT, {
            "name": self.name,
            "timestamp": datetime.now().isoformat(),
            "status": "alive"
        })

# ═══════════════════════════════════════════════════════════════════════════════
# MODULE-LEVEL IPC INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════
_ipc_client: Optional[UQFFIPCClient] = None

def get_ipc_client(name: str = "Python") -> UQFFIPCClient:
    """Get or create the module-level IPC client"""
    global _ipc_client
    if _ipc_client is None:
        _ipc_client = UQFFIPCClient(name)
        _ipc_client.connect()
    return _ipc_client

def ipc_connected() -> bool:
    """Check if IPC is connected"""
    global _ipc_client
    return _ipc_client is not None and _ipc_client.connected

# ═══════════════════════════════════════════════════════════════════════════════
# TEST / DEMO
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    
    print("UQFF IPC Bridge Test")
    print("=" * 40)
    
    client = UQFFIPCClient("TestPython")
    connected = client.connect()
    
    print(f"Connected: {connected}")
    print(f"Channel: {type(client.channel).__name__ if client.channel else 'None'}")
    
    if connected:
        # Send test message
        result = client.send_physics_result({
            "source": "Python Test",
            "F_U_Bi_i": 1.83e71,
            "timestamp": datetime.now().isoformat()
        })
        print(f"Send result: {result}")
        
        client.disconnect()
    
    print("Test complete")
