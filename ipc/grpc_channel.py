"""
grpc_channel.py - gRPC IPC Channel for UQFF Pipeline
=====================================================

Cross-platform IPC using gRPC (Protocol Buffers).
Works on Windows, Linux, and macOS.

Requires:
    pip install grpcio grpcio-tools

Proto compilation:
    python -m grpc_tools.protoc -I. --python_out=. --grpc_python_out=. uqff_service.proto

Usage:
    from ipc.grpc_channel import GrpcChannel, start_grpc_server
    
    # Client
    channel = GrpcChannel('localhost:50051')
    result = channel.calculate_field('Sagittarius A*', r=1e12)
    
    # Server
    server = start_grpc_server(port=50051)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

import os
import sys
import json
import logging
import threading
from typing import Dict, Any, Optional, List, Iterator
from dataclasses import dataclass
from datetime import datetime
from concurrent import futures

logger = logging.getLogger("UQFF_gRPC")

# Check if grpcio is available
try:
    import grpc
    GRPC_AVAILABLE = True
except ImportError:
    GRPC_AVAILABLE = False
    logger.warning("grpcio not installed - gRPC channel unavailable")

# Generated protobuf imports (compile with grpc_tools.protoc)
try:
    from . import uqff_service_pb2
    from . import uqff_service_pb2_grpc
    PROTO_AVAILABLE = True
except ImportError:
    PROTO_AVAILABLE = False
    # Define stub classes if proto not compiled
    uqff_service_pb2 = None
    uqff_service_pb2_grpc = None


@dataclass
class FieldResult:
    """Result from field calculation"""
    system_name: str
    F_U_Bi_i: float
    compressed_g: float
    resonance_a: float
    ug1_i: float
    ug2_i: float
    ug3_i: float
    ug4_i: float
    ubi_i: float
    validation: float
    equation: str
    timestamp: str


class GrpcChannel:
    """
    gRPC client channel for UQFF physics calculations
    
    Connects to a running GrpcServer and enables remote physics
    calculations across processes or machines.
    """
    
    DEFAULT_PORT = 50051
    
    def __init__(self, 
                 address: str = 'localhost:50051',
                 timeout: float = 30.0):
        """
        Initialize gRPC channel
        
        Args:
            address: Server address (host:port)
            timeout: RPC timeout in seconds
        """
        if not GRPC_AVAILABLE:
            raise ImportError("grpcio not installed. Run: pip install grpcio")
        
        self.address = address
        self.timeout = timeout
        self._channel: Optional[grpc.Channel] = None
        self._stub = None
        self._connected = False
    
    def connect(self) -> bool:
        """Connect to gRPC server"""
        try:
            self._channel = grpc.insecure_channel(self.address)
            
            # Wait for channel to be ready
            try:
                grpc.channel_ready_future(self._channel).result(timeout=5)
            except grpc.FutureTimeoutError:
                logger.warning(f"gRPC server not responding: {self.address}")
                return False
            
            if PROTO_AVAILABLE:
                self._stub = uqff_service_pb2_grpc.PhysicsServiceStub(self._channel)
            
            self._connected = True
            logger.info(f"gRPC connected: {self.address}")
            return True
            
        except Exception as e:
            logger.error(f"gRPC connect error: {e}")
            return False
    
    def disconnect(self) -> None:
        """Disconnect from gRPC server"""
        if self._channel:
            self._channel.close()
        self._channel = None
        self._stub = None
        self._connected = False
    
    def is_connected(self) -> bool:
        """Check if connected"""
        return self._connected and self._channel is not None
    
    def calculate_field(self, 
                        system_name: str,
                        r: float = 1e12,
                        t: float = 0.0,
                        theta: float = 0.0,
                        flags: int = 0) -> Optional[FieldResult]:
        """
        Calculate unified field for a system
        
        Args:
            system_name: System identifier (e.g., 'Sagittarius A*')
            r: Radial distance [m]
            t: Time [s]
            theta: Angle [rad]
            flags: Calculation options
        
        Returns:
            FieldResult or None on error
        """
        if not self._connected or not PROTO_AVAILABLE:
            return None
        
        try:
            request = uqff_service_pb2.FieldRequest(
                system_name=system_name,
                r=r,
                t=t,
                theta=theta,
                flags=flags
            )
            
            response = self._stub.CalculateField(request, timeout=self.timeout)
            
            return FieldResult(
                system_name=system_name,
                F_U_Bi_i=response.F_U_Bi_i,
                compressed_g=response.compressed_g,
                resonance_a=response.resonance_a,
                ug1_i=response.ug1_i,
                ug2_i=response.ug2_i,
                ug3_i=response.ug3_i,
                ug4_i=response.ug4_i,
                ubi_i=response.ubi_i,
                validation=response.validation,
                equation=response.equation,
                timestamp=datetime.now().isoformat()
            )
            
        except Exception as e:
            logger.error(f"Calculate field error: {e}")
            return None
    
    def calculate_batch(self, systems: List[Dict[str, Any]]) -> List[Optional[FieldResult]]:
        """
        Calculate fields for multiple systems
        
        Args:
            systems: List of dicts with 'name', 'r', 't', 'theta' keys
        
        Returns:
            List of FieldResult (or None) for each system
        """
        if not self._connected or not PROTO_AVAILABLE:
            return [None] * len(systems)
        
        try:
            request = uqff_service_pb2.BatchRequest()
            for s in systems:
                req = request.requests.add()
                req.system_name = s.get('name', '')
                req.r = s.get('r', 1e12)
                req.t = s.get('t', 0.0)
                req.theta = s.get('theta', 0.0)
            
            response = self._stub.CalculateBatch(request, timeout=self.timeout * len(systems))
            
            results = []
            for resp in response.responses:
                results.append(FieldResult(
                    system_name=resp.system_name,
                    F_U_Bi_i=resp.F_U_Bi_i,
                    compressed_g=resp.compressed_g,
                    resonance_a=resp.resonance_a,
                    ug1_i=resp.ug1_i,
                    ug2_i=resp.ug2_i,
                    ug3_i=resp.ug3_i,
                    ug4_i=resp.ug4_i,
                    ubi_i=resp.ubi_i,
                    validation=resp.validation,
                    equation=resp.equation,
                    timestamp=datetime.now().isoformat()
                ))
            return results
            
        except Exception as e:
            logger.error(f"Batch calculate error: {e}")
            return [None] * len(systems)
    
    def stream_updates(self, 
                       system_name: str,
                       update_rate_hz: float = 60.0) -> Iterator[FieldResult]:
        """
        Stream real-time field updates (for VR)
        
        Args:
            system_name: System to stream
            update_rate_hz: Update rate in Hz
        
        Yields:
            FieldResult for each update
        """
        if not self._connected or not PROTO_AVAILABLE:
            return
        
        try:
            request = uqff_service_pb2.StreamRequest(
                system_name=system_name,
                update_rate_hz=update_rate_hz
            )
            
            for response in self._stub.StreamFieldUpdates(request):
                yield FieldResult(
                    system_name=system_name,
                    F_U_Bi_i=response.F_U_Bi_i,
                    compressed_g=response.compressed_g,
                    resonance_a=response.resonance_a,
                    ug1_i=response.ug1_i,
                    ug2_i=response.ug2_i,
                    ug3_i=response.ug3_i,
                    ug4_i=response.ug4_i,
                    ubi_i=response.ubi_i,
                    validation=response.validation,
                    equation=response.equation,
                    timestamp=datetime.now().isoformat()
                )
                
        except Exception as e:
            logger.error(f"Stream error: {e}")
    
    def get_status(self) -> Dict[str, Any]:
        """Get server status"""
        if not self._connected or not PROTO_AVAILABLE:
            return {'connected': False}
        
        try:
            request = uqff_service_pb2.StatusRequest()
            response = self._stub.GetStatus(request, timeout=5)
            
            return {
                'connected': True,
                'version': response.version,
                'uptime_seconds': response.uptime_seconds,
                'calculations': response.total_calculations,
                'systems': response.system_count,
                'status': response.status
            }
            
        except Exception as e:
            return {'connected': False, 'error': str(e)}
    
    def __enter__(self):
        self.connect()
        return self
    
    def __exit__(self, *args):
        self.disconnect()


class PhysicsServiceImpl:
    """
    gRPC server implementation for PhysicsService
    
    Wraps the actual physics calculations from CondensedPhysics.py
    and MAIN_1_CoAnQi.cpp.
    """
    
    def __init__(self):
        self.start_time = datetime.now()
        self.total_calculations = 0
        self._lock = threading.Lock()
        
        # Try to import physics calculators
        try:
            from CondensedPhysics import solve
            self._solve = solve
        except ImportError:
            self._solve = None
            logger.warning("CondensedPhysics not available for gRPC server")
    
    def CalculateField(self, request, context):
        """Handle single field calculation"""
        if not PROTO_AVAILABLE:
            context.set_code(grpc.StatusCode.UNAVAILABLE)
            return None
        
        with self._lock:
            self.total_calculations += 1
        
        # Perform calculation
        if self._solve:
            try:
                result = self._solve({
                    'system': request.system_name,
                    'r': request.r,
                    't': request.t,
                    'theta': request.theta
                })
                
                response = uqff_service_pb2.FieldResponse(
                    system_name=request.system_name,
                    F_U_Bi_i=result.get('F_U_Bi_i', 0),
                    compressed_g=result.get('compressed_g', 0),
                    resonance_a=result.get('resonance_a', 0),
                    ug1_i=result.get('ug1_i', 0),
                    ug2_i=result.get('ug2_i', 0),
                    ug3_i=result.get('ug3_i', 0),
                    ug4_i=result.get('ug4_i', 0),
                    ubi_i=result.get('ubi_i', 0),
                    validation=result.get('validation', 0),
                    equation=result.get('equation', '')
                )
                return response
                
            except Exception as e:
                logger.error(f"Calculation error: {e}")
                context.set_code(grpc.StatusCode.INTERNAL)
                context.set_details(str(e))
                return uqff_service_pb2.FieldResponse()
        else:
            context.set_code(grpc.StatusCode.UNAVAILABLE)
            context.set_details("Physics calculator not available")
            return uqff_service_pb2.FieldResponse()
    
    def GetStatus(self, request, context):
        """Return server status"""
        if not PROTO_AVAILABLE:
            return None
        
        uptime = (datetime.now() - self.start_time).total_seconds()
        
        return uqff_service_pb2.StatusResponse(
            version="3.0.0",
            uptime_seconds=uptime,
            total_calculations=self.total_calculations,
            system_count=106,  # From index.js
            status="running"
        )


class GrpcServer:
    """
    gRPC server for UQFF physics service
    """
    
    def __init__(self, 
                 port: int = 50051,
                 max_workers: int = 10):
        """
        Initialize gRPC server
        
        Args:
            port: Port to listen on
            max_workers: Maximum thread pool workers
        """
        if not GRPC_AVAILABLE:
            raise ImportError("grpcio not installed")
        
        self.port = port
        self.max_workers = max_workers
        self._server: Optional[grpc.Server] = None
        self._running = False
    
    def start(self) -> bool:
        """Start the gRPC server"""
        if self._running:
            return True
        
        if not PROTO_AVAILABLE:
            logger.error("Proto files not compiled - cannot start server")
            return False
        
        try:
            self._server = grpc.server(
                futures.ThreadPoolExecutor(max_workers=self.max_workers)
            )
            
            service = PhysicsServiceImpl()
            uqff_service_pb2_grpc.add_PhysicsServiceServicer_to_server(
                service, self._server
            )
            
            self._server.add_insecure_port(f'[::]:{self.port}')
            self._server.start()
            
            self._running = True
            logger.info(f"gRPC server started on port {self.port}")
            return True
            
        except Exception as e:
            logger.error(f"gRPC server start error: {e}")
            return False
    
    def stop(self, grace: float = 5.0) -> None:
        """Stop the gRPC server"""
        if self._server:
            self._server.stop(grace)
        self._running = False
        logger.info("gRPC server stopped")
    
    def wait(self) -> None:
        """Wait for server termination"""
        if self._server:
            self._server.wait_for_termination()
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, *args):
        self.stop()


def start_grpc_server(port: int = 50051) -> Optional[GrpcServer]:
    """
    Convenience function to start gRPC server
    
    Args:
        port: Port to listen on
    
    Returns:
        Running GrpcServer or None on error
    """
    if not GRPC_AVAILABLE:
        logger.error("grpcio not installed")
        return None
    
    server = GrpcServer(port=port)
    if server.start():
        return server
    return None


def compile_proto() -> bool:
    """
    Compile proto file to Python modules
    
    Run this once to generate uqff_service_pb2.py and uqff_service_pb2_grpc.py
    """
    try:
        from grpc_tools import protoc
        
        proto_path = os.path.dirname(__file__)
        proto_file = os.path.join(proto_path, 'uqff_service.proto')
        
        result = protoc.main([
            'grpc_tools.protoc',
            f'-I{proto_path}',
            f'--python_out={proto_path}',
            f'--grpc_python_out={proto_path}',
            proto_file
        ])
        
        if result == 0:
            logger.info("Proto compilation successful")
            return True
        else:
            logger.error(f"Proto compilation failed: {result}")
            return False
            
    except ImportError:
        logger.error("grpcio-tools not installed. Run: pip install grpcio-tools")
        return False
    except Exception as e:
        logger.error(f"Proto compilation error: {e}")
        return False


# ═══════════════════════════════════════════════════════════════════════════════
# TEST
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("gRPC Channel Test")
    print("=" * 40)
    
    if not GRPC_AVAILABLE:
        print("grpcio not installed - run: pip install grpcio")
        sys.exit(1)
    
    if not PROTO_AVAILABLE:
        print("Proto not compiled - attempting compilation...")
        if not compile_proto():
            print("Failed to compile proto. Run manually:")
            print("  python -m grpc_tools.protoc -I. --python_out=. --grpc_python_out=. uqff_service.proto")
            sys.exit(1)
    
    # Test client connection
    channel = GrpcChannel('localhost:50051')
    if channel.connect():
        status = channel.get_status()
        print(f"Server status: {status}")
        channel.disconnect()
    else:
        print("Could not connect to gRPC server")
        print("Start server first with: python -c 'from ipc.grpc_channel import *; s = start_grpc_server(); s.wait()'")
    
    print("\nTest complete")
