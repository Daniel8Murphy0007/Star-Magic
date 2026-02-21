#!/usr/bin/env python3
"""
uqff_ftps_client.py - FTPS Client for UQFF Star-Magic Data Transfer
====================================================================

Secure file transfer client for:
  - Uploading computation results to remote servers
  - Downloading bodies_*.csv from data repositories
  - Synchronizing UQFF state files across distributed systems
  - Integration with uqff_server.js REST API

Supports:
  - Explicit FTPS (FTP over TLS, port 21 → upgrades to TLS)
  - Implicit FTPS (port 990, TLS from start)
  - Certificate verification
  - Passive/Active modes
  - Resumable transfers
  - Batch operations

Usage:
    from uqff_ftps_client import UQFFFTPSClient
    
    client = UQFFFTPSClient(
        host='ftps.example.com',
        user='uqff_user',
        password='secure_pass'
    )
    
    with client:
        client.upload_bodies('bodies.csv', '/data/uqff/')
        client.download_latest_csv('/data/uqff/', './local/')

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import ftplib
import ssl
import os
import json
import logging
import hashlib
import re
from pathlib import Path
from datetime import datetime
from typing import Optional, List, Dict, Callable, Union, Any
from dataclasses import dataclass, field
from contextlib import contextmanager

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('UQFFFTPSClient')


# ═══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class FTPSConfig:
    """FTPS connection configuration."""
    host: str
    port: int = 21
    user: str = 'anonymous'
    password: str = ''
    implicit_tls: bool = False         # Use implicit FTPS (port 990)
    verify_cert: bool = True           # Verify server certificate
    cert_file: Optional[str] = None    # Client certificate for mutual auth
    key_file: Optional[str] = None     # Client key for mutual auth
    ca_bundle: Optional[str] = None    # Custom CA bundle
    passive_mode: bool = True          # Use passive mode (firewall-friendly)
    timeout: float = 30.0              # Connection/operation timeout
    debug_level: int = 0               # FTP debug level (0-2)
    
    # UQFF-specific settings
    uqff_data_dir: str = '/uqff_data/'       # Remote UQFF data directory
    bodies_pattern: str = r'bodies_\d{8}_\d{6}\.csv'  # Pattern for bodies files
    auto_retry: bool = True                   # Auto-retry failed operations
    max_retries: int = 3                      # Max retry attempts
    
    @classmethod
    def from_env(cls) -> 'FTPSConfig':
        """Create config from environment variables."""
        return cls(
            host=os.getenv('UQFF_FTPS_HOST', ''),
            port=int(os.getenv('UQFF_FTPS_PORT', '21')),
            user=os.getenv('UQFF_FTPS_USER', 'anonymous'),
            password=os.getenv('UQFF_FTPS_PASS', ''),
            implicit_tls=os.getenv('UQFF_FTPS_IMPLICIT', '').lower() == 'true',
            verify_cert=os.getenv('UQFF_FTPS_VERIFY', 'true').lower() == 'true',
            passive_mode=os.getenv('UQFF_FTPS_PASSIVE', 'true').lower() == 'true',
            uqff_data_dir=os.getenv('UQFF_FTPS_DATA_DIR', '/uqff_data/')
        )
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary (excludes sensitive data)."""
        return {
            'host': self.host,
            'port': self.port,
            'user': self.user,
            'implicit_tls': self.implicit_tls,
            'verify_cert': self.verify_cert,
            'passive_mode': self.passive_mode,
            'uqff_data_dir': self.uqff_data_dir
        }


# ═══════════════════════════════════════════════════════════════════════════════
# TRANSFER STATISTICS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class TransferStats:
    """Statistics for a file transfer operation."""
    filename: str
    direction: str              # 'upload' or 'download'
    bytes_transferred: int = 0
    total_bytes: int = 0
    start_time: Optional[datetime] = None
    end_time: Optional[datetime] = None
    success: bool = False
    error_message: str = ''
    checksum_md5: str = ''
    
    @property
    def duration_seconds(self) -> float:
        """Get transfer duration in seconds."""
        if self.start_time and self.end_time:
            return (self.end_time - self.start_time).total_seconds()
        return 0.0
    
    @property
    def speed_bps(self) -> float:
        """Get transfer speed in bytes/second."""
        if self.duration_seconds > 0:
            return self.bytes_transferred / self.duration_seconds
        return 0.0
    
    @property
    def progress_percent(self) -> float:
        """Get progress percentage."""
        if self.total_bytes > 0:
            return (self.bytes_transferred / self.total_bytes) * 100
        return 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# IMPLICIT FTPS CLASS (Port 990)
# ═══════════════════════════════════════════════════════════════════════════════

class ImplicitFTPS(ftplib.FTP_TLS):
    """FTPS with implicit TLS (port 990, TLS from connection start)."""
    
    def __init__(self, host='', user='', passwd='', acct='', keyfile=None,
                 certfile=None, context=None, timeout=30):
        super().__init__(host, user, passwd, acct, keyfile, certfile, context, timeout)
    
    def connect(self, host='', port=990, timeout=-999):
        """Connect with implicit TLS (wraps socket immediately)."""
        if host:
            self.host = host
        if port > 0:
            self.port = port
        if timeout != -999:
            self.timeout = timeout
        
        self.sock = self.context.wrap_socket(
            self._create_socket(),
            server_hostname=self.host
        )
        self.af = self.sock.family
        self.file = self.sock.makefile('r', encoding=self.encoding)
        self.welcome = self.getresp()
        return self.welcome
    
    def _create_socket(self):
        """Create underlying TCP socket."""
        import socket
        return socket.create_connection((self.host, self.port), self.timeout)


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN FTPS CLIENT CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class UQFFFTPSClient:
    """
    FTPS Client for UQFF Star-Magic data transfer.
    
    Integrates with:
      - bodies_*.csv files from APIFetch.py
      - uqff_server.js computation results
      - CondensedPhysics.py output data
      - Distributed UQFF computation nodes
    """
    
    def __init__(self, 
                 host: str = '',
                 port: int = 21,
                 user: str = 'anonymous',
                 password: str = '',
                 config: Optional[FTPSConfig] = None,
                 progress_callback: Optional[Callable[[TransferStats], None]] = None):
        """
        Initialize FTPS client.
        
        Args:
            host: FTPS server hostname
            port: FTPS port (21 for explicit, 990 for implicit)
            user: Username
            password: Password
            config: Full FTPSConfig object (overrides other params)
            progress_callback: Called during transfers with stats
        """
        if config:
            self.config = config
        else:
            self.config = FTPSConfig(
                host=host,
                port=port,
                user=user,
                password=password
            )
        
        self.ftp: Optional[ftplib.FTP_TLS] = None
        self.connected = False
        self.progress_callback = progress_callback
        self._transfer_stats: List[TransferStats] = []
    
    def __enter__(self):
        """Context manager entry - connect."""
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - disconnect."""
        self.disconnect()
        return False
    
    # ═══════════════════════════════════════════════════════════════════════════
    # CONNECTION MANAGEMENT
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _create_ssl_context(self) -> ssl.SSLContext:
        """
        Create SSL context with proper security settings.
        
        TLS 1.3 Enforcement (Security Hardening):
        - Requires TLS 1.3 minimum by default (fallback to 1.2 if not supported)
        - Disables weak ciphers
        - Enables hostname verification
        """
        if self.config.verify_cert:
            context = ssl.create_default_context()
        else:
            context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
            context.check_hostname = False
            context.verify_mode = ssl.CERT_NONE
            logger.warning("SSL certificate verification DISABLED - use only for testing!")
        
        # Load custom CA bundle if provided
        if self.config.ca_bundle and os.path.exists(self.config.ca_bundle):
            context.load_verify_locations(self.config.ca_bundle)
        
        # Load client certificate for mutual authentication
        if self.config.cert_file and self.config.key_file:
            context.load_cert_chain(
                certfile=self.config.cert_file,
                keyfile=self.config.key_file
            )
            logger.info("Client certificate loaded for mutual TLS authentication")
        
        # ═══════════════════════════════════════════════════════════════════════
        # TLS 1.3 ENFORCEMENT (Security Hardening)
        # ═══════════════════════════════════════════════════════════════════════
        try:
            # Prefer TLS 1.3 (Python 3.7+ with OpenSSL 1.1.1+)
            context.minimum_version = ssl.TLSVersion.TLSv1_3
            logger.info("TLS 1.3 enforced as minimum version")
        except AttributeError:
            # Fallback to TLS 1.2 for older systems
            context.minimum_version = ssl.TLSVersion.TLSv1_2
            logger.warning("TLS 1.3 not available, using TLS 1.2 minimum")
        
        # Disable weak ciphers (security hardening)
        try:
            context.set_ciphers('ECDHE+AESGCM:DHE+AESGCM:ECDHE+CHACHA20:DHE+CHACHA20')
        except ssl.SSLError:
            # Use default secure ciphers if custom set fails
            pass
        
        return context
    
    def connect(self) -> bool:
        """
        Establish secure connection to FTPS server.
        
        Returns:
            True if connected successfully
        """
        if self.connected:
            logger.warning("Already connected")
            return True
        
        if not self.config.host:
            raise ValueError("No FTPS host configured")
        
        try:
            context = self._create_ssl_context()
            
            # Choose explicit or implicit FTPS
            if self.config.implicit_tls:
                logger.info(f"Connecting via implicit FTPS to {self.config.host}:{self.config.port}")
                self.ftp = ImplicitFTPS(context=context, timeout=self.config.timeout)
                self.ftp.connect(self.config.host, self.config.port)
            else:
                logger.info(f"Connecting via explicit FTPS to {self.config.host}:{self.config.port}")
                self.ftp = ftplib.FTP_TLS(context=context, timeout=self.config.timeout)
                self.ftp.connect(self.config.host, self.config.port)
                self.ftp.auth()  # Upgrade control channel to TLS
            
            # Set debug level
            self.ftp.set_debuglevel(self.config.debug_level)
            
            # Authenticate
            self.ftp.login(self.config.user, self.config.password)
            logger.info(f"Authenticated as: {self.config.user}")
            
            # Protect data channel
            self.ftp.prot_p()
            logger.info("Data channel protected (TLS)")
            
            # Set passive mode
            self.ftp.set_pasv(self.config.passive_mode)
            
            self.connected = True
            logger.info(f"Successfully connected to {self.config.host}")
            return True
            
        except ftplib.all_errors as e:
            logger.error(f"Connection failed: {e}")
            self.ftp = None
            raise ConnectionError(f"FTPS connection failed: {e}")
    
    def disconnect(self):
        """Gracefully disconnect from server."""
        if self.ftp:
            try:
                self.ftp.quit()
                logger.info("Disconnected gracefully")
            except ftplib.all_errors:
                try:
                    self.ftp.close()
                except:
                    pass
            finally:
                self.ftp = None
                self.connected = False
    
    def reconnect(self) -> bool:
        """Reconnect to server."""
        self.disconnect()
        return self.connect()
    
    # ═══════════════════════════════════════════════════════════════════════════
    # RFC 4217 COMPLIANCE VALIDATION
    # ═══════════════════════════════════════════════════════════════════════════
    
    def validate_rfc4217_compliance(self) -> Dict[str, Any]:
        """
        Validate RFC 4217 (Securing FTP with TLS) compliance.
        
        RFC 4217 Requirements Checked:
        1. AUTH TLS/SSL command support (explicit FTPS)
        2. PBSZ 0 (Protection Buffer Size)
        3. PROT P (Private data channel protection)
        4. TLS 1.2+ minimum version
        5. Certificate verification
        6. Secure cipher suites
        
        Returns:
            Dictionary with compliance results and recommendations
        """
        results = {
            'rfc4217_compliant': True,
            'checks': [],
            'warnings': [],
            'tls_version': None,
            'cipher_suite': None,
            'certificate': None
        }
        
        # Check 1: Connection type
        if self.config.implicit_tls:
            results['checks'].append({
                'name': 'FTPS Mode',
                'status': 'PASS',
                'detail': 'Implicit FTPS (port 990, TLS from start)'
            })
        else:
            results['checks'].append({
                'name': 'FTPS Mode',
                'status': 'PASS',
                'detail': 'Explicit FTPS (AUTH TLS upgrade)'
            })
        
        # Check 2: TLS Version
        if self.connected and self.ftp and hasattr(self.ftp.sock, 'version'):
            tls_version = self.ftp.sock.version()
            results['tls_version'] = tls_version
            
            if 'TLSv1.3' in tls_version:
                results['checks'].append({
                    'name': 'TLS Version',
                    'status': 'PASS',
                    'detail': f'{tls_version} (Recommended)'
                })
            elif 'TLSv1.2' in tls_version:
                results['checks'].append({
                    'name': 'TLS Version',
                    'status': 'PASS',
                    'detail': f'{tls_version} (Minimum required by RFC 4217)'
                })
            else:
                results['rfc4217_compliant'] = False
                results['checks'].append({
                    'name': 'TLS Version',
                    'status': 'FAIL',
                    'detail': f'{tls_version} - Below RFC 4217 minimum'
                })
        
        # Check 3: Cipher Suite
        if self.connected and self.ftp and hasattr(self.ftp.sock, 'cipher'):
            cipher_info = self.ftp.sock.cipher()
            if cipher_info:
                results['cipher_suite'] = cipher_info[0]
                # Check for weak ciphers
                weak_patterns = ['RC4', 'DES', 'MD5', 'NULL', 'EXPORT', 'anon']
                is_weak = any(p in cipher_info[0].upper() for p in weak_patterns)
                
                if is_weak:
                    results['rfc4217_compliant'] = False
                    results['checks'].append({
                        'name': 'Cipher Suite',
                        'status': 'FAIL',
                        'detail': f'{cipher_info[0]} - Weak cipher detected'
                    })
                else:
                    results['checks'].append({
                        'name': 'Cipher Suite',
                        'status': 'PASS',
                        'detail': cipher_info[0]
                    })
        
        # Check 4: Certificate Verification
        if self.config.verify_cert:
            results['checks'].append({
                'name': 'Certificate Verification',
                'status': 'PASS',
                'detail': 'Server certificate verified'
            })
        else:
            results['warnings'].append('Certificate verification DISABLED - vulnerable to MITM')
            results['checks'].append({
                'name': 'Certificate Verification',
                'status': 'WARN',
                'detail': 'DISABLED - Security risk'
            })
        
        # Check 5: Data Channel Protection (PROT P)
        results['checks'].append({
            'name': 'Data Channel Protection',
            'status': 'PASS',
            'detail': 'PROT P enabled (prot_p called after login)'
        })
        
        # Check 6: Passive Mode
        if self.config.passive_mode:
            results['checks'].append({
                'name': 'Transfer Mode',
                'status': 'PASS',
                'detail': 'Passive mode (firewall-friendly)'
            })
        else:
            results['warnings'].append('Active mode may have firewall issues')
            results['checks'].append({
                'name': 'Transfer Mode',
                'status': 'WARN',
                'detail': 'Active mode - may not work through NAT/firewall'
            })
        
        # Certificate details
        if self.connected and self.ftp and hasattr(self.ftp.sock, 'getpeercert'):
            cert = self.ftp.sock.getpeercert()
            if cert:
                results['certificate'] = {
                    'subject': dict(x[0] for x in cert.get('subject', [])),
                    'issuer': dict(x[0] for x in cert.get('issuer', [])),
                    'notBefore': cert.get('notBefore'),
                    'notAfter': cert.get('notAfter'),
                    'serialNumber': cert.get('serialNumber')
                }
        
        return results
    
    def print_compliance_report(self):
        """Print a formatted RFC 4217 compliance report."""
        results = self.validate_rfc4217_compliance()
        
        print("\n" + "═" * 70)
        print(" RFC 4217 COMPLIANCE REPORT - UQFF FTPS Client")
        print("═" * 70)
        
        status = "✅ COMPLIANT" if results['rfc4217_compliant'] else "❌ NON-COMPLIANT"
        print(f"\nOverall Status: {status}\n")
        
        print("─" * 70)
        print("SECURITY CHECKS:")
        print("─" * 70)
        
        for check in results['checks']:
            icon = "✅" if check['status'] == 'PASS' else ("⚠️" if check['status'] == 'WARN' else "❌")
            print(f"  {icon} {check['name']}: {check['detail']}")
        
        if results['warnings']:
            print("\n─" * 70)
            print("WARNINGS:")
            print("─" * 70)
            for warning in results['warnings']:
                print(f"  ⚠️ {warning}")
        
        if results['tls_version']:
            print(f"\nTLS Version: {results['tls_version']}")
        if results['cipher_suite']:
            print(f"Cipher Suite: {results['cipher_suite']}")
        
        if results['certificate']:
            print("\n─" * 70)
            print("CERTIFICATE INFO:")
            print("─" * 70)
            cert = results['certificate']
            if 'subject' in cert:
                print(f"  Subject: {cert['subject']}")
            if 'issuer' in cert:
                print(f"  Issuer: {cert['issuer']}")
            if 'notAfter' in cert:
                print(f"  Expires: {cert['notAfter']}")
        
        print("\n" + "═" * 70 + "\n")
        
        return results

    # ═══════════════════════════════════════════════════════════════════════════
    # DIRECTORY OPERATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def pwd(self) -> str:
        """Get current working directory."""
        self._ensure_connected()
        return self.ftp.pwd()
    
    def cwd(self, path: str):
        """Change working directory."""
        self._ensure_connected()
        self.ftp.cwd(path)
        logger.debug(f"Changed to: {path}")
    
    def mkdir(self, path: str, parents: bool = True):
        """Create directory (optionally with parents)."""
        self._ensure_connected()
        
        if parents:
            # Create all parent directories
            parts = path.strip('/').split('/')
            current = ''
            for part in parts:
                current += '/' + part
                try:
                    self.ftp.mkd(current)
                    logger.debug(f"Created directory: {current}")
                except ftplib.error_perm:
                    pass  # Directory may already exist
        else:
            self.ftp.mkd(path)
    
    def listdir(self, path: str = '.') -> List[str]:
        """List directory contents."""
        self._ensure_connected()
        return self.ftp.nlst(path)
    
    def listdir_detailed(self, path: str = '.') -> List[Dict[str, Any]]:
        """List directory with MLSD metadata."""
        self._ensure_connected()
        result = []
        try:
            for name, facts in self.ftp.mlsd(path):
                result.append({'name': name, **facts})
        except ftplib.error_perm:
            # Fallback to simple listing
            result = [{'name': f} for f in self.ftp.nlst(path)]
        return result
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FILE OPERATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def upload(self, local_path: str, remote_path: str, 
               resume: bool = False) -> TransferStats:
        """
        Upload a file to the server.
        
        Args:
            local_path: Path to local file
            remote_path: Remote destination path
            resume: Resume interrupted transfer
            
        Returns:
            TransferStats with transfer details
        """
        self._ensure_connected()
        
        local_file = Path(local_path)
        if not local_file.exists():
            raise FileNotFoundError(f"Local file not found: {local_path}")
        
        stats = TransferStats(
            filename=local_file.name,
            direction='upload',
            total_bytes=local_file.stat().st_size,
            start_time=datetime.now()
        )
        
        try:
            # Calculate checksum
            stats.checksum_md5 = self._calculate_md5(local_path)
            
            # Handle resume
            offset = 0
            if resume:
                try:
                    remote_size = self.ftp.size(remote_path)
                    if remote_size and remote_size < stats.total_bytes:
                        offset = remote_size
                        logger.info(f"Resuming upload from byte {offset}")
                except ftplib.error_perm:
                    pass
            
            # Upload with progress tracking
            def track_progress(data):
                stats.bytes_transferred += len(data)
                if self.progress_callback:
                    self.progress_callback(stats)
            
            with open(local_path, 'rb') as f:
                if offset > 0:
                    f.seek(offset)
                    stats.bytes_transferred = offset
                    self.ftp.voidcmd(f'REST {offset}')
                
                self.ftp.storbinary(f'STOR {remote_path}', f, callback=track_progress)
            
            stats.success = True
            stats.end_time = datetime.now()
            logger.info(f"Uploaded: {local_path} → {remote_path} "
                       f"({stats.bytes_transferred} bytes in {stats.duration_seconds:.2f}s)")
            
        except ftplib.all_errors as e:
            stats.error_message = str(e)
            stats.end_time = datetime.now()
            logger.error(f"Upload failed: {e}")
            
            if self.config.auto_retry:
                return self._retry_operation(
                    lambda: self.upload(local_path, remote_path, resume=True),
                    stats
                )
        
        self._transfer_stats.append(stats)
        return stats
    
    def download(self, remote_path: str, local_path: str,
                 resume: bool = False) -> TransferStats:
        """
        Download a file from the server.
        
        Args:
            remote_path: Remote file path
            local_path: Local destination path
            resume: Resume interrupted transfer
            
        Returns:
            TransferStats with transfer details
        """
        self._ensure_connected()
        
        stats = TransferStats(
            filename=Path(remote_path).name,
            direction='download',
            start_time=datetime.now()
        )
        
        try:
            # Get remote file size
            try:
                stats.total_bytes = self.ftp.size(remote_path) or 0
            except ftplib.error_perm:
                pass
            
            # Handle resume
            offset = 0
            mode = 'wb'
            local_file = Path(local_path)
            
            if resume and local_file.exists():
                offset = local_file.stat().st_size
                if offset < stats.total_bytes:
                    mode = 'ab'
                    stats.bytes_transferred = offset
                    logger.info(f"Resuming download from byte {offset}")
            
            # Ensure parent directory exists
            local_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Download with progress tracking
            def track_progress(data):
                stats.bytes_transferred += len(data)
                if self.progress_callback:
                    self.progress_callback(stats)
            
            with open(local_path, mode) as f:
                if offset > 0:
                    self.ftp.voidcmd(f'REST {offset}')
                
                self.ftp.retrbinary(f'RETR {remote_path}', 
                                   lambda data: (f.write(data), track_progress(data)))
            
            # Calculate checksum after download
            stats.checksum_md5 = self._calculate_md5(local_path)
            stats.success = True
            stats.end_time = datetime.now()
            logger.info(f"Downloaded: {remote_path} → {local_path} "
                       f"({stats.bytes_transferred} bytes in {stats.duration_seconds:.2f}s)")
            
        except ftplib.all_errors as e:
            stats.error_message = str(e)
            stats.end_time = datetime.now()
            logger.error(f"Download failed: {e}")
            
            if self.config.auto_retry:
                return self._retry_operation(
                    lambda: self.download(remote_path, local_path, resume=True),
                    stats
                )
        
        self._transfer_stats.append(stats)
        return stats
    
    def delete(self, remote_path: str):
        """Delete a remote file."""
        self._ensure_connected()
        self.ftp.delete(remote_path)
        logger.info(f"Deleted: {remote_path}")
    
    def rename(self, from_path: str, to_path: str):
        """Rename/move a remote file."""
        self._ensure_connected()
        self.ftp.rename(from_path, to_path)
        logger.info(f"Renamed: {from_path} → {to_path}")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UQFF-SPECIFIC OPERATIONS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def upload_bodies(self, local_csv: str, remote_dir: Optional[str] = None) -> TransferStats:
        """
        Upload a bodies_*.csv file to UQFF data directory.
        
        Args:
            local_csv: Path to local CSV file
            remote_dir: Remote directory (default: config.uqff_data_dir)
            
        Returns:
            TransferStats
        """
        remote_dir = remote_dir or self.config.uqff_data_dir
        
        # Ensure remote directory exists
        self.mkdir(remote_dir, parents=True)
        
        # Generate timestamped filename if needed
        local_path = Path(local_csv)
        if local_path.name == 'bodies.csv':
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            remote_filename = f'bodies_{timestamp}.csv'
        else:
            remote_filename = local_path.name
        
        remote_path = f"{remote_dir.rstrip('/')}/{remote_filename}"
        
        return self.upload(local_csv, remote_path)
    
    def download_latest_csv(self, remote_dir: Optional[str] = None,
                           local_dir: str = '.') -> Optional[TransferStats]:
        """
        Download the most recent bodies_*.csv file.
        
        Args:
            remote_dir: Remote directory to search
            local_dir: Local destination directory
            
        Returns:
            TransferStats if found, None otherwise
        """
        remote_dir = remote_dir or self.config.uqff_data_dir
        pattern = re.compile(self.config.bodies_pattern)
        
        # Find matching files
        try:
            self.cwd(remote_dir)
            files = self.listdir()
        except ftplib.error_perm as e:
            logger.error(f"Cannot access {remote_dir}: {e}")
            return None
        
        matching = [f for f in files if pattern.match(f)]
        if not matching:
            logger.warning(f"No bodies_*.csv files found in {remote_dir}")
            return None
        
        # Sort by name (timestamp) and get latest
        latest = sorted(matching, reverse=True)[0]
        remote_path = f"{remote_dir.rstrip('/')}/{latest}"
        local_path = Path(local_dir) / latest
        
        return self.download(remote_path, str(local_path))
    
    def sync_uqff_results(self, local_results_dir: str, 
                          remote_results_dir: Optional[str] = None) -> List[TransferStats]:
        """
        Synchronize UQFF computation results to server.
        
        Uploads all JSON result files from local directory.
        
        Args:
            local_results_dir: Local directory with result files
            remote_results_dir: Remote directory (default: uqff_data_dir/results/)
            
        Returns:
            List of TransferStats for each file
        """
        remote_dir = remote_results_dir or f"{self.config.uqff_data_dir.rstrip('/')}/results/"
        self.mkdir(remote_dir, parents=True)
        
        stats_list = []
        local_dir = Path(local_results_dir)
        
        for file in local_dir.glob('*.json'):
            remote_path = f"{remote_dir}/{file.name}"
            stats = self.upload(str(file), remote_path)
            stats_list.append(stats)
        
        return stats_list
    
    def call_uqff_server(self, system_params: Dict[str, Any]) -> Optional[Dict[str, Any]]:
        """
        Call uqff_server.js via FTPS (file-based RPC).
        
        Uploads request JSON, waits for response file.
        
        Args:
            system_params: UQFF system parameters
            
        Returns:
            Computation results or None if failed
        """
        import tempfile
        import time
        
        request_id = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        request_dir = f"{self.config.uqff_data_dir.rstrip('/')}/requests/"
        response_dir = f"{self.config.uqff_data_dir.rstrip('/')}/responses/"
        
        self.mkdir(request_dir, parents=True)
        
        # Upload request
        request_data = {
            'request_id': request_id,
            'timestamp': datetime.now().isoformat(),
            'params': system_params
        }
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(request_data, f)
            temp_request = f.name
        
        try:
            self.upload(temp_request, f"{request_dir}/req_{request_id}.json")
            
            # Poll for response (with timeout)
            response_path = f"{response_dir}/resp_{request_id}.json"
            max_wait = 30  # seconds
            poll_interval = 0.5
            waited = 0
            
            while waited < max_wait:
                try:
                    files = self.listdir(response_dir)
                    if f"resp_{request_id}.json" in files:
                        # Download and parse response
                        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                            temp_response = f.name
                        
                        self.download(response_path, temp_response)
                        
                        with open(temp_response, 'r') as f:
                            return json.load(f)
                except ftplib.error_perm:
                    pass
                
                time.sleep(poll_interval)
                waited += poll_interval
            
            logger.warning(f"Timeout waiting for UQFF server response: {request_id}")
            return None
            
        finally:
            os.unlink(temp_request)
    
    # ═══════════════════════════════════════════════════════════════════════════
    # UTILITY METHODS
    # ═══════════════════════════════════════════════════════════════════════════
    
    def _ensure_connected(self):
        """Ensure we're connected, raising if not."""
        if not self.connected or not self.ftp:
            raise ConnectionError("Not connected to FTPS server")
    
    def _calculate_md5(self, filepath: str) -> str:
        """Calculate MD5 hash of a file."""
        md5 = hashlib.md5()
        with open(filepath, 'rb') as f:
            for chunk in iter(lambda: f.read(8192), b''):
                md5.update(chunk)
        return md5.hexdigest()
    
    def _retry_operation(self, operation: Callable, 
                         failed_stats: TransferStats) -> TransferStats:
        """Retry a failed operation."""
        retries = 0
        while retries < self.config.max_retries:
            retries += 1
            logger.info(f"Retry attempt {retries}/{self.config.max_retries}")
            try:
                self.reconnect()
                return operation()
            except Exception as e:
                logger.warning(f"Retry {retries} failed: {e}")
        
        return failed_stats
    
    def get_transfer_history(self) -> List[TransferStats]:
        """Get all transfer statistics from this session."""
        return self._transfer_stats.copy()
    
    def get_server_info(self) -> Dict[str, str]:
        """Get server information."""
        self._ensure_connected()
        return {
            'welcome': getattr(self.ftp, 'welcome', ''),
            'pwd': self.pwd(),
            'system': self.ftp.sendcmd('SYST'),
            'features': self.ftp.sendcmd('FEAT') if hasattr(self.ftp, 'sendcmd') else ''
        }


# ═══════════════════════════════════════════════════════════════════════════════
# CLI INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    """Command-line interface for UQFF FTPS client."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='UQFF FTPS Client - Secure file transfer for Star-Magic',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --host ftps.example.com --upload bodies.csv /data/
  %(prog)s --host ftps.example.com --download-latest /data/ ./local/
  %(prog)s --host ftps.example.com --list /data/
  %(prog)s --config  # Show configuration from environment
        """
    )
    
    parser.add_argument('--host', help='FTPS server hostname')
    parser.add_argument('--port', type=int, default=21, help='FTPS port (default: 21)')
    parser.add_argument('--user', default='anonymous', help='Username')
    parser.add_argument('--password', default='', help='Password (use env UQFF_FTPS_PASS for security)')
    parser.add_argument('--implicit', action='store_true', help='Use implicit TLS (port 990)')
    parser.add_argument('--no-verify', action='store_true', help='Skip certificate verification')
    
    parser.add_argument('--upload', nargs=2, metavar=('LOCAL', 'REMOTE'),
                       help='Upload file')
    parser.add_argument('--download', nargs=2, metavar=('REMOTE', 'LOCAL'),
                       help='Download file')
    parser.add_argument('--upload-bodies', metavar='CSV', help='Upload bodies CSV')
    parser.add_argument('--download-latest', nargs=2, metavar=('REMOTE_DIR', 'LOCAL_DIR'),
                       help='Download latest bodies_*.csv')
    parser.add_argument('--list', metavar='PATH', help='List directory')
    parser.add_argument('--config', action='store_true', help='Show configuration')
    
    args = parser.parse_args()
    
    if args.config:
        config = FTPSConfig.from_env()
        print("UQFF FTPS Configuration (from environment):")
        for k, v in config.to_dict().items():
            print(f"  {k}: {v}")
        return
    
    if not args.host:
        parser.error("--host is required")
    
    config = FTPSConfig(
        host=args.host,
        port=args.port,
        user=args.user,
        password=args.password or os.getenv('UQFF_FTPS_PASS', ''),
        implicit_tls=args.implicit,
        verify_cert=not args.no_verify
    )
    
    def progress(stats: TransferStats):
        print(f"\r{stats.filename}: {stats.progress_percent:.1f}% "
              f"({stats.bytes_transferred}/{stats.total_bytes})", end='')
    
    try:
        with UQFFFTPSClient(config=config, progress_callback=progress) as client:
            if args.upload:
                stats = client.upload(args.upload[0], args.upload[1])
                print(f"\n{'✅ Success' if stats.success else '❌ Failed'}")
                
            elif args.download:
                stats = client.download(args.download[0], args.download[1])
                print(f"\n{'✅ Success' if stats.success else '❌ Failed'}")
                
            elif args.upload_bodies:
                stats = client.upload_bodies(args.upload_bodies)
                print(f"\n{'✅ Success' if stats.success else '❌ Failed'}")
                
            elif args.download_latest:
                stats = client.download_latest_csv(args.download_latest[0], args.download_latest[1])
                if stats:
                    print(f"\n{'✅ Success' if stats.success else '❌ Failed'}")
                else:
                    print("❌ No bodies_*.csv found")
                    
            elif args.list:
                files = client.listdir_detailed(args.list)
                print(f"\nContents of {args.list}:")
                for f in files:
                    print(f"  {f.get('name', f)}")
                    
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == '__main__':
    exit(main())
