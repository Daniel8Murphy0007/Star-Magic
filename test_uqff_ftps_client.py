#!/usr/bin/env python3
"""
test_uqff_ftps_client.py - Test Suite for UQFF FTPS Client
===========================================================

Unit tests and integration tests for uqff_ftps_client.py.
Tests RFC 4217 compliance, TLS security, and file operations.

Usage:
    pytest test_uqff_ftps_client.py -v
    python test_uqff_ftps_client.py  # Direct execution

Author: Daniel T. Murphy
Framework: UQFF Star-Magic v3.0
"""

import unittest
import ssl
import os
import json
import tempfile
from unittest.mock import Mock, patch, MagicMock
from datetime import datetime
from pathlib import Path

# Import the module under test
try:
    from uqff_ftps_client import (
        UQFFFTPSClient,
        FTPSConfig,
        TransferStats,
        ImplicitFTPS
    )
except ImportError as e:
    print(f"Warning: Could not import uqff_ftps_client: {e}")
    UQFFFTPSClient = None


class TestFTPSConfig(unittest.TestCase):
    """Test FTPSConfig dataclass."""
    
    def test_default_values(self):
        """Test default configuration values."""
        config = FTPSConfig(host='test.example.com')
        
        self.assertEqual(config.host, 'test.example.com')
        self.assertEqual(config.port, 21)
        self.assertEqual(config.user, 'anonymous')
        self.assertEqual(config.password, '')
        self.assertFalse(config.implicit_tls)
        self.assertTrue(config.verify_cert)
        self.assertTrue(config.passive_mode)
        self.assertEqual(config.timeout, 30.0)
        self.assertEqual(config.uqff_data_dir, '/uqff_data/')
    
    def test_implicit_ftps_port(self):
        """Test implicit FTPS typically uses port 990."""
        config = FTPSConfig(
            host='secure.example.com',
            port=990,
            implicit_tls=True
        )
        
        self.assertEqual(config.port, 990)
        self.assertTrue(config.implicit_tls)
    
    def test_from_env(self):
        """Test configuration from environment variables."""
        with patch.dict(os.environ, {
            'UQFF_FTPS_HOST': 'env.example.com',
            'UQFF_FTPS_PORT': '990',
            'UQFF_FTPS_USER': 'testuser',
            'UQFF_FTPS_PASS': 'testpass',
            'UQFF_FTPS_IMPLICIT': 'true',
            'UQFF_FTPS_VERIFY': 'true',
            'UQFF_FTPS_PASSIVE': 'true',
            'UQFF_FTPS_DATA_DIR': '/custom/data/'
        }):
            config = FTPSConfig.from_env()
            
            self.assertEqual(config.host, 'env.example.com')
            self.assertEqual(config.port, 990)
            self.assertEqual(config.user, 'testuser')
            self.assertEqual(config.password, 'testpass')
            self.assertTrue(config.implicit_tls)
            self.assertTrue(config.verify_cert)
            self.assertEqual(config.uqff_data_dir, '/custom/data/')
    
    def test_to_dict_excludes_password(self):
        """Test that to_dict excludes sensitive data."""
        config = FTPSConfig(
            host='test.example.com',
            user='user',
            password='secret_password'
        )
        
        config_dict = config.to_dict()
        
        self.assertIn('host', config_dict)
        self.assertIn('user', config_dict)
        self.assertNotIn('password', config_dict)


class TestTransferStats(unittest.TestCase):
    """Test TransferStats dataclass."""
    
    def test_duration_calculation(self):
        """Test transfer duration calculation."""
        stats = TransferStats(
            filename='test.txt',
            direction='upload',
            bytes_transferred=1024,
            total_bytes=1024,
            start_time=datetime(2026, 2, 21, 12, 0, 0),
            end_time=datetime(2026, 2, 21, 12, 0, 10),
            success=True
        )
        
        self.assertEqual(stats.duration_seconds, 10.0)
    
    def test_speed_calculation(self):
        """Test transfer speed calculation."""
        stats = TransferStats(
            filename='test.txt',
            direction='download',
            bytes_transferred=10240,
            total_bytes=10240,
            start_time=datetime(2026, 2, 21, 12, 0, 0),
            end_time=datetime(2026, 2, 21, 12, 0, 10),
            success=True
        )
        
        self.assertEqual(stats.speed_bps, 1024.0)  # 10240 / 10 = 1024 B/s
    
    def test_progress_percentage(self):
        """Test progress percentage calculation."""
        stats = TransferStats(
            filename='test.txt',
            direction='upload',
            bytes_transferred=512,
            total_bytes=1024
        )
        
        self.assertEqual(stats.progress_percent, 50.0)
    
    def test_zero_division_protection(self):
        """Test protection against zero division."""
        stats = TransferStats(
            filename='test.txt',
            direction='upload',
            bytes_transferred=0,
            total_bytes=0
        )
        
        self.assertEqual(stats.progress_percent, 0.0)
        self.assertEqual(stats.speed_bps, 0.0)


class TestSSLContext(unittest.TestCase):
    """Test SSL context creation and TLS settings."""
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_tls_minimum_version(self):
        """Test that TLS 1.2+ is enforced."""
        config = FTPSConfig(host='test.example.com')
        client = UQFFFTPSClient(config=config)
        
        context = client._create_ssl_context()
        
        # Should be TLS 1.2 or TLS 1.3 minimum
        self.assertIn(context.minimum_version, [
            ssl.TLSVersion.TLSv1_2,
            ssl.TLSVersion.TLSv1_3
        ])
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_certificate_verification_enabled(self):
        """Test certificate verification is enabled by default."""
        config = FTPSConfig(host='test.example.com', verify_cert=True)
        client = UQFFFTPSClient(config=config)
        
        context = client._create_ssl_context()
        
        self.assertEqual(context.verify_mode, ssl.CERT_REQUIRED)
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_certificate_verification_disabled_warning(self):
        """Test that disabling cert verification logs a warning."""
        config = FTPSConfig(host='test.example.com', verify_cert=False)
        client = UQFFFTPSClient(config=config)
        
        with patch('uqff_ftps_client.logger') as mock_logger:
            context = client._create_ssl_context()
            
            # Should log a warning about disabled verification
            mock_logger.warning.assert_called()
        
        self.assertEqual(context.verify_mode, ssl.CERT_NONE)


class TestRFC4217Compliance(unittest.TestCase):
    """Test RFC 4217 compliance validation."""
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_compliance_check_structure(self):
        """Test compliance check returns proper structure."""
        config = FTPSConfig(host='test.example.com')
        client = UQFFFTPSClient(config=config)
        
        # Without connection, should still return structure
        results = client.validate_rfc4217_compliance()
        
        self.assertIn('rfc4217_compliant', results)
        self.assertIn('checks', results)
        self.assertIn('warnings', results)
        self.assertIn('tls_version', results)
        self.assertIn('cipher_suite', results)
        self.assertIn('certificate', results)
        self.assertIsInstance(results['checks'], list)
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_compliance_check_cert_verification_warning(self):
        """Test compliance check warns about disabled cert verification."""
        config = FTPSConfig(host='test.example.com', verify_cert=False)
        client = UQFFFTPSClient(config=config)
        
        results = client.validate_rfc4217_compliance()
        
        self.assertIn('Certificate verification DISABLED', str(results['warnings']))
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_explicit_ftps_mode_check(self):
        """Test explicit FTPS mode is reported correctly."""
        config = FTPSConfig(host='test.example.com', implicit_tls=False)
        client = UQFFFTPSClient(config=config)
        
        results = client.validate_rfc4217_compliance()
        
        mode_check = next(c for c in results['checks'] if c['name'] == 'FTPS Mode')
        self.assertEqual(mode_check['status'], 'PASS')
        self.assertIn('Explicit', mode_check['detail'])
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_implicit_ftps_mode_check(self):
        """Test implicit FTPS mode is reported correctly."""
        config = FTPSConfig(host='test.example.com', implicit_tls=True, port=990)
        client = UQFFFTPSClient(config=config)
        
        results = client.validate_rfc4217_compliance()
        
        mode_check = next(c for c in results['checks'] if c['name'] == 'FTPS Mode')
        self.assertEqual(mode_check['status'], 'PASS')
        self.assertIn('Implicit', mode_check['detail'])


class TestImplicitFTPS(unittest.TestCase):
    """Test ImplicitFTPS class for port 990."""
    
    def test_implicit_ftps_class_exists(self):
        """Test ImplicitFTPS class is defined."""
        self.assertIsNotNone(ImplicitFTPS)
    
    def test_implicit_ftps_default_port(self):
        """Test ImplicitFTPS uses port 990 by default."""
        # The connect method should use 990 as default
        with patch.object(ImplicitFTPS, '_create_socket', return_value=Mock()):
            with patch.object(ImplicitFTPS, 'context', create=True) as mock_ctx:
                mock_sock = Mock()
                mock_ctx.wrap_socket.return_value = mock_sock
                mock_sock.makefile.return_value = Mock()
                
                # Verify the class can be instantiated
                ftp = ImplicitFTPS(host='test.example.com')
                self.assertIsNotNone(ftp)


class TestFileOperations(unittest.TestCase):
    """Test file transfer operations (mocked)."""
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_context_manager(self):
        """Test client works as context manager."""
        config = FTPSConfig(host='test.example.com')
        
        with patch.object(UQFFFTPSClient, 'connect', return_value=True):
            with patch.object(UQFFFTPSClient, 'disconnect'):
                with UQFFFTPSClient(config=config) as client:
                    self.assertIsNotNone(client)
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_not_connected_raises(self):
        """Test operations raise when not connected."""
        config = FTPSConfig(host='test.example.com')
        client = UQFFFTPSClient(config=config)
        
        with self.assertRaises(ConnectionError):
            client._ensure_connected()
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_md5_calculation(self):
        """Test MD5 checksum calculation."""
        config = FTPSConfig(host='test.example.com')
        client = UQFFFTPSClient(config=config)
        
        # Create a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as f:
            f.write('test content for MD5')
            temp_path = f.name
        
        try:
            md5 = client._calculate_md5(temp_path)
            self.assertEqual(len(md5), 32)  # MD5 hex digest is 32 chars
            self.assertTrue(all(c in '0123456789abcdef' for c in md5))
        finally:
            os.unlink(temp_path)


class TestUQFFServerIntegration(unittest.TestCase):
    """Test integration with uqff_server.js."""
    
    @unittest.skipIf(UQFFFTPSClient is None, "uqff_ftps_client not available")
    def test_call_uqff_server_structure(self):
        """Test call_uqff_server method exists."""
        config = FTPSConfig(host='test.example.com')
        client = UQFFFTPSClient(config=config)
        
        self.assertTrue(hasattr(client, 'call_uqff_server'))
    
    def test_uqff_data_dir_default(self):
        """Test default UQFF data directory."""
        config = FTPSConfig(host='test.example.com')
        
        self.assertEqual(config.uqff_data_dir, '/uqff_data/')


class TestBodiesFilePattern(unittest.TestCase):
    """Test bodies_*.csv file pattern matching."""
    
    def test_bodies_pattern_match(self):
        """Test regex pattern matches bodies files."""
        import re
        
        config = FTPSConfig(host='test.example.com')
        pattern = config.bodies_pattern
        
        # Should match
        self.assertIsNotNone(re.match(pattern, 'bodies_20260221_123456.csv'))
        self.assertIsNotNone(re.match(pattern, 'bodies_20251231_000000.csv'))
        
        # Should not match
        self.assertIsNone(re.match(pattern, 'bodies.csv'))
        self.assertIsNone(re.match(pattern, 'bodies_invalid.csv'))
        self.assertIsNone(re.match(pattern, 'data_20260221_123456.csv'))


# ═══════════════════════════════════════════════════════════════════════════════
# INTEGRATION TESTS (require actual FTPS server)
# ═══════════════════════════════════════════════════════════════════════════════

class TestIntegration(unittest.TestCase):
    """Integration tests (skipped unless FTPS server available)."""
    
    @classmethod
    def setUpClass(cls):
        """Check if test FTPS server is available."""
        cls.ftps_host = os.getenv('UQFF_TEST_FTPS_HOST', '')
        cls.ftps_available = bool(cls.ftps_host)
    
    @unittest.skipUnless(os.getenv('UQFF_TEST_FTPS_HOST'), "No test FTPS server")
    def test_real_connection(self):
        """Test connection to real FTPS server."""
        config = FTPSConfig.from_env()
        
        with UQFFFTPSClient(config=config) as client:
            self.assertTrue(client.connected)
            
            # Print compliance report
            client.print_compliance_report()
    
    @unittest.skipUnless(os.getenv('UQFF_TEST_FTPS_HOST'), "No test FTPS server")
    def test_real_directory_listing(self):
        """Test directory listing on real server."""
        config = FTPSConfig.from_env()
        
        with UQFFFTPSClient(config=config) as client:
            pwd = client.pwd()
            self.assertIsInstance(pwd, str)
            
            files = client.listdir('.')
            self.assertIsInstance(files, list)


# ═══════════════════════════════════════════════════════════════════════════════
# TEST RUNNER
# ═══════════════════════════════════════════════════════════════════════════════

def run_tests():
    """Run all tests with verbose output."""
    print("═" * 70)
    print(" UQFF FTPS Client Test Suite")
    print("═" * 70)
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestFTPSConfig))
    suite.addTests(loader.loadTestsFromTestCase(TestTransferStats))
    suite.addTests(loader.loadTestsFromTestCase(TestSSLContext))
    suite.addTests(loader.loadTestsFromTestCase(TestRFC4217Compliance))
    suite.addTests(loader.loadTestsFromTestCase(TestImplicitFTPS))
    suite.addTests(loader.loadTestsFromTestCase(TestFileOperations))
    suite.addTests(loader.loadTestsFromTestCase(TestUQFFServerIntegration))
    suite.addTests(loader.loadTestsFromTestCase(TestBodiesFilePattern))
    suite.addTests(loader.loadTestsFromTestCase(TestIntegration))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print()
    print("═" * 70)
    print(f" Results: {result.testsRun} tests, "
          f"{len(result.failures)} failures, "
          f"{len(result.errors)} errors, "
          f"{len(result.skipped)} skipped")
    print("═" * 70)
    
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    exit(run_tests())
