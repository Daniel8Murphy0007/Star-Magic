"""
file_lock.py - Lock-Free File Access Patterns
==============================================

Provides safe file access patterns for concurrent access from
multiple programs in the UQFF pipeline.

Patterns:
- Advisory file locking (cross-platform)
- Atomic write-then-rename
- Safe concurrent reads
- Lockfile-based serialization

Usage:
    from ipc.file_lock import FileLock, atomic_write, safe_read
    
    # Atomic write (safe for concurrent readers)
    atomic_write('data.json', json.dumps(data))
    
    # Safe read (handles partial writes)
    data = safe_read('data.json')
    
    # Exclusive lock for writes
    with FileLock('important.csv') as lock:
        lock.write(csv_data)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

import os
import sys
import time
import logging
import tempfile
import threading
from pathlib import Path
from typing import Optional, Union, Any
from contextlib import contextmanager

logger = logging.getLogger("UQFF_FileLock")


class FileLock:
    """
    Cross-platform advisory file lock
    
    Uses:
    - Windows: msvcrt.locking()
    - Unix: fcntl.flock()
    
    Falls back to lockfile pattern if native locking unavailable.
    """
    
    def __init__(self, 
                 filepath: Union[str, Path],
                 timeout: float = 10.0,
                 mode: str = 'exclusive'):
        """
        Initialize FileLock
        
        Args:
            filepath: Path to file to lock
            timeout: Lock acquisition timeout in seconds
            mode: 'exclusive' (write) or 'shared' (read)
        """
        self.filepath = Path(filepath)
        self.timeout = timeout
        self.mode = mode
        
        self._lockfile = self.filepath.with_suffix(self.filepath.suffix + '.lock')
        self._fd: Optional[int] = None
        self._locked = False
        self._use_native = self._check_native_support()
    
    def _check_native_support(self) -> bool:
        """Check if native locking is available"""
        if sys.platform == 'win32':
            try:
                import msvcrt
                return True
            except ImportError:
                return False
        else:
            try:
                import fcntl
                return True
            except ImportError:
                return False
    
    def acquire(self) -> bool:
        """
        Acquire the file lock
        
        Returns:
            True if lock acquired, False on timeout
        """
        if self._locked:
            return True
        
        if self._use_native:
            return self._acquire_native()
        else:
            return self._acquire_lockfile()
    
    def release(self) -> None:
        """Release the file lock"""
        if not self._locked:
            return
        
        if self._use_native:
            self._release_native()
        else:
            self._release_lockfile()
        
        self._locked = False
    
    def _acquire_native(self) -> bool:
        """Acquire using native OS locking"""
        start = time.time()
        
        while time.time() - start < self.timeout:
            try:
                # Open or create the file
                if not self.filepath.exists():
                    self.filepath.touch()
                
                if sys.platform == 'win32':
                    import msvcrt
                    self._fd = os.open(str(self.filepath), os.O_RDWR)
                    
                    lock_mode = msvcrt.LK_NBLCK if self.mode == 'exclusive' else msvcrt.LK_NBRLCK
                    try:
                        msvcrt.locking(self._fd, lock_mode, 1)
                        self._locked = True
                        return True
                    except IOError:
                        os.close(self._fd)
                        self._fd = None
                else:
                    import fcntl
                    self._fd = os.open(str(self.filepath), os.O_RDWR)
                    
                    lock_mode = fcntl.LOCK_EX if self.mode == 'exclusive' else fcntl.LOCK_SH
                    lock_mode |= fcntl.LOCK_NB  # Non-blocking
                    
                    try:
                        fcntl.flock(self._fd, lock_mode)
                        self._locked = True
                        return True
                    except IOError:
                        os.close(self._fd)
                        self._fd = None
                
            except Exception as e:
                logger.debug(f"Lock attempt failed: {e}")
            
            time.sleep(0.1)
        
        logger.warning(f"Lock acquisition timeout: {self.filepath}")
        return False
    
    def _release_native(self) -> None:
        """Release native OS lock"""
        if self._fd is None:
            return
        
        try:
            if sys.platform == 'win32':
                import msvcrt
                msvcrt.locking(self._fd, msvcrt.LK_UNLCK, 1)
            else:
                import fcntl
                fcntl.flock(self._fd, fcntl.LOCK_UN)
            
            os.close(self._fd)
        except Exception as e:
            logger.error(f"Lock release error: {e}")
        finally:
            self._fd = None
    
    def _acquire_lockfile(self) -> bool:
        """Acquire using lockfile pattern"""
        start = time.time()
        
        while time.time() - start < self.timeout:
            try:
                # Create lockfile atomically
                fd = os.open(
                    str(self._lockfile),
                    os.O_CREAT | os.O_EXCL | os.O_WRONLY
                )
                os.write(fd, f"{os.getpid()}:{threading.current_thread().ident}\n".encode())
                os.close(fd)
                self._locked = True
                return True
                
            except FileExistsError:
                # Check if lock is stale (holder dead)
                if self._is_stale_lock():
                    try:
                        self._lockfile.unlink()
                    except:
                        pass
                else:
                    time.sleep(0.1)
            except Exception as e:
                logger.debug(f"Lockfile attempt failed: {e}")
                time.sleep(0.1)
        
        return False
    
    def _release_lockfile(self) -> None:
        """Release lockfile"""
        try:
            self._lockfile.unlink()
        except Exception as e:
            logger.error(f"Lockfile release error: {e}")
    
    def _is_stale_lock(self) -> bool:
        """Check if existing lock is stale"""
        try:
            if not self._lockfile.exists():
                return False
            
            # Check lock age
            age = time.time() - self._lockfile.stat().st_mtime
            if age > 60:  # Stale if older than 60 seconds
                return True
            
            # Check if PID is alive
            content = self._lockfile.read_text().strip()
            if ':' in content:
                pid = int(content.split(':')[0])
                
                if sys.platform == 'win32':
                    import ctypes
                    kernel32 = ctypes.windll.kernel32
                    handle = kernel32.OpenProcess(0x1000, 0, pid)
                    if handle == 0:
                        return True
                    kernel32.CloseHandle(handle)
                else:
                    try:
                        os.kill(pid, 0)  # Check if process exists
                    except OSError:
                        return True
            
            return False
            
        except Exception:
            return True  # Assume stale on error
    
    def __enter__(self):
        if not self.acquire():
            raise TimeoutError(f"Could not acquire lock: {self.filepath}")
        return self
    
    def __exit__(self, *args):
        self.release()


def atomic_write(filepath: Union[str, Path], 
                 content: Union[str, bytes],
                 encoding: str = 'utf-8') -> bool:
    """
    Write file atomically using temp file + rename
    
    Safe for concurrent readers - they will either see the old
    or new file, never a partial write.
    
    Args:
        filepath: Target file path
        content: Content to write (str or bytes)
        encoding: Encoding for string content
    
    Returns:
        True if write succeeded
    """
    filepath = Path(filepath)
    
    try:
        # Write to temp file in same directory
        fd, tmp_path = tempfile.mkstemp(
            suffix='.tmp',
            prefix=filepath.stem + '_',
            dir=filepath.parent
        )
        
        try:
            if isinstance(content, str):
                os.write(fd, content.encode(encoding))
            else:
                os.write(fd, content)
            os.close(fd)
            
            # Atomic rename (same filesystem)
            Path(tmp_path).replace(filepath)
            return True
            
        except Exception:
            os.close(fd)
            try:
                os.unlink(tmp_path)
            except:
                pass
            raise
            
    except Exception as e:
        logger.error(f"Atomic write failed: {e}")
        return False


def safe_read(filepath: Union[str, Path],
              encoding: str = 'utf-8',
              default: Any = None,
              retries: int = 3) -> Any:
    """
    Read file safely, handling concurrent writes
    
    Retries on read failures (e.g., during atomic rename).
    Returns default if file doesn't exist.
    
    Args:
        filepath: File path to read
        encoding: Encoding for text files
        default: Default value if file missing/unreadable
        retries: Number of retry attempts
    
    Returns:
        File content or default
    """
    filepath = Path(filepath)
    
    if not filepath.exists():
        return default
    
    for attempt in range(retries):
        try:
            with open(filepath, 'r', encoding=encoding) as f:
                return f.read()
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(0.1)  # Brief delay before retry
            else:
                logger.error(f"Safe read failed after {retries} attempts: {e}")
                return default
    
    return default


def safe_read_json(filepath: Union[str, Path],
                   default: Any = None) -> Any:
    """
    Read and parse JSON file safely
    
    Args:
        filepath: Path to JSON file
        default: Default value if missing/invalid
    
    Returns:
        Parsed JSON or default
    """
    import json
    
    content = safe_read(filepath, default=None)
    if content is None:
        return default
    
    try:
        return json.loads(content)
    except json.JSONDecodeError as e:
        logger.error(f"JSON parse error: {e}")
        return default


def atomic_write_json(filepath: Union[str, Path],
                      data: Any,
                      indent: int = 2) -> bool:
    """
    Write JSON file atomically
    
    Args:
        filepath: Target file path
        data: Data to serialize
        indent: JSON indentation
    
    Returns:
        True if write succeeded
    """
    import json
    
    try:
        content = json.dumps(data, indent=indent, default=str)
        return atomic_write(filepath, content)
    except Exception as e:
        logger.error(f"JSON write error: {e}")
        return False


@contextmanager
def locked_file(filepath: Union[str, Path],
                mode: str = 'r',
                timeout: float = 10.0):
    """
    Context manager for locked file access
    
    Usage:
        with locked_file('data.csv', 'w') as f:
            f.write(csv_data)
    
    Args:
        filepath: File to access
        mode: File open mode
        timeout: Lock timeout
    """
    filepath = Path(filepath)
    lock_mode = 'exclusive' if 'w' in mode or 'a' in mode else 'shared'
    
    lock = FileLock(filepath, timeout=timeout, mode=lock_mode)
    
    if not lock.acquire():
        raise TimeoutError(f"Could not acquire lock: {filepath}")
    
    try:
        with open(filepath, mode) as f:
            yield f
    finally:
        lock.release()


# ═══════════════════════════════════════════════════════════════════════════════
# TEST
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    import json
    
    print("FileLock Test")
    print("=" * 40)
    
    test_file = Path(".test_lock_file.json")
    
    # Test atomic write
    data = {'test': 'value', 'timestamp': time.time()}
    print(f"Atomic write: {atomic_write_json(test_file, data)}")
    
    # Test safe read
    read_data = safe_read_json(test_file)
    print(f"Safe read: {read_data}")
    
    # Test file lock
    with FileLock(test_file) as lock:
        print("Lock acquired, doing work...")
        time.sleep(0.5)
    print("Lock released")
    
    # Cleanup
    test_file.unlink()
    
    print("\nTest complete")
