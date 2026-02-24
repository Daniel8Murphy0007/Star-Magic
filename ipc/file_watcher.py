"""
file_watcher.py - File System Watcher for UQFF Pipeline
========================================================

Watches for changes to bodies_*.csv files from APIFetch.py and triggers
recalculation in the physics calculators.

Usage:
    from ipc.file_watcher import BodiesFileWatcher, watch_bodies_directory
    
    # Simple usage
    watcher = watch_bodies_directory(callback=on_new_file)
    
    # Advanced usage
    watcher = BodiesFileWatcher()
    watcher.on_new_file = lambda path: print(f"New file: {path}")
    watcher.start()

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF Star-Magic v3.0
"""

import os
import sys
import time
import logging
import threading
import fnmatch
from pathlib import Path
from typing import Callable, Optional, List, Set
from datetime import datetime
from dataclasses import dataclass, field

logger = logging.getLogger("UQFF_FileWatcher")


@dataclass
class FileEvent:
    """Represents a file system event"""
    path: Path
    event_type: str  # 'created', 'modified', 'deleted'
    timestamp: datetime = field(default_factory=datetime.now)
    size: int = 0


class FileWatcher:
    """
    Cross-platform file watcher using polling
    
    For production, consider watchdog library for efficiency,
    but this works without external dependencies.
    """
    
    def __init__(self, 
                 directory: str = ".",
                 patterns: List[str] = None,
                 poll_interval: float = 1.0,
                 recursive: bool = False):
        """
        Initialize FileWatcher
        
        Args:
            directory: Directory to watch
            patterns: File patterns to match (e.g., ['bodies_*.csv'])
            poll_interval: Polling interval in seconds
            recursive: Watch subdirectories
        """
        self.directory = Path(directory).resolve()
        self.patterns = patterns or ['*']
        self.poll_interval = poll_interval
        self.recursive = recursive
        
        self._running = False
        self._thread: Optional[threading.Thread] = None
        self._known_files: dict = {}  # path -> (mtime, size)
        self._callbacks: List[Callable[[FileEvent], None]] = []
        
        # Specific callbacks
        self.on_created: Optional[Callable[[Path], None]] = None
        self.on_modified: Optional[Callable[[Path], None]] = None
        self.on_deleted: Optional[Callable[[Path], None]] = None
    
    def add_callback(self, callback: Callable[[FileEvent], None]) -> None:
        """Add a callback for all file events"""
        self._callbacks.append(callback)
    
    def _matches_pattern(self, filename: str) -> bool:
        """Check if filename matches any pattern"""
        return any(fnmatch.fnmatch(filename, p) for p in self.patterns)
    
    def _scan_directory(self) -> dict:
        """Scan directory and return file info dict"""
        files = {}
        
        if self.recursive:
            for root, dirs, filenames in os.walk(self.directory):
                for name in filenames:
                    if self._matches_pattern(name):
                        path = Path(root) / name
                        try:
                            stat = path.stat()
                            files[str(path)] = (stat.st_mtime, stat.st_size)
                        except OSError:
                            pass
        else:
            for path in self.directory.iterdir():
                if path.is_file() and self._matches_pattern(path.name):
                    try:
                        stat = path.stat()
                        files[str(path)] = (stat.st_mtime, stat.st_size)
                    except OSError:
                        pass
        
        return files
    
    def _emit_event(self, event: FileEvent) -> None:
        """Emit file event to callbacks"""
        # Call specific callbacks
        if event.event_type == 'created' and self.on_created:
            self.on_created(event.path)
        elif event.event_type == 'modified' and self.on_modified:
            self.on_modified(event.path)
        elif event.event_type == 'deleted' and self.on_deleted:
            self.on_deleted(event.path)
        
        # Call general callbacks
        for callback in self._callbacks:
            try:
                callback(event)
            except Exception as e:
                logger.error(f"Callback error: {e}")
    
    def _poll_loop(self) -> None:
        """Main polling loop"""
        while self._running:
            try:
                current_files = self._scan_directory()
                
                # Check for new/modified files
                for path, (mtime, size) in current_files.items():
                    if path not in self._known_files:
                        # New file
                        event = FileEvent(
                            path=Path(path),
                            event_type='created',
                            size=size
                        )
                        self._emit_event(event)
                        logger.info(f"New file: {path}")
                    elif self._known_files[path] != (mtime, size):
                        # Modified file
                        event = FileEvent(
                            path=Path(path),
                            event_type='modified',
                            size=size
                        )
                        self._emit_event(event)
                        logger.debug(f"Modified: {path}")
                
                # Check for deleted files
                for path in list(self._known_files.keys()):
                    if path not in current_files:
                        event = FileEvent(
                            path=Path(path),
                            event_type='deleted'
                        )
                        self._emit_event(event)
                        logger.info(f"Deleted: {path}")
                
                self._known_files = current_files
                
            except Exception as e:
                logger.error(f"Poll error: {e}")
            
            time.sleep(self.poll_interval)
    
    def start(self) -> None:
        """Start watching in background thread"""
        if self._running:
            return
        
        # Initial scan
        self._known_files = self._scan_directory()
        
        self._running = True
        self._thread = threading.Thread(target=self._poll_loop, daemon=True)
        self._thread.start()
        logger.info(f"FileWatcher started: {self.directory} ({self.patterns})")
    
    def stop(self) -> None:
        """Stop watching"""
        self._running = False
        if self._thread:
            self._thread.join(timeout=2.0)
        logger.info("FileWatcher stopped")
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, *args):
        self.stop()


class BodiesFileWatcher(FileWatcher):
    """
    Specialized watcher for bodies_*.csv files from APIFetch.py
    
    Automatically triggers physics recalculation when new data arrives.
    """
    
    def __init__(self, 
                 directory: str = ".",
                 poll_interval: float = 2.0):
        super().__init__(
            directory=directory,
            patterns=['bodies_*.csv'],
            poll_interval=poll_interval,
            recursive=False
        )
        
        # Physics callback
        self.on_new_data: Optional[Callable[[Path, dict], None]] = None
        
        # Override on_created
        self.on_created = self._handle_new_bodies_file
    
    def _parse_bodies_csv(self, path: Path) -> dict:
        """Parse bodies CSV file into parameters dict"""
        import csv
        
        params = {}
        try:
            with open(path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                rows = list(reader)
                if rows:
                    params = {
                        'file': str(path),
                        'timestamp': path.stem.split('_', 1)[1] if '_' in path.stem else '',
                        'body_count': len(rows),
                        'bodies': rows
                    }
        except Exception as e:
            logger.error(f"Failed to parse {path}: {e}")
        
        return params
    
    def _handle_new_bodies_file(self, path: Path) -> None:
        """Handle new bodies_*.csv file"""
        logger.info(f"New bodies file detected: {path}")
        
        if self.on_new_data:
            params = self._parse_bodies_csv(path)
            if params:
                self.on_new_data(path, params)


def watch_bodies_directory(
    directory: str = ".",
    callback: Callable[[Path, dict], None] = None,
    poll_interval: float = 2.0
) -> BodiesFileWatcher:
    """
    Convenience function to watch for bodies_*.csv files
    
    Args:
        directory: Directory to watch (default: current)
        callback: Function called with (path, params) on new file
        poll_interval: Polling interval in seconds
    
    Returns:
        Started BodiesFileWatcher instance
    
    Example:
        def on_new_bodies(path, params):
            print(f"New data: {params['body_count']} bodies")
            # Trigger physics calculation...
        
        watcher = watch_bodies_directory(callback=on_new_bodies)
    """
    watcher = BodiesFileWatcher(directory, poll_interval)
    watcher.on_new_data = callback
    watcher.start()
    return watcher


# ═══════════════════════════════════════════════════════════════════════════════
# TEST
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    def on_event(event: FileEvent):
        print(f"[{event.event_type}] {event.path} ({event.size} bytes)")
    
    print("FileWatcher Test - watching current directory for *.csv")
    print("Press Ctrl+C to stop")
    
    watcher = FileWatcher(patterns=['*.csv'])
    watcher.add_callback(on_event)
    
    try:
        watcher.start()
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        watcher.stop()
        print("Stopped")
