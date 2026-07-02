"""
Temporary file management utilities.
Provides context managers for safe temporary file handling.
"""

import tempfile
import shutil
import os
from contextlib import contextmanager


@contextmanager
def temp_file_manager(cleanup=True):
    """
    Context manager for temporary directory management.
    
    Args:
        cleanup: Whether to clean up temp directory on exit (default: True)
        
    Yields:
        str: Path to temporary directory
        
    Example:
        with temp_file_manager() as temp_dir:
            file_path = os.path.join(temp_dir, 'test.txt')
            # Use file_path
        # Directory automatically cleaned up
    """
    temp_dir = tempfile.mkdtemp()
    try:
        yield temp_dir
    finally:
        if cleanup:
            try:
                shutil.rmtree(temp_dir, ignore_errors=True)
            except Exception as e:
                import logging
                logging.getLogger(__name__).debug("Temp dir cleanup failed: %s", e)
