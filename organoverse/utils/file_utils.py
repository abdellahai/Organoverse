"""
File handling utilities for Organoverse workbench.
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional, List, Union
import gzip
import hashlib

from .exceptions import OrganoverseError


def create_temp_dir(prefix: str = "organoverse_", base_dir: Optional[str] = None) -> str:
    """
    Create a temporary directory.
    
    Args:
        prefix: Prefix for directory name
        base_dir: Base directory for temp dir
    
    Returns:
        Path to created temporary directory
    """
    try:
        temp_dir = tempfile.mkdtemp(prefix=prefix, dir=base_dir)
        return temp_dir
    except Exception as e:
        raise OrganoverseError(f"Failed to create temporary directory: {str(e)}")


def cleanup_temp_files(temp_paths: List[str], force: bool = False):
    """
    Clean up temporary files and directories.
    
    Args:
        temp_paths: List of paths to clean up
        force: Force removal even if errors occur
    """
    for path in temp_paths:
        try:
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                shutil.rmtree(path)
        except Exception as e:
            if not force:
                raise OrganoverseError(f"Failed to cleanup {path}: {str(e)}")


def ensure_directory(directory: Union[str, Path]) -> Path:
    """
    Ensure directory exists, create if necessary.
    
    Args:
        directory: Directory path
    
    Returns:
        Path object
    """
    dir_path = Path(directory)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def copy_file(src: str, dst: str, overwrite: bool = False) -> str:
    """
    Copy file with error handling.
    
    Args:
        src: Source file path
        dst: Destination file path
        overwrite: Whether to overwrite existing files
    
    Returns:
        Destination file path
    """
    src_path = Path(src)
    dst_path = Path(dst)
    
    if not src_path.exists():
        raise OrganoverseError(f"Source file does not exist: {src}")
    
    if dst_path.exists() and not overwrite:
        raise OrganoverseError(f"Destination file exists: {dst}")
    
    try:
        # Ensure destination directory exists
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Copy file
        shutil.copy2(src, dst)
        return str(dst_path)
        
    except Exception as e:
        raise OrganoverseError(f"Failed to copy file: {str(e)}")


def move_file(src: str, dst: str, overwrite: bool = False) -> str:
    """
    Move file with error handling.
    
    Args:
        src: Source file path
        dst: Destination file path
        overwrite: Whether to overwrite existing files
    
    Returns:
        Destination file path
    """
    src_path = Path(src)
    dst_path = Path(dst)
    
    if not src_path.exists():
        raise OrganoverseError(f"Source file does not exist: {src}")
    
    if dst_path.exists() and not overwrite:
        raise OrganoverseError(f"Destination file exists: {dst}")
    
    try:
        # Ensure destination directory exists
        dst_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Move file
        shutil.move(src, dst)
        return str(dst_path)
        
    except Exception as e:
        raise OrganoverseError(f"Failed to move file: {str(e)}")


def get_file_size(file_path: str) -> int:
    """
    Get file size in bytes.
    
    Args:
        file_path: Path to file
    
    Returns:
        File size in bytes
    """
    try:
        return os.path.getsize(file_path)
    except Exception as e:
        raise OrganoverseError(f"Failed to get file size: {str(e)}")


def get_file_hash(file_path: str, algorithm: str = "md5") -> str:
    """
    Calculate file hash.
    
    Args:
        file_path: Path to file
        algorithm: Hash algorithm (md5, sha1, sha256)
    
    Returns:
        File hash as hexadecimal string
    """
    try:
        hash_obj = hashlib.new(algorithm)
        
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_obj.update(chunk)
        
        return hash_obj.hexdigest()
        
    except Exception as e:
        raise OrganoverseError(f"Failed to calculate file hash: {str(e)}")


def compress_file(file_path: str, output_path: Optional[str] = None, remove_original: bool = False) -> str:
    """
    Compress file using gzip.
    
    Args:
        file_path: Path to file to compress
        output_path: Output path (defaults to input + .gz)
        remove_original: Whether to remove original file
    
    Returns:
        Path to compressed file
    """
    if output_path is None:
        output_path = file_path + ".gz"
    
    try:
        with open(file_path, 'rb') as f_in:
            with gzip.open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        if remove_original:
            os.remove(file_path)
        
        return output_path
        
    except Exception as e:
        raise OrganoverseError(f"Failed to compress file: {str(e)}")


def decompress_file(file_path: str, output_path: Optional[str] = None, remove_original: bool = False) -> str:
    """
    Decompress gzip file.
    
    Args:
        file_path: Path to compressed file
        output_path: Output path (defaults to input without .gz)
        remove_original: Whether to remove original file
    
    Returns:
        Path to decompressed file
    """
    if output_path is None:
        if file_path.endswith('.gz'):
            output_path = file_path[:-3]
        else:
            output_path = file_path + ".decompressed"
    
    try:
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        if remove_original:
            os.remove(file_path)
        
        return output_path
        
    except Exception as e:
        raise OrganoverseError(f"Failed to decompress file: {str(e)}")


def is_compressed(file_path: str) -> bool:
    """
    Check if file is gzip compressed.
    
    Args:
        file_path: Path to file
    
    Returns:
        True if file is compressed
    """
    try:
        with open(file_path, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'
    except Exception:
        return False


def get_available_space(directory: str) -> int:
    """
    Get available disk space in bytes.
    
    Args:
        directory: Directory path
    
    Returns:
        Available space in bytes
    """
    try:
        statvfs = os.statvfs(directory)
        return statvfs.f_frsize * statvfs.f_bavail
    except Exception as e:
        raise OrganoverseError(f"Failed to get available space: {str(e)}")


def find_files(directory: str, pattern: str = "*", recursive: bool = True) -> List[str]:
    """
    Find files matching pattern.
    
    Args:
        directory: Directory to search
        pattern: File pattern (glob style)
        recursive: Whether to search recursively
    
    Returns:
        List of matching file paths
    """
    try:
        dir_path = Path(directory)
        
        if recursive:
            files = list(dir_path.rglob(pattern))
        else:
            files = list(dir_path.glob(pattern))
        
        return [str(f) for f in files if f.is_file()]
        
    except Exception as e:
        raise OrganoverseError(f"Failed to find files: {str(e)}")