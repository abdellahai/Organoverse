"""
Input validation utilities for Organoverse workbench.
"""

import os
import re
from pathlib import Path
from typing import Optional, List
from Bio import SeqIO

from .exceptions import ValidationError


def validate_fastq_files(reads_1: str, reads_2: Optional[str] = None) -> bool:
    """
    Validate FASTQ input files.
    
    Args:
        reads_1: Path to forward reads file
        reads_2: Optional path to reverse reads file
    
    Returns:
        True if files are valid
    
    Raises:
        ValidationError: If files are invalid
    """
    # Check if files exist
    if not os.path.exists(reads_1):
        raise ValidationError(f"Forward reads file not found: {reads_1}")
    
    if reads_2 and not os.path.exists(reads_2):
        raise ValidationError(f"Reverse reads file not found: {reads_2}")
    
    # Check file extensions
    valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    
    if not any(reads_1.endswith(ext) for ext in valid_extensions):
        raise ValidationError(f"Invalid file extension for forward reads: {reads_1}")
    
    if reads_2 and not any(reads_2.endswith(ext) for ext in valid_extensions):
        raise ValidationError(f"Invalid file extension for reverse reads: {reads_2}")
    
    # Check if files are readable and contain valid FASTQ data
    try:
        _validate_fastq_format(reads_1)
        if reads_2:
            _validate_fastq_format(reads_2)
    except Exception as e:
        raise ValidationError(f"Invalid FASTQ format: {str(e)}")
    
    return True


def _validate_fastq_format(fastq_path: str, max_records: int = 100) -> bool:
    """
    Validate FASTQ file format by reading a few records.
    
    Args:
        fastq_path: Path to FASTQ file
        max_records: Maximum number of records to check
    
    Returns:
        True if format is valid
    """
    try:
        if fastq_path.endswith('.gz'):
            import gzip
            handle = gzip.open(fastq_path, 'rt')
        else:
            handle = open(fastq_path, 'r')
        
        with handle:
            record_count = 0
            for record in SeqIO.parse(handle, "fastq"):
                record_count += 1
                if record_count >= max_records:
                    break
            
            if record_count == 0:
                raise ValidationError("No valid FASTQ records found")
        
        return True
        
    except Exception as e:
        raise ValidationError(f"FASTQ validation failed: {str(e)}")


def validate_species_name(species: str) -> bool:
    """
    Validate species name format.
    
    Args:
        species: Species name (e.g., "Arabidopsis thaliana")
    
    Returns:
        True if species name is valid
    
    Raises:
        ValidationError: If species name is invalid
    """
    if not species or not isinstance(species, str):
        raise ValidationError("Species name must be a non-empty string")
    
    species = species.strip()
    
    if len(species) < 3:
        raise ValidationError("Species name too short")
    
    # Check for basic binomial nomenclature pattern
    # Should contain at least genus and species (two words)
    words = species.split()
    if len(words) < 2:
        raise ValidationError("Species name should contain at least genus and species")
    
    # Check for valid characters (letters, spaces, hyphens, periods)
    if not re.match(r'^[A-Za-z\s\-\.]+$', species):
        raise ValidationError("Species name contains invalid characters")
    
    # Check that genus starts with capital letter
    if not words[0][0].isupper():
        raise ValidationError("Genus name should start with a capital letter")
    
    return True


def validate_file_path(file_path: str, must_exist: bool = True) -> bool:
    """
    Validate file path.
    
    Args:
        file_path: Path to validate
        must_exist: Whether file must exist
    
    Returns:
        True if path is valid
    
    Raises:
        ValidationError: If path is invalid
    """
    if not file_path:
        raise ValidationError("File path cannot be empty")
    
    path = Path(file_path)
    
    if must_exist and not path.exists():
        raise ValidationError(f"File does not exist: {file_path}")
    
    if must_exist and not path.is_file():
        raise ValidationError(f"Path is not a file: {file_path}")
    
    # Check if parent directory exists for output files
    if not must_exist and not path.parent.exists():
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            raise ValidationError(f"Cannot create parent directory: {str(e)}")
    
    return True


def validate_organelle_type(organelle: str) -> bool:
    """
    Validate organelle type.
    
    Args:
        organelle: Organelle type
    
    Returns:
        True if valid
    
    Raises:
        ValidationError: If invalid
    """
    valid_organelles = ['chloroplast', 'mitochondrion', 'both']
    
    if organelle not in valid_organelles:
        raise ValidationError(f"Invalid organelle type: {organelle}. Must be one of {valid_organelles}")
    
    return True


def validate_kmer_size(kmer_size: int) -> bool:
    """
    Validate k-mer size parameter.
    
    Args:
        kmer_size: K-mer size
    
    Returns:
        True if valid
    
    Raises:
        ValidationError: If invalid
    """
    if not isinstance(kmer_size, int):
        raise ValidationError("K-mer size must be an integer")
    
    if kmer_size < 15 or kmer_size > 127:
        raise ValidationError("K-mer size must be between 15 and 127")
    
    if kmer_size % 2 == 0:
        raise ValidationError("K-mer size must be odd")
    
    return True


def validate_coverage_cutoff(coverage: float) -> bool:
    """
    Validate coverage cutoff parameter.
    
    Args:
        coverage: Coverage cutoff
    
    Returns:
        True if valid
    
    Raises:
        ValidationError: If invalid
    """
    if not isinstance(coverage, (int, float)):
        raise ValidationError("Coverage cutoff must be a number")
    
    if coverage <= 0:
        raise ValidationError("Coverage cutoff must be positive")
    
    if coverage > 1000:
        raise ValidationError("Coverage cutoff seems too high (>1000x)")
    
    return True


def validate_threads(threads: int) -> bool:
    """
    Validate number of threads.
    
    Args:
        threads: Number of threads
    
    Returns:
        True if valid
    
    Raises:
        ValidationError: If invalid
    """
    if not isinstance(threads, int):
        raise ValidationError("Number of threads must be an integer")
    
    if threads < 1:
        raise ValidationError("Number of threads must be at least 1")
    
    if threads > 128:
        raise ValidationError("Number of threads seems too high (>128)")
    
    return True


def validate_memory(memory: int) -> bool:
    """
    Validate memory limit.
    
    Args:
        memory: Memory limit in GB
    
    Returns:
        True if valid
    
    Raises:
        ValidationError: If invalid
    """
    if not isinstance(memory, int):
        raise ValidationError("Memory limit must be an integer")
    
    if memory < 1:
        raise ValidationError("Memory limit must be at least 1 GB")
    
    if memory > 1024:
        raise ValidationError("Memory limit seems too high (>1024 GB)")
    
    return True