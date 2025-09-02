"""
Sequence processing utilities for Organoverse workbench.
"""

import os
import gzip
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Iterator, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

from .exceptions import OrganoverseError


def parse_fasta(fasta_path: str) -> List[SeqRecord]:
    """
    Parse FASTA file and return list of SeqRecord objects.
    
    Args:
        fasta_path: Path to FASTA file
    
    Returns:
        List of SeqRecord objects
    """
    try:
        records = []
        
        if fasta_path.endswith('.gz'):
            with gzip.open(fasta_path, 'rt') as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        else:
            with open(fasta_path, 'r') as handle:
                records = list(SeqIO.parse(handle, "fasta"))
        
        return records
        
    except Exception as e:
        raise OrganoverseError(f"Failed to parse FASTA file {fasta_path}: {str(e)}")


def write_fasta(sequences: List[SeqRecord], output_path: str, compress: bool = False):
    """
    Write sequences to FASTA file.
    
    Args:
        sequences: List of SeqRecord objects
        output_path: Output file path
        compress: Whether to compress output
    """
    try:
        # Ensure output directory exists
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        
        if compress or output_path.endswith('.gz'):
            with gzip.open(output_path, 'wt') as handle:
                SeqIO.write(sequences, handle, "fasta")
        else:
            with open(output_path, 'w') as handle:
                SeqIO.write(sequences, handle, "fasta")
                
    except Exception as e:
        raise OrganoverseError(f"Failed to write FASTA file {output_path}: {str(e)}")


def parse_fastq(fastq_path: str, max_records: Optional[int] = None) -> List[SeqRecord]:
    """
    Parse FASTQ file and return list of SeqRecord objects.
    
    Args:
        fastq_path: Path to FASTQ file
        max_records: Maximum number of records to read
    
    Returns:
        List of SeqRecord objects
    """
    try:
        records = []
        
        if fastq_path.endswith('.gz'):
            handle = gzip.open(fastq_path, 'rt')
        else:
            handle = open(fastq_path, 'r')
        
        with handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                records.append(record)
                if max_records and i + 1 >= max_records:
                    break
        
        return records
        
    except Exception as e:
        raise OrganoverseError(f"Failed to parse FASTQ file {fastq_path}: {str(e)}")


def calculate_n50(sequences: List[Union[str, SeqRecord]]) -> int:
    """
    Calculate N50 statistic for a set of sequences.
    
    Args:
        sequences: List of sequences (strings or SeqRecord objects)
    
    Returns:
        N50 value
    """
    try:
        # Extract sequence lengths
        if isinstance(sequences[0], SeqRecord):
            lengths = [len(seq.seq) for seq in sequences]
        else:
            lengths = [len(seq) for seq in sequences]
        
        # Sort lengths in descending order
        lengths.sort(reverse=True)
        
        # Calculate total length
        total_length = sum(lengths)
        target_length = total_length / 2
        
        # Find N50
        cumulative_length = 0
        for length in lengths:
            cumulative_length += length
            if cumulative_length >= target_length:
                return length
        
        return 0
        
    except Exception as e:
        raise OrganoverseError(f"Failed to calculate N50: {str(e)}")


def calculate_gc_content(sequence: Union[str, Seq]) -> float:
    """
    Calculate GC content of a sequence.
    
    Args:
        sequence: DNA sequence
    
    Returns:
        GC content as fraction (0.0 to 1.0)
    """
    try:
        seq_str = str(sequence).upper()
        gc_count = seq_str.count('G') + seq_str.count('C')
        total_count = len(seq_str)
        
        if total_count == 0:
            return 0.0
        
        return gc_count / total_count
        
    except Exception as e:
        raise OrganoverseError(f"Failed to calculate GC content: {str(e)}")


def reverse_complement(sequence: Union[str, Seq]) -> str:
    """
    Get reverse complement of DNA sequence.
    
    Args:
        sequence: DNA sequence
    
    Returns:
        Reverse complement sequence
    """
    try:
        if isinstance(sequence, str):
            seq_obj = Seq(sequence)
        else:
            seq_obj = sequence
        
        return str(seq_obj.reverse_complement())
        
    except Exception as e:
        raise OrganoverseError(f"Failed to get reverse complement: {str(e)}")


def translate_sequence(sequence: Union[str, Seq], table: int = 1) -> str:
    """
    Translate DNA sequence to protein.
    
    Args:
        sequence: DNA sequence
        table: Genetic code table number
    
    Returns:
        Protein sequence
    """
    try:
        if isinstance(sequence, str):
            seq_obj = Seq(sequence)
        else:
            seq_obj = sequence
        
        return str(seq_obj.translate(table=table))
        
    except Exception as e:
        raise OrganoverseError(f"Failed to translate sequence: {str(e)}")


def find_orfs(sequence: Union[str, Seq], min_length: int = 100) -> List[Tuple[int, int, int, str]]:
    """
    Find open reading frames in a sequence.
    
    Args:
        sequence: DNA sequence
        min_length: Minimum ORF length in nucleotides
    
    Returns:
        List of tuples (start, end, frame, protein_sequence)
    """
    try:
        seq_str = str(sequence).upper()
        orfs = []
        
        # Check all 6 reading frames
        for strand in [1, -1]:
            if strand == -1:
                seq_str = reverse_complement(seq_str)
            
            for frame in range(3):
                frame_seq = seq_str[frame:]
                
                # Find start and stop codons
                start_codons = ['ATG']
                stop_codons = ['TAA', 'TAG', 'TGA']
                
                i = 0
                while i < len(frame_seq) - 2:
                    codon = frame_seq[i:i+3]
                    
                    if codon in start_codons:
                        # Found start codon, look for stop codon
                        j = i + 3
                        while j < len(frame_seq) - 2:
                            stop_codon = frame_seq[j:j+3]
                            if stop_codon in stop_codons:
                                orf_length = j + 3 - i
                                if orf_length >= min_length:
                                    orf_seq = frame_seq[i:j+3]
                                    protein = translate_sequence(orf_seq)
                                    
                                    if strand == 1:
                                        start_pos = frame + i
                                        end_pos = frame + j + 3
                                        frame_num = frame + 1
                                    else:
                                        start_pos = len(sequence) - (frame + j + 3)
                                        end_pos = len(sequence) - (frame + i)
                                        frame_num = -(frame + 1)
                                    
                                    orfs.append((start_pos, end_pos, frame_num, protein))
                                break
                            j += 3
                    i += 3
        
        return orfs
        
    except Exception as e:
        raise OrganoverseError(f"Failed to find ORFs: {str(e)}")


def count_kmers(sequence: Union[str, Seq], k: int) -> Dict[str, int]:
    """
    Count k-mers in a sequence.
    
    Args:
        sequence: DNA sequence
        k: K-mer size
    
    Returns:
        Dictionary of k-mer counts
    """
    try:
        seq_str = str(sequence).upper()
        kmer_counts = {}
        
        for i in range(len(seq_str) - k + 1):
            kmer = seq_str[i:i+k]
            if 'N' not in kmer:  # Skip k-mers with ambiguous bases
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
        
        return kmer_counts
        
    except Exception as e:
        raise OrganoverseError(f"Failed to count k-mers: {str(e)}")


def calculate_sequence_complexity(sequence: Union[str, Seq], window_size: int = 50) -> float:
    """
    Calculate sequence complexity using Shannon entropy.
    
    Args:
        sequence: DNA sequence
        window_size: Window size for complexity calculation
    
    Returns:
        Average complexity score (0.0 to 2.0)
    """
    try:
        seq_str = str(sequence).upper()
        complexities = []
        
        for i in range(0, len(seq_str) - window_size + 1, window_size):
            window = seq_str[i:i+window_size]
            
            # Count nucleotides
            counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
            for base in window:
                if base in counts:
                    counts[base] += 1
            
            # Calculate Shannon entropy
            total = sum(counts.values())
            if total == 0:
                complexity = 0.0
            else:
                entropy = 0.0
                for count in counts.values():
                    if count > 0:
                        p = count / total
                        entropy -= p * np.log2(p)
                complexity = entropy
            
            complexities.append(complexity)
        
        return np.mean(complexities) if complexities else 0.0
        
    except Exception as e:
        raise OrganoverseError(f"Failed to calculate sequence complexity: {str(e)}")


def mask_low_complexity_regions(sequence: Union[str, Seq], threshold: float = 1.0, 
                               window_size: int = 50, mask_char: str = 'N') -> str:
    """
    Mask low complexity regions in a sequence.
    
    Args:
        sequence: DNA sequence
        threshold: Complexity threshold below which to mask
        window_size: Window size for complexity calculation
        mask_char: Character to use for masking
    
    Returns:
        Masked sequence
    """
    try:
        seq_str = str(sequence).upper()
        masked_seq = list(seq_str)
        
        for i in range(0, len(seq_str) - window_size + 1, window_size):
            window = seq_str[i:i+window_size]
            complexity = calculate_sequence_complexity(window, window_size)
            
            if complexity < threshold:
                for j in range(i, min(i + window_size, len(seq_str))):
                    masked_seq[j] = mask_char
        
        return ''.join(masked_seq)
        
    except Exception as e:
        raise OrganoverseError(f"Failed to mask low complexity regions: {str(e)}")


def split_sequences_by_length(sequences: List[SeqRecord], min_length: int, 
                             max_length: int) -> Tuple[List[SeqRecord], List[SeqRecord]]:
    """
    Split sequences into two groups based on length.
    
    Args:
        sequences: List of SeqRecord objects
        min_length: Minimum length threshold
        max_length: Maximum length threshold
    
    Returns:
        Tuple of (sequences_in_range, sequences_out_of_range)
    """
    try:
        in_range = []
        out_of_range = []
        
        for seq in sequences:
            seq_len = len(seq.seq)
            if min_length <= seq_len <= max_length:
                in_range.append(seq)
            else:
                out_of_range.append(seq)
        
        return in_range, out_of_range
        
    except Exception as e:
        raise OrganoverseError(f"Failed to split sequences by length: {str(e)}")