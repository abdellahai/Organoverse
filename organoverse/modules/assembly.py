"""
Assembly Module for Organoverse workbench.

This module performs k-mer based read extraction and multi-reference assembly
using existing bioinformatics tools enhanced with AI guidance.
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import shutil

from .base_module import BaseModule
from ..utils.sequence_utils import parse_fasta, write_fasta, count_kmers
from ..utils.file_utils import create_temp_dir, cleanup_temp_files
from ..utils.exceptions import ModuleError, AssemblyError


class AssemblyModule(BaseModule):
    """
    K-mer based read extraction and multi-reference assembly module.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Assembly module.
        
        Args:
            config: Module configuration
        """
        super().__init__(config, "assembly")
        
        # Assembly parameters
        self.kmer_size = config.get('kmer_size', 31)
        self.coverage_cutoff = config.get('coverage_cutoff', 10.0)
        self.min_contig_length = config.get('min_contig_length', 500)
        self.assembler = config.get('assembler', 'spades')
        
        # Validate assembler availability
        self._check_assembler_availability()
    
    def run(self, reads_1: str, reads_2: Optional[str], references: List[Dict[str, Any]], 
            organelle: str) -> Dict[str, Any]:
        """
        Run assembly pipeline.
        
        Args:
            reads_1: Forward reads file
            reads_2: Reverse reads file (optional)
            references: List of reference genomes
            organelle: Target organelle type
        
        Returns:
            Assembly results
        """
        self.logger.info(f"Starting assembly for {organelle}")
        
        try:
            results = {
                'organelle': organelle,
                'assembler': self.assembler,
                'kmer_size': self.kmer_size,
                'coverage_cutoff': self.coverage_cutoff
            }
            
            # Step 1: Extract organellar reads
            self.logger.info("Extracting organellar reads")
            extracted_reads = self._extract_organellar_reads(reads_1, reads_2, references)
            results['extracted_reads'] = extracted_reads
            
            # Step 2: Run assembly
            self.logger.info(f"Running {self.assembler} assembly")
            assembly_result = self._run_assembly(extracted_reads, references)
            results['assembly_result'] = assembly_result
            
            # Step 3: Filter and validate contigs
            self.logger.info("Filtering and validating contigs")
            filtered_contigs = self._filter_contigs(assembly_result['contigs'])
            results['filtered_contigs'] = filtered_contigs
            
            # Step 4: Select best assembly
            self.logger.info("Selecting best assembly")
            best_assembly = self._select_best_assembly(filtered_contigs, references)
            results['best_assembly'] = best_assembly
            
            # Step 5: Save final assembly
            final_assembly_path = self.module_output_dir / f"{organelle}_assembly.fasta"
            write_fasta(best_assembly['contigs'], str(final_assembly_path))
            results['assembly_path'] = str(final_assembly_path)
            
            # Save results
            self.save_results(results)
            
            self.logger.info(f"Assembly completed. Final assembly: {final_assembly_path}")
            
            return results
            
        except Exception as e:
            self.logger.exception("Assembly failed")
            raise ModuleError(f"Assembly failed: {str(e)}")
    
    def _check_assembler_availability(self):
        """Check if the selected assembler is available."""
        assembler_commands = {
            'spades': 'spades.py',
            'novoplasty': 'NOVOPlasty.pl',
            'getorganelle': 'get_organelle_from_reads.py'
        }
        
        command = assembler_commands.get(self.assembler)
        if not command:
            raise ModuleError(f"Unknown assembler: {self.assembler}")
        
        try:
            result = subprocess.run(['which', command], capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.warning(f"Assembler {command} not found in PATH. Assembly may fail.")
        except Exception as e:
            self.logger.warning(f"Could not check assembler availability: {str(e)}")
    
    def _extract_organellar_reads(self, reads_1: str, reads_2: Optional[str], 
                                references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Extract organellar reads using k-mer matching."""
        try:
            # Build reference k-mer set
            reference_kmers = self._build_reference_kmers(references)
            
            # Extract matching reads
            extracted_1 = self.temp_dir / "extracted_R1.fastq"
            extracted_2 = self.temp_dir / "extracted_R2.fastq" if reads_2 else None
            
            stats = self._extract_reads_by_kmers(
                reads_1, reads_2, reference_kmers, str(extracted_1), 
                str(extracted_2) if extracted_2 else None
            )
            
            return {
                'extracted_reads_1': str(extracted_1),
                'extracted_reads_2': str(extracted_2) if extracted_2 else None,
                'extraction_stats': stats,
                'reference_kmers_count': len(reference_kmers)
            }
            
        except Exception as e:
            raise AssemblyError(f"Read extraction failed: {str(e)}")
    
    def _build_reference_kmers(self, references: List[Dict[str, Any]]) -> set:
        """Build k-mer set from reference genomes."""
        all_kmers = set()
        
        for ref in references:
            try:
                sequences = parse_fasta(ref['sequence_file'])
                for seq_record in sequences:
                    kmers = count_kmers(seq_record.seq, self.kmer_size)
                    all_kmers.update(kmers.keys())
            except Exception as e:
                self.logger.warning(f"Failed to process reference {ref['accession']}: {str(e)}")
                continue
        
        self.logger.info(f"Built k-mer set with {len(all_kmers)} unique k-mers")
        return all_kmers
    
    def _extract_reads_by_kmers(self, reads_1: str, reads_2: Optional[str], 
                              reference_kmers: set, output_1: str, 
                              output_2: Optional[str]) -> Dict[str, int]:
        """Extract reads that match reference k-mers."""
        from Bio import SeqIO
        import gzip
        
        stats = {
            'total_reads': 0,
            'extracted_reads': 0,
            'extraction_rate': 0.0
        }
        
        # Open input files
        if reads_1.endswith('.gz'):
            handle_1 = gzip.open(reads_1, 'rt')
        else:
            handle_1 = open(reads_1, 'r')
        
        handle_2 = None
        if reads_2:
            if reads_2.endswith('.gz'):
                handle_2 = gzip.open(reads_2, 'rt')
            else:
                handle_2 = open(reads_2, 'r')
        
        # Open output files
        out_1 = open(output_1, 'w')
        out_2 = open(output_2, 'w') if output_2 else None
        
        try:
            if reads_2:
                # Paired-end reads
                for record_1, record_2 in zip(SeqIO.parse(handle_1, "fastq"), 
                                            SeqIO.parse(handle_2, "fastq")):
                    stats['total_reads'] += 1
                    
                    # Check if either read matches reference k-mers
                    if self._read_matches_kmers(record_1.seq, reference_kmers) or \
                       self._read_matches_kmers(record_2.seq, reference_kmers):
                        SeqIO.write(record_1, out_1, "fastq")
                        SeqIO.write(record_2, out_2, "fastq")
                        stats['extracted_reads'] += 1
            else:
                # Single-end reads
                for record in SeqIO.parse(handle_1, "fastq"):
                    stats['total_reads'] += 1
                    
                    if self._read_matches_kmers(record.seq, reference_kmers):
                        SeqIO.write(record, out_1, "fastq")
                        stats['extracted_reads'] += 1
        
        finally:
            # Close all files
            handle_1.close()
            if handle_2:
                handle_2.close()
            out_1.close()
            if out_2:
                out_2.close()
        
        # Calculate extraction rate
        if stats['total_reads'] > 0:
            stats['extraction_rate'] = stats['extracted_reads'] / stats['total_reads']
        
        self.logger.info(f"Extracted {stats['extracted_reads']}/{stats['total_reads']} reads "
                        f"({stats['extraction_rate']:.2%})")
        
        return stats
    
    def _read_matches_kmers(self, sequence: str, reference_kmers: set, 
                          min_matches: int = 3) -> bool:
        """Check if a read matches reference k-mers."""
        seq_str = str(sequence).upper()
        matches = 0
        
        for i in range(len(seq_str) - self.kmer_size + 1):
            kmer = seq_str[i:i + self.kmer_size]
            if 'N' not in kmer and kmer in reference_kmers:
                matches += 1
                if matches >= min_matches:
                    return True
        
        return False
    
    def _run_assembly(self, extracted_reads: Dict[str, Any], 
                     references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Run genome assembly using the selected assembler."""
        if self.assembler == 'spades':
            return self._run_spades_assembly(extracted_reads, references)
        elif self.assembler == 'novoplasty':
            return self._run_novoplasty_assembly(extracted_reads, references)
        elif self.assembler == 'getorganelle':
            return self._run_getorganelle_assembly(extracted_reads, references)
        else:
            raise AssemblyError(f"Unsupported assembler: {self.assembler}")
    
    def _run_spades_assembly(self, extracted_reads: Dict[str, Any], 
                           references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Run SPAdes assembly."""
        try:
            output_dir = self.temp_dir / "spades_output"
            
            # Build SPAdes command
            cmd = [
                'spades.py',
                '--meta',  # Use metagenomic mode
                '-k', str(self.kmer_size),
                '--cov-cutoff', str(self.coverage_cutoff),
                '-t', str(self.threads),
                '-m', str(self.memory),
                '-o', str(output_dir)
            ]
            
            # Add input reads
            cmd.extend(['-1', extracted_reads['extracted_reads_1']])
            if extracted_reads['extracted_reads_2']:
                cmd.extend(['-2', extracted_reads['extracted_reads_2']])
            else:
                cmd.extend(['-s', extracted_reads['extracted_reads_1']])
            
            # Run SPAdes
            self.logger.debug(f"Running command: {' '.join(cmd)}")
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.temp_dir)
            
            if result.returncode != 0:
                raise AssemblyError(f"SPAdes failed: {result.stderr}")
            
            # Parse results
            contigs_file = output_dir / "contigs.fasta"
            if not contigs_file.exists():
                raise AssemblyError("SPAdes did not produce contigs.fasta")
            
            contigs = parse_fasta(str(contigs_file))
            
            return {
                'assembler': 'spades',
                'contigs_file': str(contigs_file),
                'contigs': contigs,
                'num_contigs': len(contigs),
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
        except Exception as e:
            raise AssemblyError(f"SPAdes assembly failed: {str(e)}")
    
    def _run_novoplasty_assembly(self, extracted_reads: Dict[str, Any], 
                                references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Run NOVOPlasty assembly (simplified implementation)."""
        # This is a placeholder - full NOVOPlasty integration would require
        # creating configuration files and handling its specific requirements
        self.logger.warning("NOVOPlasty assembly not fully implemented, falling back to SPAdes")
        return self._run_spades_assembly(extracted_reads, references)
    
    def _run_getorganelle_assembly(self, extracted_reads: Dict[str, Any], 
                                  references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Run GetOrganelle assembly (simplified implementation)."""
        # This is a placeholder - full GetOrganelle integration would require
        # proper configuration and handling of its specific parameters
        self.logger.warning("GetOrganelle assembly not fully implemented, falling back to SPAdes")
        return self._run_spades_assembly(extracted_reads, references)
    
    def _filter_contigs(self, contigs: List) -> List:
        """Filter contigs by length and quality."""
        filtered = []
        
        for contig in contigs:
            if len(contig.seq) >= self.min_contig_length:
                # Additional quality filters could be added here
                filtered.append(contig)
        
        self.logger.info(f"Filtered {len(contigs)} contigs to {len(filtered)} contigs")
        return filtered
    
    def _select_best_assembly(self, contigs: List, references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Select the best assembly based on various criteria."""
        if not contigs:
            raise AssemblyError("No contigs available for selection")
        
        # Sort contigs by length (descending)
        sorted_contigs = sorted(contigs, key=lambda x: len(x.seq), reverse=True)
        
        # Calculate assembly statistics
        total_length = sum(len(contig.seq) for contig in sorted_contigs)
        n50 = self._calculate_n50([len(contig.seq) for contig in sorted_contigs])
        
        return {
            'contigs': sorted_contigs,
            'num_contigs': len(sorted_contigs),
            'total_length': total_length,
            'n50': n50,
            'longest_contig': len(sorted_contigs[0].seq) if sorted_contigs else 0
        }
    
    def _calculate_n50(self, lengths: List[int]) -> int:
        """Calculate N50 statistic."""
        if not lengths:
            return 0
        
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target_length = total_length / 2
        
        cumulative_length = 0
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= target_length:
                return length
        
        return 0