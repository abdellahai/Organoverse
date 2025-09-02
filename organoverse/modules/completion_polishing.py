"""
Completion and Polishing Module for Organoverse workbench.

This module performs AI-powered assembly completion and polishing using
LSTM gap filling and iterative improvement algorithms.
"""

import os
import re
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

from .base_module import BaseModule
from ..models.lstm_gap_filler import LSTMGapFiller
from ..utils.sequence_utils import parse_fasta, write_fasta, calculate_n50
from ..utils.exceptions import ModuleError


class CompletionPolishingModule(BaseModule):
    """
    AI-powered assembly completion and polishing module.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Completion and Polishing module.
        
        Args:
            config: Module configuration
        """
        super().__init__(config, "completion_polishing")
        
        # Model configuration
        self.lstm_model_path = config.get('lstm_model_path')
        self.device = config.get('device', 'cpu')
        self.polish_rounds = config.get('polish_rounds', 3)
        self.max_iterations = config.get('max_iterations', 10)
        
        # Gap filling parameters
        self.min_gap_size = config.get('min_gap_size', 10)
        self.max_gap_size = config.get('max_gap_size', 1000)
        self.context_length = config.get('context_length', 100)
        
        # Initialize LSTM model if available
        self.lstm_model = None
        if self.lstm_model_path and os.path.exists(self.lstm_model_path):
            try:
                self.lstm_model = LSTMGapFiller(self.lstm_model_path, device=self.device)
                self.logger.info("LSTM gap filler loaded successfully")
            except Exception as e:
                self.logger.warning(f"Failed to load LSTM model: {str(e)}")
    
    def run(self, assembly: str, reads_1: str, reads_2: Optional[str], 
            references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Run completion and polishing pipeline.
        
        Args:
            assembly: Path to initial assembly
            reads_1: Forward reads file
            reads_2: Reverse reads file (optional)
            references: List of reference genomes
        
        Returns:
            Completion and polishing results
        """
        self.logger.info("Starting assembly completion and polishing")
        
        try:
            results = {
                'initial_assembly': assembly,
                'polish_rounds': self.polish_rounds,
                'max_iterations': self.max_iterations
            }
            
            # Load initial assembly
            current_assembly = parse_fasta(assembly)
            if not current_assembly:
                raise ModuleError("No contigs found in assembly")
            
            results['initial_stats'] = self._calculate_assembly_stats(current_assembly)
            
            # Step 1: Gap identification and filling
            if self.lstm_model:
                self.logger.info("Identifying and filling gaps with LSTM")
                gap_filled_assembly = self._fill_gaps_with_lstm(current_assembly, references)
                results['gap_filling'] = self._calculate_assembly_stats(gap_filled_assembly)
                current_assembly = gap_filled_assembly
            
            # Step 2: Iterative polishing
            self.logger.info("Starting iterative polishing")
            polished_assembly = self._iterative_polishing(current_assembly, reads_1, reads_2, references)
            results['polishing'] = self._calculate_assembly_stats(polished_assembly)
            
            # Step 3: Final validation
            self.logger.info("Performing final validation")
            validation_results = self._validate_final_assembly(polished_assembly, references)
            results['validation'] = validation_results
            
            # Save final assembly
            final_assembly_path = self.module_output_dir / "polished_assembly.fasta"
            write_fasta(polished_assembly, str(final_assembly_path))
            results['polished_assembly'] = str(final_assembly_path)
            
            # Save results
            self.save_results(results)
            
            self.logger.info(f"Completion and polishing completed. Final assembly: {final_assembly_path}")
            
            return results
            
        except Exception as e:
            self.logger.exception("Completion and polishing failed")
            raise ModuleError(f"Completion and polishing failed: {str(e)}")
    
    def _calculate_assembly_stats(self, contigs: List) -> Dict[str, Any]:
        """Calculate assembly statistics."""
        if not contigs:
            return {'num_contigs': 0, 'total_length': 0, 'n50': 0, 'longest_contig': 0}
        
        lengths = [len(contig.seq) for contig in contigs]
        
        return {
            'num_contigs': len(contigs),
            'total_length': sum(lengths),
            'n50': calculate_n50([str(contig.seq) for contig in contigs]),
            'longest_contig': max(lengths),
            'shortest_contig': min(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths)
        }
    
    def _fill_gaps_with_lstm(self, contigs: List, references: List[Dict[str, Any]]) -> List:
        """Fill gaps in assembly using LSTM model."""
        try:
            filled_contigs = []
            
            for contig in contigs:
                # Identify gaps (represented as N's)
                gaps = self._identify_gaps(str(contig.seq))
                
                if not gaps:
                    filled_contigs.append(contig)
                    continue
                
                # Fill gaps using LSTM
                filled_sequence = self._fill_sequence_gaps(str(contig.seq), gaps)
                
                # Create new contig with filled sequence
                from Bio.SeqRecord import SeqRecord
                from Bio.Seq import Seq
                
                filled_contig = SeqRecord(
                    Seq(filled_sequence),
                    id=contig.id + "_gap_filled",
                    description=contig.description + " [gap filled]"
                )
                filled_contigs.append(filled_contig)
            
            return filled_contigs
            
        except Exception as e:
            self.logger.warning(f"Gap filling failed: {str(e)}")
            return contigs
    
    def _identify_gaps(self, sequence: str) -> List[Tuple[int, int]]:
        """Identify gaps (N regions) in sequence."""
        gaps = []
        
        # Find all N regions
        for match in re.finditer(r'N+', sequence.upper()):
            start, end = match.span()
            gap_size = end - start
            
            # Only consider gaps within size limits
            if self.min_gap_size <= gap_size <= self.max_gap_size:
                gaps.append((start, end))
        
        return gaps
    
    def _fill_sequence_gaps(self, sequence: str, gaps: List[Tuple[int, int]]) -> str:
        """Fill gaps in a sequence using LSTM model."""
        filled_sequence = sequence
        
        # Process gaps from right to left to maintain indices
        for start, end in reversed(gaps):
            gap_size = end - start
            
            # Extract context
            left_start = max(0, start - self.context_length)
            left_context = sequence[left_start:start]
            
            right_end = min(len(sequence), end + self.context_length)
            right_context = sequence[end:right_end]
            
            # Predict gap sequence
            try:
                predicted_sequences = self.lstm_model.predict_gap_sequences(
                    [left_context], [right_context], gap_size
                )
                
                if predicted_sequences:
                    # Score predictions and select best
                    scores = self.lstm_model.score_gap_predictions(predicted_sequences)
                    best_idx = np.argmax(scores)
                    gap_fill = predicted_sequences[best_idx]
                    
                    # Replace gap with prediction
                    filled_sequence = filled_sequence[:start] + gap_fill + filled_sequence[end:]
                    
                    self.logger.debug(f"Filled gap at {start}-{end} with {len(gap_fill)} bp sequence")
            
            except Exception as e:
                self.logger.warning(f"Failed to fill gap at {start}-{end}: {str(e)}")
                continue
        
        return filled_sequence
    
    def _iterative_polishing(self, contigs: List, reads_1: str, reads_2: Optional[str], 
                           references: List[Dict[str, Any]]) -> List:
        """Perform iterative polishing of assembly."""
        current_contigs = contigs
        
        for round_num in range(self.polish_rounds):
            self.logger.info(f"Polishing round {round_num + 1}/{self.polish_rounds}")
            
            try:
                # Polish contigs
                polished_contigs = self._polish_contigs(current_contigs, references)
                
                # Check for improvement
                current_stats = self._calculate_assembly_stats(current_contigs)
                polished_stats = self._calculate_assembly_stats(polished_contigs)
                
                improvement = polished_stats['n50'] - current_stats['n50']
                
                if improvement > 0:
                    self.logger.info(f"Round {round_num + 1}: N50 improved by {improvement} bp")
                    current_contigs = polished_contigs
                else:
                    self.logger.info(f"Round {round_num + 1}: No improvement, stopping")
                    break
            
            except Exception as e:
                self.logger.warning(f"Polishing round {round_num + 1} failed: {str(e)}")
                break
        
        return current_contigs
    
    def _polish_contigs(self, contigs: List, references: List[Dict[str, Any]]) -> List:
        """Polish contigs using reference-based correction."""
        polished_contigs = []
        
        for contig in contigs:
            try:
                # Simple polishing: correct obvious errors
                polished_seq = self._correct_sequence_errors(str(contig.seq), references)
                
                # Create polished contig
                from Bio.SeqRecord import SeqRecord
                from Bio.Seq import Seq
                
                polished_contig = SeqRecord(
                    Seq(polished_seq),
                    id=contig.id + "_polished",
                    description=contig.description + " [polished]"
                )
                polished_contigs.append(polished_contig)
            
            except Exception as e:
                self.logger.warning(f"Failed to polish contig {contig.id}: {str(e)}")
                polished_contigs.append(contig)
        
        return polished_contigs
    
    def _correct_sequence_errors(self, sequence: str, references: List[Dict[str, Any]]) -> str:
        """Correct obvious sequence errors using references."""
        corrected_sequence = sequence
        
        # This is a simplified implementation
        # In practice, this would involve more sophisticated error correction
        
        # Remove ambiguous bases
        corrected_sequence = re.sub(r'[^ATGC]', 'N', corrected_sequence.upper())
        
        # Additional corrections could be implemented here based on reference alignment
        
        return corrected_sequence
    
    def _validate_final_assembly(self, contigs: List, references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Validate the final assembly."""
        validation = {
            'num_contigs': len(contigs),
            'total_length': sum(len(contig.seq) for contig in contigs),
            'quality_checks': []
        }
        
        # Check 1: Reasonable genome size
        expected_sizes = {
            'chloroplast': (120000, 200000),
            'mitochondrion': (200000, 2000000)
        }
        
        for organelle, (min_size, max_size) in expected_sizes.items():
            if min_size <= validation['total_length'] <= max_size:
                validation['quality_checks'].append(f"Size appropriate for {organelle}")
                validation['likely_organelle'] = organelle
                break
        
        # Check 2: Contiguity
        if validation['num_contigs'] == 1:
            validation['quality_checks'].append("Single contig assembly (excellent contiguity)")
        elif validation['num_contigs'] <= 5:
            validation['quality_checks'].append("Low contig count (good contiguity)")
        else:
            validation['quality_checks'].append("High contig count (fragmented assembly)")
        
        # Check 3: N content
        total_ns = sum(str(contig.seq).upper().count('N') for contig in contigs)
        n_percentage = (total_ns / validation['total_length']) * 100 if validation['total_length'] > 0 else 0
        
        if n_percentage < 1:
            validation['quality_checks'].append("Low N content (good quality)")
        elif n_percentage < 5:
            validation['quality_checks'].append("Moderate N content")
        else:
            validation['quality_checks'].append("High N content (poor quality)")
        
        validation['n_percentage'] = n_percentage
        
        return validation