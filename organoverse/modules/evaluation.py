"""
Evaluation Module for Organoverse workbench.

This module performs comprehensive assembly evaluation and comparison
with reference genomes.
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

from .base_module import BaseModule
from ..utils.sequence_utils import parse_fasta, calculate_n50, calculate_gc_content
from ..utils.exceptions import ModuleError


class EvaluationModule(BaseModule):
    """
    Comprehensive assembly evaluation and comparison module.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Evaluation module.
        
        Args:
            config: Module configuration
        """
        super().__init__(config, "evaluation")
        
        # Evaluation parameters
        self.min_alignment_length = config.get('min_alignment_length', 1000)
        self.identity_threshold = config.get('identity_threshold', 0.95)
    
    def run(self, assembly: str, references: List[Dict[str, Any]], 
            organelle: str) -> Dict[str, Any]:
        """
        Run comprehensive assembly evaluation.
        
        Args:
            assembly: Path to assembly file
            references: List of reference genomes
            organelle: Target organelle type
        
        Returns:
            Evaluation results
        """
        self.logger.info("Starting assembly evaluation")
        
        try:
            results = {
                'assembly_file': assembly,
                'organelle': organelle,
                'num_references': len(references)
            }
            
            # Load assembly
            contigs = parse_fasta(assembly)
            if not contigs:
                raise ModuleError("No contigs found in assembly")
            
            # Step 1: Basic assembly statistics
            self.logger.info("Calculating basic assembly statistics")
            basic_stats = self._calculate_basic_statistics(contigs)
            results['basic_statistics'] = basic_stats
            
            # Step 2: Sequence composition analysis
            self.logger.info("Analyzing sequence composition")
            composition_stats = self._analyze_sequence_composition(contigs)
            results['composition_analysis'] = composition_stats
            
            # Step 3: Reference comparison
            if references:
                self.logger.info("Comparing with reference genomes")
                reference_comparison = self._compare_with_references(contigs, references)
                results['reference_comparison'] = reference_comparison
            
            # Step 4: Quality assessment
            self.logger.info("Assessing assembly quality")
            quality_assessment = self._assess_assembly_quality(contigs, organelle)
            results['quality_assessment'] = quality_assessment
            
            # Step 5: Generate overall score
            overall_score = self._calculate_overall_score(results)
            results['overall_score'] = overall_score
            
            # Save results
            self.save_results(results)
            
            self.logger.info(f"Evaluation completed. Overall score: {overall_score:.3f}")
            
            return results
            
        except Exception as e:
            self.logger.exception("Assembly evaluation failed")
            raise ModuleError(f"Assembly evaluation failed: {str(e)}")
    
    def _calculate_basic_statistics(self, contigs: List) -> Dict[str, Any]:
        """Calculate basic assembly statistics."""
        if not contigs:
            return {}
        
        lengths = [len(contig.seq) for contig in contigs]
        
        stats = {
            'num_contigs': len(contigs),
            'total_length': sum(lengths),
            'longest_contig': max(lengths),
            'shortest_contig': min(lengths),
            'mean_length': np.mean(lengths),
            'median_length': np.median(lengths),
            'n50': calculate_n50([str(contig.seq) for contig in contigs]),
            'n90': self._calculate_nx(lengths, 0.9)
        }
        
        # Calculate L50 (number of contigs needed to cover 50% of assembly)
        stats['l50'] = self._calculate_lx(lengths, 0.5)
        stats['l90'] = self._calculate_lx(lengths, 0.9)
        
        return stats
    
    def _calculate_nx(self, lengths: List[int], x: float) -> int:
        """Calculate Nx statistic (e.g., N50, N90)."""
        if not lengths:
            return 0
        
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target_length = total_length * x
        
        cumulative_length = 0
        for length in sorted_lengths:
            cumulative_length += length
            if cumulative_length >= target_length:
                return length
        
        return 0
    
    def _calculate_lx(self, lengths: List[int], x: float) -> int:
        """Calculate Lx statistic (number of contigs for x% coverage)."""
        if not lengths:
            return 0
        
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        target_length = total_length * x
        
        cumulative_length = 0
        for i, length in enumerate(sorted_lengths):
            cumulative_length += length
            if cumulative_length >= target_length:
                return i + 1
        
        return len(lengths)
    
    def _analyze_sequence_composition(self, contigs: List) -> Dict[str, Any]:
        """Analyze sequence composition."""
        composition = {
            'gc_content': [],
            'n_content': [],
            'length_distribution': {}
        }
        
        all_sequence = ""
        for contig in contigs:
            seq_str = str(contig.seq).upper()
            all_sequence += seq_str
            
            # GC content per contig
            gc_content = calculate_gc_content(contig.seq)
            composition['gc_content'].append(gc_content)
            
            # N content per contig
            n_count = seq_str.count('N')
            n_content = n_count / len(seq_str) if len(seq_str) > 0 else 0
            composition['n_content'].append(n_content)
        
        # Overall composition
        composition['overall_gc_content'] = calculate_gc_content(all_sequence)
        composition['overall_n_content'] = all_sequence.count('N') / len(all_sequence) if all_sequence else 0
        
        # Statistics
        if composition['gc_content']:
            composition['gc_content_stats'] = {
                'mean': np.mean(composition['gc_content']),
                'std': np.std(composition['gc_content']),
                'min': np.min(composition['gc_content']),
                'max': np.max(composition['gc_content'])
            }
        
        if composition['n_content']:
            composition['n_content_stats'] = {
                'mean': np.mean(composition['n_content']),
                'std': np.std(composition['n_content']),
                'min': np.min(composition['n_content']),
                'max': np.max(composition['n_content'])
            }
        
        return composition
    
    def _compare_with_references(self, contigs: List, references: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Compare assembly with reference genomes."""
        comparison_results = {
            'reference_comparisons': [],
            'best_reference': None,
            'overall_similarity': 0.0
        }
        
        best_similarity = 0.0
        best_reference = None
        
        for ref in references:
            try:
                # Load reference sequence
                ref_sequences = parse_fasta(ref['sequence_file'])
                if not ref_sequences:
                    continue
                
                ref_sequence = ref_sequences[0]
                
                # Compare assembly with this reference
                similarity_metrics = self._calculate_similarity_metrics(contigs, ref_sequence)
                
                comparison = {
                    'reference_accession': ref['accession'],
                    'reference_organism': ref.get('organism', ''),
                    'reference_length': len(ref_sequence.seq),
                    'similarity_metrics': similarity_metrics
                }
                
                comparison_results['reference_comparisons'].append(comparison)
                
                # Track best reference
                overall_sim = similarity_metrics.get('overall_similarity', 0)
                if overall_sim > best_similarity:
                    best_similarity = overall_sim
                    best_reference = comparison
            
            except Exception as e:
                self.logger.warning(f"Failed to compare with reference {ref['accession']}: {str(e)}")
                continue
        
        comparison_results['best_reference'] = best_reference
        comparison_results['overall_similarity'] = best_similarity
        
        return comparison_results
    
    def _calculate_similarity_metrics(self, contigs: List, reference) -> Dict[str, float]:
        """Calculate similarity metrics between assembly and reference."""
        # This is a simplified implementation
        # In practice, you would use tools like BLAST, minimap2, or MUMmer
        
        assembly_length = sum(len(contig.seq) for contig in contigs)
        reference_length = len(reference.seq)
        
        # Length similarity
        length_ratio = min(assembly_length, reference_length) / max(assembly_length, reference_length)
        
        # GC content similarity
        assembly_gc = np.mean([calculate_gc_content(contig.seq) for contig in contigs])
        reference_gc = calculate_gc_content(reference.seq)
        gc_similarity = 1.0 - abs(assembly_gc - reference_gc)
        
        # Simple k-mer based similarity (simplified)
        kmer_similarity = self._calculate_kmer_similarity(contigs, reference)
        
        # Overall similarity (weighted average)
        overall_similarity = (length_ratio * 0.3 + gc_similarity * 0.2 + kmer_similarity * 0.5)
        
        return {
            'length_similarity': length_ratio,
            'gc_similarity': gc_similarity,
            'kmer_similarity': kmer_similarity,
            'overall_similarity': overall_similarity,
            'assembly_length': assembly_length,
            'reference_length': reference_length,
            'length_difference': abs(assembly_length - reference_length)
        }
    
    def _calculate_kmer_similarity(self, contigs: List, reference, k: int = 21) -> float:
        """Calculate k-mer based similarity (simplified implementation)."""
        try:
            from ..utils.sequence_utils import count_kmers
            
            # Get k-mers from assembly
            assembly_kmers = set()
            for contig in contigs:
                contig_kmers = count_kmers(contig.seq, k)
                assembly_kmers.update(contig_kmers.keys())
            
            # Get k-mers from reference
            reference_kmers = set(count_kmers(reference.seq, k).keys())
            
            # Calculate Jaccard similarity
            if not assembly_kmers and not reference_kmers:
                return 1.0
            
            intersection = len(assembly_kmers.intersection(reference_kmers))
            union = len(assembly_kmers.union(reference_kmers))
            
            return intersection / union if union > 0 else 0.0
            
        except Exception as e:
            self.logger.warning(f"K-mer similarity calculation failed: {str(e)}")
            return 0.0
    
    def _assess_assembly_quality(self, contigs: List, organelle: str) -> Dict[str, Any]:
        """Assess overall assembly quality."""
        quality = {
            'quality_score': 0.0,
            'quality_factors': [],
            'warnings': [],
            'recommendations': []
        }
        
        # Factor 1: Contiguity
        num_contigs = len(contigs)
        if num_contigs == 1:
            quality['quality_factors'].append(('Contiguity', 1.0, 'Single contig assembly'))
        elif num_contigs <= 3:
            quality['quality_factors'].append(('Contiguity', 0.8, 'Low contig count'))
        elif num_contigs <= 10:
            quality['quality_factors'].append(('Contiguity', 0.6, 'Moderate fragmentation'))
        else:
            quality['quality_factors'].append(('Contiguity', 0.3, 'High fragmentation'))
            quality['warnings'].append('Assembly is highly fragmented')
        
        # Factor 2: Size appropriateness
        total_length = sum(len(contig.seq) for contig in contigs)
        expected_sizes = {
            'chloroplast': (120000, 200000),
            'mitochondrion': (200000, 2000000)
        }
        
        if organelle in expected_sizes:
            min_size, max_size = expected_sizes[organelle]
            if min_size <= total_length <= max_size:
                quality['quality_factors'].append(('Size', 1.0, 'Appropriate size'))
            elif total_length < min_size:
                size_score = total_length / min_size
                quality['quality_factors'].append(('Size', size_score, 'Smaller than expected'))
                quality['warnings'].append('Assembly smaller than expected')
            else:
                size_score = max_size / total_length
                quality['quality_factors'].append(('Size', size_score, 'Larger than expected'))
                quality['warnings'].append('Assembly larger than expected')
        
        # Factor 3: N content
        total_ns = sum(str(contig.seq).upper().count('N') for contig in contigs)
        n_percentage = (total_ns / total_length) * 100 if total_length > 0 else 0
        
        if n_percentage < 1:
            quality['quality_factors'].append(('Completeness', 1.0, 'Low N content'))
        elif n_percentage < 5:
            quality['quality_factors'].append(('Completeness', 0.8, 'Moderate N content'))
        else:
            quality['quality_factors'].append(('Completeness', 0.5, 'High N content'))
            quality['warnings'].append('High N content indicates gaps')
        
        # Calculate overall quality score
        if quality['quality_factors']:
            scores = [factor[1] for factor in quality['quality_factors']]
            quality['quality_score'] = np.mean(scores)
        
        # Generate recommendations
        if quality['quality_score'] < 0.5:
            quality['recommendations'].append('Consider additional polishing or gap filling')
        if num_contigs > 5:
            quality['recommendations'].append('Consider scaffolding to reduce fragmentation')
        if n_percentage > 5:
            quality['recommendations'].append('Consider gap filling to reduce N content')
        
        return quality
    
    def _calculate_overall_score(self, results: Dict[str, Any]) -> float:
        """Calculate overall assembly score."""
        score_components = []
        
        # Basic statistics component
        basic_stats = results.get('basic_statistics', {})
        if basic_stats:
            # Contiguity score (fewer contigs is better)
            num_contigs = basic_stats.get('num_contigs', 1)
            contiguity_score = 1.0 / num_contigs if num_contigs > 0 else 0
            score_components.append(min(contiguity_score, 1.0) * 0.3)
        
        # Quality assessment component
        quality_assessment = results.get('quality_assessment', {})
        if quality_assessment:
            quality_score = quality_assessment.get('quality_score', 0.5)
            score_components.append(quality_score * 0.4)
        
        # Reference comparison component
        reference_comparison = results.get('reference_comparison', {})
        if reference_comparison:
            similarity_score = reference_comparison.get('overall_similarity', 0.5)
            score_components.append(similarity_score * 0.3)
        
        # Calculate weighted average
        if score_components:
            overall_score = sum(score_components)
        else:
            overall_score = 0.5  # Default neutral score
        
        return min(max(overall_score, 0.0), 1.0)  # Clamp to [0, 1]