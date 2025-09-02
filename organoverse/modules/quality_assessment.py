"""
Quality Assessment Module for Organoverse workbench.

This module performs AI-powered quality assessment of sequencing reads,
including contamination detection and k-mer analysis.
"""

import os
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
import pandas as pd
from Bio import SeqIO

from .base_module import BaseModule
from ..utils.sequence_utils import parse_fastq, count_kmers, calculate_gc_content
from ..utils.exceptions import ModuleError


class QualityAssessmentModule(BaseModule):
    """
    AI-powered quality assessment module using CNN and k-mer analysis.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Quality Assessment module.
        
        Args:
            config: Module configuration
        """
        super().__init__(config, "quality_assessment")
        
        # Model configuration
        self.model_path = config.get('model_path')
        self.device = config.get('device', 'cpu')
        self.batch_size = config.get('batch_size', 32)
        self.confidence_threshold = config.get('confidence_threshold', 0.8)
        
        # K-mer analysis parameters
        self.kmer_size = config.get('kmer_size', 21)
        self.sample_size = config.get('sample_size', 10000)  # Number of reads to sample
        
        # Initialize CNN model if available
        self.cnn_model = None
        if self.model_path and os.path.exists(self.model_path):
            try:
                from ..models.cnn_quality_classifier import CNNQualityClassifier
                self.cnn_model = CNNQualityClassifier(self.model_path, device=self.device)
                self.logger.info("CNN quality classifier loaded successfully")
            except Exception as e:
                self.logger.warning(f"Failed to load CNN model: {str(e)}")
    
    def run(self, reads_1: str, reads_2: Optional[str] = None, species: str = "") -> Dict[str, Any]:
        """
        Run quality assessment on sequencing reads.
        
        Args:
            reads_1: Path to forward reads
            reads_2: Path to reverse reads (optional)
            species: Species name for context
        
        Returns:
            Quality assessment results
        """
        self.logger.info("Starting quality assessment")
        
        try:
            results = {
                'species': species,
                'reads_1': reads_1,
                'reads_2': reads_2,
                'timestamp': str(pd.Timestamp.now())
            }
            
            # Step 1: Basic read statistics
            self.logger.info("Calculating basic read statistics")
            read_stats = self._calculate_read_statistics(reads_1, reads_2)
            results['read_statistics'] = read_stats
            
            # Step 2: K-mer analysis
            self.logger.info("Performing k-mer analysis")
            kmer_analysis = self._perform_kmer_analysis(reads_1, reads_2)
            results['kmer_analysis'] = kmer_analysis
            
            # Step 3: CNN-based quality classification (if model available)
            if self.cnn_model:
                self.logger.info("Running CNN quality classification")
                cnn_results = self._run_cnn_classification(reads_1, reads_2)
                results['cnn_classification'] = cnn_results
            
            # Step 4: Contamination detection
            self.logger.info("Detecting contamination")
            contamination_results = self._detect_contamination(reads_1, reads_2, kmer_analysis)
            results['contamination'] = contamination_results
            
            # Step 5: Overall quality score
            quality_score = self._calculate_overall_quality_score(results)
            results['quality_score'] = quality_score
            
            # Step 6: Generate recommendations
            recommendations = self._generate_recommendations(results)
            results['recommendations'] = recommendations
            
            # Save results
            self.save_results(results)
            
            self.logger.info(f"Quality assessment completed. Overall score: {quality_score:.3f}")
            
            return results
            
        except Exception as e:
            self.logger.exception("Quality assessment failed")
            raise ModuleError(f"Quality assessment failed: {str(e)}")
    
    def _calculate_read_statistics(self, reads_1: str, reads_2: Optional[str]) -> Dict[str, Any]:
        """Calculate basic statistics for sequencing reads."""
        stats = {
            'forward_reads': {},
            'reverse_reads': {} if reads_2 else None,
            'paired_end': reads_2 is not None
        }
        
        # Analyze forward reads
        stats['forward_reads'] = self._analyze_single_file(reads_1)
        
        # Analyze reverse reads if present
        if reads_2:
            stats['reverse_reads'] = self._analyze_single_file(reads_2)
        
        # Calculate combined statistics
        total_reads = stats['forward_reads']['num_reads']
        if reads_2:
            total_reads += stats['reverse_reads']['num_reads']
        
        stats['total_reads'] = total_reads
        stats['estimated_coverage'] = self._estimate_coverage(stats)
        
        return stats
    
    def _analyze_single_file(self, fastq_path: str) -> Dict[str, Any]:
        """Analyze a single FASTQ file."""
        try:
            # Sample reads for analysis
            sampled_reads = []
            read_count = 0
            
            if fastq_path.endswith('.gz'):
                import gzip
                handle = gzip.open(fastq_path, 'rt')
            else:
                handle = open(fastq_path, 'r')
            
            with handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if len(sampled_reads) < self.sample_size:
                        sampled_reads.append(record)
                    read_count += 1
                    
                    # Stop after counting enough reads
                    if read_count > self.sample_size * 10:
                        break
            
            # Calculate statistics
            lengths = [len(read.seq) for read in sampled_reads]
            qualities = []
            gc_contents = []
            
            for read in sampled_reads:
                # Calculate average quality
                if hasattr(read, 'letter_annotations') and 'phred_quality' in read.letter_annotations:
                    avg_qual = np.mean(read.letter_annotations['phred_quality'])
                    qualities.append(avg_qual)
                
                # Calculate GC content
                gc_content = calculate_gc_content(read.seq)
                gc_contents.append(gc_content)
            
            return {
                'num_reads': read_count,
                'sampled_reads': len(sampled_reads),
                'avg_length': np.mean(lengths),
                'min_length': min(lengths),
                'max_length': max(lengths),
                'avg_quality': np.mean(qualities) if qualities else None,
                'avg_gc_content': np.mean(gc_contents),
                'length_distribution': {
                    'q25': np.percentile(lengths, 25),
                    'q50': np.percentile(lengths, 50),
                    'q75': np.percentile(lengths, 75)
                }
            }
            
        except Exception as e:
            raise ModuleError(f"Failed to analyze FASTQ file {fastq_path}: {str(e)}")
    
    def _perform_kmer_analysis(self, reads_1: str, reads_2: Optional[str]) -> Dict[str, Any]:
        """Perform k-mer analysis on reads."""
        try:
            all_kmers = {}
            
            # Analyze forward reads
            forward_kmers = self._count_kmers_in_file(reads_1)
            for kmer, count in forward_kmers.items():
                all_kmers[kmer] = all_kmers.get(kmer, 0) + count
            
            # Analyze reverse reads if present
            if reads_2:
                reverse_kmers = self._count_kmers_in_file(reads_2)
                for kmer, count in reverse_kmers.items():
                    all_kmers[kmer] = all_kmers.get(kmer, 0) + count
            
            # Calculate k-mer statistics
            kmer_counts = list(all_kmers.values())
            unique_kmers = len(all_kmers)
            total_kmers = sum(kmer_counts)
            
            # Identify high-frequency k-mers (potential repeats)
            high_freq_threshold = np.percentile(kmer_counts, 95)
            high_freq_kmers = {k: v for k, v in all_kmers.items() if v >= high_freq_threshold}
            
            return {
                'kmer_size': self.kmer_size,
                'unique_kmers': unique_kmers,
                'total_kmers': total_kmers,
                'avg_kmer_frequency': np.mean(kmer_counts),
                'kmer_frequency_distribution': {
                    'min': min(kmer_counts),
                    'max': max(kmer_counts),
                    'q25': np.percentile(kmer_counts, 25),
                    'q50': np.percentile(kmer_counts, 50),
                    'q75': np.percentile(kmer_counts, 75),
                    'q95': np.percentile(kmer_counts, 95)
                },
                'high_frequency_kmers': len(high_freq_kmers),
                'repetitive_content_estimate': len(high_freq_kmers) / unique_kmers if unique_kmers > 0 else 0
            }
            
        except Exception as e:
            raise ModuleError(f"K-mer analysis failed: {str(e)}")
    
    def _count_kmers_in_file(self, fastq_path: str) -> Dict[str, int]:
        """Count k-mers in a single FASTQ file."""
        kmer_counts = {}
        
        if fastq_path.endswith('.gz'):
            import gzip
            handle = gzip.open(fastq_path, 'rt')
        else:
            handle = open(fastq_path, 'r')
        
        with handle:
            read_count = 0
            for record in SeqIO.parse(handle, "fastq"):
                if read_count >= self.sample_size:
                    break
                
                sequence_kmers = count_kmers(record.seq, self.kmer_size)
                for kmer, count in sequence_kmers.items():
                    kmer_counts[kmer] = kmer_counts.get(kmer, 0) + count
                
                read_count += 1
        
        return kmer_counts
    
    def _run_cnn_classification(self, reads_1: str, reads_2: Optional[str]) -> Dict[str, Any]:
        """Run CNN-based quality classification."""
        try:
            # Sample sequences for classification
            sequences = []
            
            # Get sequences from forward reads
            sequences.extend(self._sample_sequences_from_file(reads_1, self.sample_size // 2))
            
            # Get sequences from reverse reads if present
            if reads_2:
                sequences.extend(self._sample_sequences_from_file(reads_2, self.sample_size // 2))
            
            # Run CNN classification
            predictions = self.cnn_model.predict_quality_scores(sequences)
            organellar_predictions = self.cnn_model.classify_organellar_reads(sequences)
            
            return {
                'num_sequences_classified': len(sequences),
                'avg_quality_score': np.mean(predictions),
                'quality_score_distribution': {
                    'min': np.min(predictions),
                    'max': np.max(predictions),
                    'q25': np.percentile(predictions, 25),
                    'q50': np.percentile(predictions, 50),
                    'q75': np.percentile(predictions, 75)
                },
                'organellar_fraction': np.mean(organellar_predictions),
                'high_quality_fraction': np.mean(predictions > self.confidence_threshold)
            }
            
        except Exception as e:
            raise ModuleError(f"CNN classification failed: {str(e)}")
    
    def _sample_sequences_from_file(self, fastq_path: str, num_sequences: int) -> List[str]:
        """Sample sequences from a FASTQ file."""
        sequences = []
        
        if fastq_path.endswith('.gz'):
            import gzip
            handle = gzip.open(fastq_path, 'rt')
        else:
            handle = open(fastq_path, 'r')
        
        with handle:
            for i, record in enumerate(SeqIO.parse(handle, "fastq")):
                if len(sequences) >= num_sequences:
                    break
                sequences.append(str(record.seq))
        
        return sequences
    
    def _detect_contamination(self, reads_1: str, reads_2: Optional[str], 
                           kmer_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Detect potential contamination in reads."""
        try:
            contamination_indicators = []
            
            # Check for unusual k-mer frequency distribution
            repetitive_content = kmer_analysis.get('repetitive_content_estimate', 0)
            if repetitive_content > 0.5:  # More than 50% repetitive
                contamination_indicators.append("High repetitive content detected")
            
            # Estimate contamination level
            contamination_level = min(repetitive_content * 100, 100)  # Convert to percentage
            
            return {
                'contamination_indicators': contamination_indicators,
                'estimated_contamination_level': contamination_level,
                'repetitive_content_fraction': repetitive_content,
                'contamination_risk': 'high' if contamination_level > 20 else 'medium' if contamination_level > 10 else 'low'
            }
            
        except Exception as e:
            raise ModuleError(f"Contamination detection failed: {str(e)}")
    
    def _estimate_coverage(self, stats: Dict[str, Any]) -> float:
        """Estimate sequencing coverage."""
        try:
            # Simplified coverage estimation
            # Assumes average organellar genome size of 150kb
            avg_genome_size = 150000  # 150kb
            
            total_bases = 0
            if 'forward_reads' in stats:
                total_bases += stats['forward_reads']['num_reads'] * stats['forward_reads']['avg_length']
            
            if stats.get('reverse_reads'):
                total_bases += stats['reverse_reads']['num_reads'] * stats['reverse_reads']['avg_length']
            
            estimated_coverage = total_bases / avg_genome_size
            return estimated_coverage
            
        except Exception as e:
            self.logger.warning(f"Coverage estimation failed: {str(e)}")
            return 0.0
    
    def _calculate_overall_quality_score(self, results: Dict[str, Any]) -> float:
        """Calculate overall quality score."""
        try:
            score_components = []
            
            # Read quality component
            if 'read_statistics' in results:
                read_stats = results['read_statistics']
                if 'forward_reads' in read_stats and read_stats['forward_reads'].get('avg_quality'):
                    # Normalize quality score (assuming Phred scale 0-40)
                    quality_score = min(read_stats['forward_reads']['avg_quality'] / 40.0, 1.0)
                    score_components.append(quality_score * 0.3)  # 30% weight
            
            # Coverage component
            coverage = results.get('read_statistics', {}).get('estimated_coverage', 0)
            coverage_score = min(coverage / 50.0, 1.0)  # Normalize to 50x coverage
            score_components.append(coverage_score * 0.2)  # 20% weight
            
            # CNN component (if available)
            if 'cnn_classification' in results:
                cnn_score = results['cnn_classification'].get('avg_quality_score', 0.5)
                score_components.append(cnn_score * 0.3)  # 30% weight
            
            # Contamination component (inverted)
            contamination_level = results.get('contamination', {}).get('estimated_contamination_level', 0)
            contamination_score = max(0, 1.0 - contamination_level / 100.0)
            score_components.append(contamination_score * 0.2)  # 20% weight
            
            # Calculate weighted average
            if score_components:
                overall_score = sum(score_components)
            else:
                overall_score = 0.5  # Default neutral score
            
            return min(max(overall_score, 0.0), 1.0)  # Clamp to [0, 1]
            
        except Exception as e:
            self.logger.warning(f"Quality score calculation failed: {str(e)}")
            return 0.5
    
    def _generate_recommendations(self, results: Dict[str, Any]) -> List[str]:
        """Generate processing recommendations based on quality assessment."""
        recommendations = []
        
        quality_score = results.get('quality_score', 0.5)
        contamination = results.get('contamination', {})
        
        if quality_score < 0.3:
            recommendations.append("Low overall quality detected. Consider quality filtering.")
        
        if contamination.get('contamination_risk') == 'high':
            recommendations.append("High contamination risk. Consider additional filtering steps.")
        
        coverage = results.get('read_statistics', {}).get('estimated_coverage', 0)
        if coverage < 10:
            recommendations.append("Low coverage detected. Results may be incomplete.")
        elif coverage > 200:
            recommendations.append("Very high coverage detected. Consider downsampling for efficiency.")
        
        if not recommendations:
            recommendations.append("Quality assessment passed. Proceed with assembly.")
        
        return recommendations