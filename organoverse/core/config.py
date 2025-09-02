"""
Configuration management for Organoverse workbench.
"""

import os
import yaml
import toml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class ModelConfig:
    """Configuration for AI/ML models."""
    cnn_quality_model: str = "models/cnn_quality_classifier.pth"
    lstm_gap_filler_model: str = "models/lstm_gap_filler.pth"
    gnn_similarity_model: str = "models/gnn_similarity_predictor.pth"
    model_device: str = "auto"  # auto, cpu, cuda
    batch_size: int = 32
    confidence_threshold: float = 0.8


@dataclass
class AssemblyConfig:
    """Configuration for assembly parameters."""
    kmer_size: int = 31
    coverage_cutoff: float = 10.0
    min_contig_length: int = 500
    max_iterations: int = 10
    assembler: str = "spades"  # spades, novoplasty, getorganelle
    polish_rounds: int = 3


@dataclass
class DatabaseConfig:
    """Configuration for database access."""
    ncbi_email: Optional[str] = None
    ncbi_api_key: Optional[str] = None
    genbank_cache_dir: str = "data/genbank_cache"
    max_references: int = 5
    similarity_threshold: float = 0.85


@dataclass
class ComputeConfig:
    """Configuration for computational resources."""
    threads: int = 4
    memory: int = 16  # GB
    gpu_enabled: bool = True
    temp_dir: str = "/tmp/organoverse"
    cleanup_temp: bool = True


@dataclass
class Config:
    """Main configuration class for Organoverse workbench."""
    
    # Input/Output
    species: str = ""
    reads_1: str = ""
    reads_2: Optional[str] = None
    organelle: str = "chloroplast"  # chloroplast, mitochondrion, both
    output_dir: str = "results"
    
    # Processing options
    skip_quality: bool = False
    skip_polishing: bool = False
    force_overwrite: bool = False
    verbose: bool = False
    
    # Sub-configurations
    models: ModelConfig = field(default_factory=ModelConfig)
    assembly: AssemblyConfig = field(default_factory=AssemblyConfig)
    database: DatabaseConfig = field(default_factory=DatabaseConfig)
    compute: ComputeConfig = field(default_factory=ComputeConfig)
    
    # Additional parameters
    custom_params: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self):
        """Post-initialization setup."""
        # Set up paths
        self._setup_paths()
        
        # Auto-detect device for models
        if self.models.model_device == "auto":
            self.models.model_device = self._detect_device()
    
    def _setup_paths(self):
        """Setup default paths relative to package location."""
        package_dir = Path(__file__).parent.parent
        
        # Update model paths if they're relative
        if not os.path.isabs(self.models.cnn_quality_model):
            self.models.cnn_quality_model = str(package_dir / self.models.cnn_quality_model)
        if not os.path.isabs(self.models.lstm_gap_filler_model):
            self.models.lstm_gap_filler_model = str(package_dir / self.models.lstm_gap_filler_model)
        if not os.path.isabs(self.models.gnn_similarity_model):
            self.models.gnn_similarity_model = str(package_dir / self.models.gnn_similarity_model)
        
        # Update cache directory
        if not os.path.isabs(self.database.genbank_cache_dir):
            self.database.genbank_cache_dir = str(package_dir / self.database.genbank_cache_dir)
    
    def _detect_device(self) -> str:
        """Auto-detect the best available device for model inference."""
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
        except ImportError:
            pass
        
        try:
            import tensorflow as tf
            if tf.config.list_physical_devices('GPU'):
                return "gpu"
        except ImportError:
            pass
        
        return "cpu"
    
    def load_from_file(self, config_path: str):
        """Load configuration from YAML or TOML file."""
        config_path = Path(config_path)
        
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        if config_path.suffix.lower() in ['.yaml', '.yml']:
            with open(config_path, 'r') as f:
                data = yaml.safe_load(f)
        elif config_path.suffix.lower() == '.toml':
            with open(config_path, 'r') as f:
                data = toml.load(f)
        else:
            raise ValueError(f"Unsupported configuration file format: {config_path.suffix}")
        
        self.update(data)
    
    def save_to_file(self, config_path: str):
        """Save current configuration to file."""
        config_path = Path(config_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to dictionary
        data = self.to_dict()
        
        if config_path.suffix.lower() in ['.yaml', '.yml']:
            with open(config_path, 'w') as f:
                yaml.dump(data, f, default_flow_style=False, indent=2)
        elif config_path.suffix.lower() == '.toml':
            with open(config_path, 'w') as f:
                toml.dump(data, f)
        else:
            raise ValueError(f"Unsupported configuration file format: {config_path.suffix}")
    
    def update(self, data: Dict[str, Any]):
        """Update configuration with dictionary data."""
        for key, value in data.items():
            if hasattr(self, key):
                if isinstance(getattr(self, key), (ModelConfig, AssemblyConfig, DatabaseConfig, ComputeConfig)):
                    # Update nested configuration objects
                    nested_config = getattr(self, key)
                    if isinstance(value, dict):
                        for nested_key, nested_value in value.items():
                            if hasattr(nested_config, nested_key):
                                setattr(nested_config, nested_key, nested_value)
                else:
                    setattr(self, key, value)
            else:
                # Store in custom_params
                self.custom_params[key] = value
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        result = {}
        
        # Basic attributes
        for key in ['species', 'reads_1', 'reads_2', 'organelle', 'output_dir',
                   'skip_quality', 'skip_polishing', 'force_overwrite', 'verbose']:
            result[key] = getattr(self, key)
        
        # Nested configurations
        result['models'] = {
            'cnn_quality_model': self.models.cnn_quality_model,
            'lstm_gap_filler_model': self.models.lstm_gap_filler_model,
            'gnn_similarity_model': self.models.gnn_similarity_model,
            'model_device': self.models.model_device,
            'batch_size': self.models.batch_size,
            'confidence_threshold': self.models.confidence_threshold
        }
        
        result['assembly'] = {
            'kmer_size': self.assembly.kmer_size,
            'coverage_cutoff': self.assembly.coverage_cutoff,
            'min_contig_length': self.assembly.min_contig_length,
            'max_iterations': self.assembly.max_iterations,
            'assembler': self.assembly.assembler,
            'polish_rounds': self.assembly.polish_rounds
        }
        
        result['database'] = {
            'ncbi_email': self.database.ncbi_email,
            'ncbi_api_key': self.database.ncbi_api_key,
            'genbank_cache_dir': self.database.genbank_cache_dir,
            'max_references': self.database.max_references,
            'similarity_threshold': self.database.similarity_threshold
        }
        
        result['compute'] = {
            'threads': self.compute.threads,
            'memory': self.compute.memory,
            'gpu_enabled': self.compute.gpu_enabled,
            'temp_dir': self.compute.temp_dir,
            'cleanup_temp': self.compute.cleanup_temp
        }
        
        # Custom parameters
        result.update(self.custom_params)
        
        return result
    
    def get_module_config(self, module_name: str) -> Dict[str, Any]:
        """Get configuration specific to a module."""
        base_config = {
            'output_dir': self.output_dir,
            'threads': self.compute.threads,
            'memory': self.compute.memory,
            'verbose': self.verbose,
            'temp_dir': self.compute.temp_dir
        }
        
        if module_name == 'quality_assessment':
            base_config.update({
                'model_path': self.models.cnn_quality_model,
                'device': self.models.model_device,
                'batch_size': self.models.batch_size,
                'confidence_threshold': self.models.confidence_threshold
            })
        
        elif module_name == 'reference_identification':
            base_config.update({
                'ncbi_email': self.database.ncbi_email,
                'ncbi_api_key': self.database.ncbi_api_key,
                'cache_dir': self.database.genbank_cache_dir,
                'max_references': self.database.max_references,
                'similarity_threshold': self.database.similarity_threshold,
                'model_path': self.models.gnn_similarity_model
            })
        
        elif module_name == 'assembly':
            base_config.update({
                'kmer_size': self.assembly.kmer_size,
                'coverage_cutoff': self.assembly.coverage_cutoff,
                'min_contig_length': self.assembly.min_contig_length,
                'assembler': self.assembly.assembler
            })
        
        elif module_name == 'completion_polishing':
            base_config.update({
                'lstm_model_path': self.models.lstm_gap_filler_model,
                'device': self.models.model_device,
                'polish_rounds': self.assembly.polish_rounds,
                'max_iterations': self.assembly.max_iterations
            })
        
        return base_config