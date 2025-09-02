"""
Base module class for Organoverse workbench modules.
"""

import os
import time
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, Any, Optional

from ..core.logger import ModuleLogger, get_module_logger
from ..utils.exceptions import ModuleError
from ..utils.file_utils import ensure_directory


class BaseModule(ABC):
    """
    Abstract base class for all Organoverse modules.
    
    Provides common functionality for logging, configuration management,
    and output handling.
    """
    
    def __init__(self, config: Dict[str, Any], module_name: str):
        """
        Initialize base module.
        
        Args:
            config: Module configuration dictionary
            module_name: Name of the module
        """
        self.config = config
        self.module_name = module_name
        self.logger = get_module_logger(module_name)
        
        # Setup output directory
        self.output_dir = Path(config.get('output_dir', 'results'))
        self.module_output_dir = self.output_dir / module_name
        ensure_directory(self.module_output_dir)
        
        # Setup temporary directory
        self.temp_dir = Path(config.get('temp_dir', '/tmp/organoverse')) / module_name
        ensure_directory(self.temp_dir)
        
        # Common parameters
        self.threads = config.get('threads', 4)
        self.memory = config.get('memory', 16)
        self.verbose = config.get('verbose', False)
        
        self.logger.info(f"Initialized {module_name} module")
        self.logger.debug(f"Output directory: {self.module_output_dir}")
        self.logger.debug(f"Temp directory: {self.temp_dir}")
    
    @abstractmethod
    def run(self, *args, **kwargs) -> Dict[str, Any]:
        """
        Run the module's main functionality.
        
        Returns:
            Dictionary containing module results
        """
        pass
    
    def validate(self) -> bool:
        """
        Validate module configuration and dependencies.
        
        Returns:
            True if module is valid
        """
        try:
            # Check output directory is writable
            test_file = self.module_output_dir / "test_write"
            test_file.touch()
            test_file.unlink()
            
            # Check temp directory is writable
            test_file = self.temp_dir / "test_write"
            test_file.touch()
            test_file.unlink()
            
            return True
            
        except Exception as e:
            self.logger.error(f"Module validation failed: {str(e)}")
            return False
    
    def cleanup(self):
        """Clean up temporary files and resources."""
        try:
            if self.temp_dir.exists():
                import shutil
                shutil.rmtree(self.temp_dir)
            self.logger.debug("Cleanup completed")
        except Exception as e:
            self.logger.warning(f"Cleanup failed: {str(e)}")
    
    def save_results(self, results: Dict[str, Any], filename: str = "results.json"):
        """
        Save results to JSON file.
        
        Args:
            results: Results dictionary
            filename: Output filename
        """
        import json
        
        output_path = self.module_output_dir / filename
        
        try:
            with open(output_path, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            self.logger.info(f"Results saved to: {output_path}")
            
        except Exception as e:
            self.logger.error(f"Failed to save results: {str(e)}")
            raise ModuleError(f"Failed to save results: {str(e)}")
    
    def load_results(self, filename: str = "results.json") -> Dict[str, Any]:
        """
        Load results from JSON file.
        
        Args:
            filename: Input filename
        
        Returns:
            Results dictionary
        """
        import json
        
        input_path = self.module_output_dir / filename
        
        try:
            with open(input_path, 'r') as f:
                results = json.load(f)
            
            self.logger.debug(f"Results loaded from: {input_path}")
            return results
            
        except Exception as e:
            self.logger.error(f"Failed to load results: {str(e)}")
            raise ModuleError(f"Failed to load results: {str(e)}")
    
    def get_runtime_info(self) -> Dict[str, Any]:
        """
        Get runtime information for the module.
        
        Returns:
            Runtime information dictionary
        """
        return {
            'module_name': self.module_name,
            'output_dir': str(self.module_output_dir),
            'temp_dir': str(self.temp_dir),
            'threads': self.threads,
            'memory': self.memory,
            'config': self.config
        }