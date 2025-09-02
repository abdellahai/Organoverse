"""
Organoverse: AI-Assisted Plant Organellar Genome Assembly Workbench

A comprehensive framework for assembling plant mitochondrial and chloroplast genomes
using artificial intelligence and k-mer analysis.
"""

__version__ = "1.0.0"
__author__ = "Organoverse Development Team"
__email__ = "contact@organoverse.org"

from .core.pipeline import OrganoverseWorkbench
from .core.config import Config
from .core.logger import setup_logger

__all__ = [
    "OrganoverseWorkbench",
    "Config", 
    "setup_logger",
    "__version__"
]