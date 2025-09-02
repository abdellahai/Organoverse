"""
Processing modules for Organoverse workbench.
"""

from .quality_assessment import QualityAssessmentModule
from .reference_identification import ReferenceIdentificationModule
from .assembly import AssemblyModule
from .completion_polishing import CompletionPolishingModule
from .evaluation import EvaluationModule

__all__ = [
    'QualityAssessmentModule',
    'ReferenceIdentificationModule', 
    'AssemblyModule',
    'CompletionPolishingModule',
    'EvaluationModule'
]