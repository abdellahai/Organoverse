"""
Main pipeline orchestrator for Organoverse workbench.
"""

import os
import time
from pathlib import Path
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

from .config import Config
from .logger import ModuleLogger, get_module_logger
from ..modules.quality_assessment import QualityAssessmentModule
from ..modules.reference_identification import ReferenceIdentificationModule
from ..modules.assembly import AssemblyModule
from ..modules.completion_polishing import CompletionPolishingModule
from ..modules.evaluation import EvaluationModule
from ..utils.exceptions import OrganoverseError, ModuleError


@dataclass
class PipelineResult:
    """Results from pipeline execution."""
    success: bool
    organelle_type: str
    assembly_path: Optional[str] = None
    assembly_length: Optional[int] = None
    quality_score: Optional[float] = None
    coverage: Optional[float] = None
    n50: Optional[int] = None
    num_contigs: Optional[int] = None
    runtime: Optional[float] = None
    error_message: Optional[str] = None
    intermediate_results: Dict[str, Any] = None


class OrganoverseWorkbench:
    """
    Main orchestrator for the Organoverse pipeline.
    
    Coordinates execution of all modules in the correct order and manages
    data flow between components.
    """
    
    def __init__(self, config: Config):
        """
        Initialize the workbench.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.logger = get_module_logger("workbench")
        
        # Initialize modules
        self.modules = {}
        self._initialize_modules()
        
        # Setup output directory
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save configuration
        config_path = self.output_dir / "config.yaml"
        self.config.save_to_file(str(config_path))
        
        self.logger.info(f"Initialized Organoverse workbench")
        self.logger.info(f"Output directory: {self.output_dir}")
        self.logger.info(f"Target species: {config.species}")
        self.logger.info(f"Target organelle(s): {config.organelle}")
    
    def _initialize_modules(self):
        """Initialize all pipeline modules."""
        try:
            self.modules['quality_assessment'] = QualityAssessmentModule(
                self.config.get_module_config('quality_assessment')
            )
            
            self.modules['reference_identification'] = ReferenceIdentificationModule(
                self.config.get_module_config('reference_identification')
            )
            
            self.modules['assembly'] = AssemblyModule(
                self.config.get_module_config('assembly')
            )
            
            self.modules['completion_polishing'] = CompletionPolishingModule(
                self.config.get_module_config('completion_polishing')
            )
            
            self.modules['evaluation'] = EvaluationModule(
                self.config.get_module_config('evaluation')
            )
            
            self.logger.info("All modules initialized successfully")
            
        except Exception as e:
            self.logger.error(f"Failed to initialize modules: {str(e)}")
            raise OrganoverseError(f"Module initialization failed: {str(e)}")
    
    def validate_modules(self) -> bool:
        """
        Validate that all modules are properly configured.
        
        Returns:
            True if all modules are valid
        """
        try:
            for module_name, module in self.modules.items():
                if hasattr(module, 'validate'):
                    if not module.validate():
                        self.logger.error(f"Module validation failed: {module_name}")
                        return False
                self.logger.debug(f"Module validated: {module_name}")
            
            self.logger.info("All modules validated successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Module validation error: {str(e)}")
            return False
    
    def run(self) -> Dict[str, PipelineResult]:
        """
        Run the complete Organoverse pipeline.
        
        Returns:
            Dictionary of results for each organelle type
        """
        start_time = time.time()
        results = {}
        
        try:
            # Validate modules before starting
            if not self.validate_modules():
                raise OrganoverseError("Module validation failed")
            
            # Determine target organelles
            organelles = self._get_target_organelles()
            
            # Run pipeline for each organelle
            for organelle in organelles:
                self.logger.info(f"Starting pipeline for {organelle}")
                
                try:
                    result = self._run_single_organelle_pipeline(organelle)
                    result.runtime = time.time() - start_time
                    results[organelle] = result
                    
                    if result.success:
                        self.logger.info(f"Pipeline completed successfully for {organelle}")
                    else:
                        self.logger.error(f"Pipeline failed for {organelle}: {result.error_message}")
                
                except Exception as e:
                    self.logger.exception(f"Pipeline error for {organelle}")
                    results[organelle] = PipelineResult(
                        success=False,
                        organelle_type=organelle,
                        error_message=str(e),
                        runtime=time.time() - start_time
                    )
            
            # Generate summary report
            self._generate_summary_report(results)
            
            return results
            
        except Exception as e:
            self.logger.exception("Critical pipeline error")
            raise OrganoverseError(f"Pipeline execution failed: {str(e)}")
    
    def _get_target_organelles(self) -> List[str]:
        """Get list of target organelles based on configuration."""
        if self.config.organelle == "both":
            return ["chloroplast", "mitochondrion"]
        else:
            return [self.config.organelle]
    
    def _run_single_organelle_pipeline(self, organelle: str) -> PipelineResult:
        """
        Run pipeline for a single organelle type.
        
        Args:
            organelle: Target organelle type
        
        Returns:
            Pipeline result
        """
        intermediate_results = {}
        
        try:
            # Update configuration for current organelle
            current_config = self.config
            current_config.organelle = organelle
            
            # Step 1: Quality Assessment
            if not self.config.skip_quality:
                self.logger.info("Step 1: Quality Assessment")
                qa_result = self.modules['quality_assessment'].run(
                    reads_1=self.config.reads_1,
                    reads_2=self.config.reads_2,
                    species=self.config.species
                )
                intermediate_results['quality_assessment'] = qa_result
                
                # Check if quality is sufficient
                if qa_result['quality_score'] < 0.5:
                    self.logger.warning(f"Low quality score: {qa_result['quality_score']}")
            
            # Step 2: Reference Identification
            self.logger.info("Step 2: Reference Identification")
            ref_result = self.modules['reference_identification'].run(
                species=self.config.species,
                organelle=organelle,
                kmer_profile=intermediate_results.get('quality_assessment', {}).get('kmer_profile')
            )
            intermediate_results['reference_identification'] = ref_result
            
            if not ref_result['references']:
                raise ModuleError("No suitable reference genomes found")
            
            # Step 3: Assembly
            self.logger.info("Step 3: Assembly")
            assembly_result = self.modules['assembly'].run(
                reads_1=self.config.reads_1,
                reads_2=self.config.reads_2,
                references=ref_result['references'],
                organelle=organelle
            )
            intermediate_results['assembly'] = assembly_result
            
            # Step 4: Completion and Polishing
            if not self.config.skip_polishing:
                self.logger.info("Step 4: Completion and Polishing")
                polish_result = self.modules['completion_polishing'].run(
                    assembly=assembly_result['assembly_path'],
                    reads_1=self.config.reads_1,
                    reads_2=self.config.reads_2,
                    references=ref_result['references']
                )
                intermediate_results['completion_polishing'] = polish_result
                final_assembly = polish_result['polished_assembly']
            else:
                final_assembly = assembly_result['assembly_path']
            
            # Step 5: Evaluation
            self.logger.info("Step 5: Evaluation")
            eval_result = self.modules['evaluation'].run(
                assembly=final_assembly,
                references=ref_result['references'],
                organelle=organelle
            )
            intermediate_results['evaluation'] = eval_result
            
            # Create successful result
            return PipelineResult(
                success=True,
                organelle_type=organelle,
                assembly_path=final_assembly,
                assembly_length=eval_result.get('assembly_length'),
                quality_score=eval_result.get('quality_score'),
                coverage=eval_result.get('coverage'),
                n50=eval_result.get('n50'),
                num_contigs=eval_result.get('num_contigs'),
                intermediate_results=intermediate_results
            )
            
        except Exception as e:
            self.logger.exception(f"Pipeline error for {organelle}")
            return PipelineResult(
                success=False,
                organelle_type=organelle,
                error_message=str(e),
                intermediate_results=intermediate_results
            )
    
    def _generate_summary_report(self, results: Dict[str, PipelineResult]):
        """Generate summary report of pipeline results."""
        report_path = self.output_dir / "summary_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("Organoverse Pipeline Summary Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Species: {self.config.species}\n")
            f.write(f"Input reads: {self.config.reads_1}")
            if self.config.reads_2:
                f.write(f", {self.config.reads_2}")
            f.write("\n\n")
            
            for organelle, result in results.items():
                f.write(f"{organelle.capitalize()} Results:\n")
                f.write("-" * 20 + "\n")
                
                if result.success:
                    f.write(f"Status: SUCCESS\n")
                    f.write(f"Assembly length: {result.assembly_length} bp\n")
                    f.write(f"Quality score: {result.quality_score:.3f}\n")
                    f.write(f"Coverage: {result.coverage:.1f}x\n")
                    f.write(f"N50: {result.n50} bp\n")
                    f.write(f"Number of contigs: {result.num_contigs}\n")
                    f.write(f"Assembly file: {result.assembly_path}\n")
                else:
                    f.write(f"Status: FAILED\n")
                    f.write(f"Error: {result.error_message}\n")
                
                f.write(f"Runtime: {result.runtime:.1f} seconds\n\n")
        
        self.logger.info(f"Summary report saved to: {report_path}")