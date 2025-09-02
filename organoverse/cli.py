#!/usr/bin/env python3
"""
Command Line Interface for Organoverse Workbench
"""

import click
import os
import sys
from pathlib import Path
from typing import Optional, List

from .core.pipeline import OrganoverseWorkbench
from .core.config import Config
from .core.logger import setup_logger
from .utils.validators import validate_fastq_files, validate_species_name
from .__init__ import __version__


@click.group()
@click.version_option(version=__version__)
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--config', '-c', type=click.Path(exists=True), help='Configuration file path')
@click.pass_context
def cli(ctx, verbose: bool, config: Optional[str]):
    """
    Organoverse: AI-Assisted Plant Organellar Genome Assembly Workbench
    
    A comprehensive framework for assembling plant mitochondrial and chloroplast 
    genomes using artificial intelligence and k-mer analysis.
    """
    ctx.ensure_object(dict)
    ctx.obj['verbose'] = verbose
    ctx.obj['config_file'] = config
    
    # Setup logging
    log_level = 'DEBUG' if verbose else 'INFO'
    setup_logger(level=log_level)


@cli.command()
@click.option('--species', '-s', required=True, type=str,
              help='Species name (e.g., "Arabidopsis thaliana")')
@click.option('--reads-1', '-1', required=True, type=click.Path(exists=True),
              help='Forward reads file (FASTQ format)')
@click.option('--reads-2', '-2', type=click.Path(exists=True),
              help='Reverse reads file (FASTQ format, for paired-end)')
@click.option('--organelle', '-o', 
              type=click.Choice(['chloroplast', 'mitochondrion', 'both']),
              default='chloroplast', help='Target organelle(s)')
@click.option('--output-dir', '-d', required=True, type=click.Path(),
              help='Output directory')
@click.option('--kmer-size', '-k', type=int, default=31,
              help='K-mer size for assembly (default: 31)')
@click.option('--coverage-cutoff', type=float, default=10.0,
              help='Coverage cutoff for assembly (default: 10.0)')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads (default: 4)')
@click.option('--memory', '-m', type=int, default=16,
              help='Memory limit in GB (default: 16)')
@click.option('--max-references', type=int, default=5,
              help='Maximum number of reference genomes to use (default: 5)')
@click.option('--skip-quality', is_flag=True,
              help='Skip quality assessment step')
@click.option('--skip-polishing', is_flag=True,
              help='Skip AI polishing step')
@click.option('--force', is_flag=True,
              help='Overwrite existing output directory')
@click.pass_context
def assemble(ctx, species: str, reads_1: str, reads_2: Optional[str], 
             organelle: str, output_dir: str, kmer_size: int, 
             coverage_cutoff: float, threads: int, memory: int,
             max_references: int, skip_quality: bool, skip_polishing: bool,
             force: bool):
    """
    Assemble organellar genomes from short-read sequencing data.
    
    This command runs the complete Organoverse pipeline including:
    1. Quality assessment using AI
    2. Reference species identification
    3. K-mer based assembly
    4. AI-powered completion and polishing
    5. Evaluation and comparison
    """
    
    try:
        # Validate inputs
        validate_species_name(species)
        validate_fastq_files(reads_1, reads_2)
        
        # Setup output directory
        output_path = Path(output_dir)
        if output_path.exists() and not force:
            click.echo(f"Error: Output directory {output_dir} already exists. Use --force to overwrite.")
            sys.exit(1)
        
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Load configuration
        config = Config()
        if ctx.obj.get('config_file'):
            config.load_from_file(ctx.obj['config_file'])
        
        # Update config with command line arguments
        config.update({
            'species': species,
            'reads_1': reads_1,
            'reads_2': reads_2,
            'organelle': organelle,
            'output_dir': str(output_path),
            'kmer_size': kmer_size,
            'coverage_cutoff': coverage_cutoff,
            'threads': threads,
            'memory': memory,
            'max_references': max_references,
            'skip_quality': skip_quality,
            'skip_polishing': skip_polishing,
            'verbose': ctx.obj.get('verbose', False)
        })
        
        # Initialize and run workbench
        workbench = OrganoverseWorkbench(config)
        
        click.echo(f"🧬 Starting Organoverse assembly for {species}")
        click.echo(f"📁 Output directory: {output_dir}")
        click.echo(f"🎯 Target organelle(s): {organelle}")
        
        # Run the complete pipeline
        results = workbench.run()
        
        # Display results
        click.echo("\n✅ Assembly completed successfully!")
        click.echo(f"📊 Results summary:")
        for org_type, result in results.items():
            if result['success']:
                click.echo(f"  {org_type}: {result['assembly_length']} bp")
                click.echo(f"    Quality score: {result['quality_score']:.2f}")
                click.echo(f"    Coverage: {result['coverage']:.1f}x")
            else:
                click.echo(f"  {org_type}: Failed - {result['error']}")
        
    except Exception as e:
        click.echo(f"❌ Error: {str(e)}", err=True)
        if ctx.obj.get('verbose'):
            import traceback
            traceback.print_exc()
        sys.exit(1)


@cli.command()
@click.option('--species', '-s', required=True, type=str,
              help='Species name for quality assessment')
@click.option('--reads-1', '-1', required=True, type=click.Path(exists=True),
              help='Forward reads file')
@click.option('--reads-2', '-2', type=click.Path(exists=True),
              help='Reverse reads file (optional)')
@click.option('--output-dir', '-d', required=True, type=click.Path(),
              help='Output directory')
@click.option('--threads', '-t', type=int, default=4,
              help='Number of threads')
@click.pass_context
def quality(ctx, species: str, reads_1: str, reads_2: Optional[str], 
           output_dir: str, threads: int):
    """
    Perform AI-based quality assessment of sequencing reads.
    """
    
    try:
        from .modules.quality_assessment import QualityAssessmentModule
        
        config = Config()
        config.update({
            'species': species,
            'reads_1': reads_1,
            'reads_2': reads_2,
            'output_dir': output_dir,
            'threads': threads
        })
        
        qa_module = QualityAssessmentModule(config)
        results = qa_module.run()
        
        click.echo("📊 Quality Assessment Results:")
        click.echo(f"  Overall Quality Score: {results['quality_score']:.2f}")
        click.echo(f"  Estimated Coverage: {results['coverage']:.1f}x")
        click.echo(f"  Contamination Level: {results['contamination']:.1f}%")
        
    except Exception as e:
        click.echo(f"❌ Error: {str(e)}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--species', '-s', required=True, type=str,
              help='Species name for reference search')
@click.option('--organelle', '-o', 
              type=click.Choice(['chloroplast', 'mitochondrion']),
              default='chloroplast', help='Target organelle')
@click.option('--max-references', type=int, default=5,
              help='Maximum number of references to retrieve')
@click.option('--output-dir', '-d', required=True, type=click.Path(),
              help='Output directory')
@click.pass_context
def references(ctx, species: str, organelle: str, max_references: int, output_dir: str):
    """
    Identify and retrieve reference genomes from GenBank/RefSeq.
    """
    
    try:
        from .modules.reference_identification import ReferenceIdentificationModule
        
        config = Config()
        config.update({
            'species': species,
            'organelle': organelle,
            'max_references': max_references,
            'output_dir': output_dir
        })
        
        ref_module = ReferenceIdentificationModule(config)
        references = ref_module.run()
        
        click.echo(f"🔍 Found {len(references)} reference genomes:")
        for i, ref in enumerate(references, 1):
            click.echo(f"  {i}. {ref['species']} ({ref['accession']}) - {ref['length']} bp")
        
    except Exception as e:
        click.echo(f"❌ Error: {str(e)}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--input-dir', '-i', required=True, type=click.Path(exists=True),
              help='Directory containing assembly results')
@click.option('--reference', '-r', type=click.Path(exists=True),
              help='Reference genome for comparison')
@click.option('--output-dir', '-d', required=True, type=click.Path(),
              help='Output directory for evaluation results')
@click.pass_context
def evaluate(ctx, input_dir: str, reference: Optional[str], output_dir: str):
    """
    Evaluate assembly quality and compare with references.
    """
    
    try:
        from .modules.evaluation import EvaluationModule
        
        config = Config()
        config.update({
            'input_dir': input_dir,
            'reference': reference,
            'output_dir': output_dir
        })
        
        eval_module = EvaluationModule(config)
        results = eval_module.run()
        
        click.echo("📈 Evaluation Results:")
        click.echo(f"  Assembly Length: {results['length']} bp")
        click.echo(f"  N50: {results['n50']} bp")
        click.echo(f"  Number of Contigs: {results['num_contigs']}")
        if reference:
            click.echo(f"  Identity to Reference: {results['identity']:.2f}%")
            click.echo(f"  Coverage of Reference: {results['coverage']:.2f}%")
        
    except Exception as e:
        click.echo(f"❌ Error: {str(e)}", err=True)
        sys.exit(1)


@cli.command()
@click.option('--model-type', 
              type=click.Choice(['cnn_quality', 'lstm_gap_filler', 'gnn_similarity']),
              required=True, help='Type of model to train')
@click.option('--training-data', '-d', required=True, type=click.Path(exists=True),
              help='Training data directory')
@click.option('--output-dir', '-o', required=True, type=click.Path(),
              help='Output directory for trained model')
@click.option('--epochs', '-e', type=int, default=100,
              help='Number of training epochs')
@click.option('--batch-size', '-b', type=int, default=32,
              help='Batch size for training')
@click.option('--learning-rate', '-lr', type=float, default=0.001,
              help='Learning rate')
@click.pass_context
def train(ctx, model_type: str, training_data: str, output_dir: str,
          epochs: int, batch_size: int, learning_rate: float):
    """
    Train AI models for quality assessment, gap filling, or similarity prediction.
    """
    
    try:
        from .models.trainer import ModelTrainer
        
        config = Config()
        config.update({
            'model_type': model_type,
            'training_data': training_data,
            'output_dir': output_dir,
            'epochs': epochs,
            'batch_size': batch_size,
            'learning_rate': learning_rate
        })
        
        trainer = ModelTrainer(config)
        trainer.train()
        
        click.echo(f"✅ Model training completed: {model_type}")
        click.echo(f"📁 Model saved to: {output_dir}")
        
    except Exception as e:
        click.echo(f"❌ Error: {str(e)}", err=True)
        sys.exit(1)


def main():
    """Main entry point for the CLI."""
    cli()


if __name__ == '__main__':
    main()