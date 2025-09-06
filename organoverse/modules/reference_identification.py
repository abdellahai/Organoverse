"""
Reference Identification Module for Organoverse workbench.

This module identifies and retrieves closest reference genomes from RefSeq
using AI-enhanced similarity prediction and k-mer analysis.
"""

import os
import requests
import time
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

from .base_module import BaseModule
from ..utils.sequence_utils import count_kmers, parse_fasta, write_fasta
from ..utils.exceptions import ModuleError, DatabaseError


class ReferenceIdentificationModule(BaseModule):
    """
    AI-enhanced reference species identification and retrieval module.
    """
    
    def __init__(self, config: Dict[str, Any]):
        """
        Initialize Reference Identification module.
        
        Args:
            config: Module configuration
        """
        super().__init__(config, "reference_identification")
        
        # Database configuration
        self.ncbi_email = config.get('ncbi_email')
        self.ncbi_api_key = config.get('ncbi_api_key')
        self.cache_dir = Path(config.get('cache_dir', 'data/genbank_cache'))
        self.max_references = config.get('max_references', 5)
        self.similarity_threshold = config.get('similarity_threshold', 0.85)
        
        # Setup NCBI Entrez
        if self.ncbi_email:
            Entrez.email = self.ncbi_email
        if self.ncbi_api_key:
            Entrez.api_key = self.ncbi_api_key
        
        # Setup cache directory
        self.cache_dir.mkdir(parents=True, exist_ok=True)
    
    def run(self, species: str, organelle: str, kmer_profile: Optional[Dict] = None) -> Dict[str, Any]:
        """
        Run reference identification and retrieval.
        
        Args:
            species: Target species name
            organelle: Target organelle type
            kmer_profile: Optional k-mer profile from quality assessment
        
        Returns:
            Reference identification results
        """
        self.logger.info(f"Starting reference identification for {species} ({organelle})")
        
        try:
            results = {
                'species': species,
                'organelle': organelle,
                'timestamp': str(pd.Timestamp.now())
            }
            
            # Step 1: Search for candidate species (iterative taxonomy lookup)
            self.logger.info("Searching for candidate reference species from RefSeq")
            candidates = self._search_candidate_species(species, organelle)
            results['candidate_species'] = candidates
            
            # Step 2: Rank candidates taxonomically
            self.logger.info("Ranking candidates using taxonomic distance")
            ranked_candidates = self._rank_candidates_taxonomically(species, candidates)
            results['ranked_candidates'] = ranked_candidates
            
            # Step 3: Retrieve reference genomes
            self.logger.info("Retrieving reference genomes")
            references = self._retrieve_reference_genomes(ranked_candidates[:self.max_references])
            results['references'] = references
            
            # Step 4: Validate and filter references
            self.logger.info("Validating reference genomes")
            validated_references = self._validate_references(references, organelle)
            results['validated_references'] = validated_references
            
            # Save results
            self.save_results(results)
            
            self.logger.info(f"Reference identification completed. Found {len(validated_references)} references")
            
            return results
            
        except Exception as e:
            self.logger.exception("Reference identification failed")
            raise ModuleError(f"Reference identification failed: {str(e)}")
    
    def _search_candidate_species(self, species: str, organelle: str) -> List[Dict[str, Any]]:
        """Search for candidate reference species in RefSeq databases using iterative taxonomy."""
        try:
            # Split to extract genus and species names
            genus_species = species.split()[:2]  # Get genus and species
            if len(genus_species) < 2:
                raise ModuleError("Species name must include genus and species")
            
            genus, species_name = genus_species[0], genus_species[1]

            # Retrieve taxonomic lineage for iterative searching
            lineage = self._get_taxonomic_lineage(species)

            # Construct iterative search terms with RefSeq filter added
            search_terms = [
                # Current taxon: full species match
                f'"{genus} {species_name}"[Organism] AND {organelle}[Title] AND srcdb_refseq[PROP]',
                # Upper taxa: genus level
                f'{genus}[Organism] AND {organelle}[Title] AND srcdb_refseq[PROP]',
            ]

            # Add progressively broader taxa
            for rank in ["family", "order", "class", "phylum"]:
                name = next((n for r, n in lineage if r == rank), None)
                if name:
                    search_terms.append(
                        f'{name}[Organism] AND {organelle}[Title] AND srcdb_refseq[PROP]'
                    )

            # Fallback: broader organelle search among RefSeq entries
            search_terms.append(
                f'{genus}[Organism] AND organelle[Title] AND srcdb_refseq[PROP]'
            )
            
            for idx, search_term in enumerate(search_terms):
                self.logger.debug(f"Searching with term: {search_term}")
                candidates = []
                
                # Search nucleotide database
                search_handle = Entrez.esearch(
                    db="nucleotide",
                    term=search_term,
                    retmax=50,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if search_results["IdList"]:
                    # Fetch details for found sequences
                    id_list = search_results["IdList"][:20]  # Limit to top 20
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=id_list,
                        rettype="docsum",
                        retmode="xml"
                    )
                    fetch_results = Entrez.read(fetch_handle)
                    fetch_handle.close()
                    
                    # Process results
                    for doc in fetch_results:
                        candidate = {
                            'accession': doc.get('AccessionVersion', ''),
                            'title': doc.get('Title', ''),
                            'organism': doc.get('Organism', ''),
                            'length': int(doc.get('Length', 0)),
                            'search_term': search_term
                        }
                        
                        # Filter by organelle type and reasonable length
                        if self._is_valid_organelle_sequence(candidate, organelle):
                            candidates.append(candidate)
                    
                    # Remove duplicates based on accession
                    seen_accessions = set()
                    unique_candidates = []
                    for candidate in candidates:
                        if candidate['accession'] not in seen_accessions:
                            unique_candidates.append(candidate)
                            seen_accessions.add(candidate['accession'])
                    
                    # If candidates are found at the current taxonomic level, return immediately
                    if unique_candidates:
                        self.logger.debug(f"Found {len(unique_candidates)} candidates using term: {search_term}")
                        return unique_candidates
                
                # Rate limiting between search term iterations
                time.sleep(0.5)
            
            # If no candidates found using any iterative strategy, return empty list
            return []
            
        except Exception as e:
            raise DatabaseError(f"Candidate species search failed: {str(e)}")

    def _get_taxonomic_lineage(self, species: str) -> List[Tuple[str, str]]:
        """Retrieve taxonomic lineage for a species from NCBI."""
        try:
            handle = Entrez.esearch(db="taxonomy", term=species, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            if not record["IdList"]:
                return []
            tax_id = record["IdList"][0]
            fetch_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
            data = Entrez.read(fetch_handle)[0]
            fetch_handle.close()
            lineage = []
            for tax in data.get("LineageEx", []):
                lineage.append((tax.get("Rank", ""), tax.get("ScientificName", "")))
            lineage.append((data.get("Rank", ""), data.get("ScientificName", "")))
            return lineage
        except Exception as e:
            self.logger.debug(f"Taxonomy lookup failed for {species}: {e}")
            return []
    
    def _is_valid_organelle_sequence(self, candidate: Dict[str, Any], organelle: str) -> bool:
        """Check if a candidate sequence is a valid organelle sequence."""
        title = candidate.get('title', '').lower()
        length = candidate.get('length', 0)
        
        # Check title contains organelle keywords
        organelle_keywords = {
            'chloroplast': ['chloroplast', 'plastid', 'plastome'],
            'mitochondrion': ['mitochondrion', 'mitochondrial', 'mitochondria']
        }
        
        keywords = organelle_keywords.get(organelle, [])
        has_keyword = any(keyword in title for keyword in keywords)
        
        # Check reasonable length ranges
        length_ranges = {
            'chloroplast': (120000, 200000),  # 120-200kb
            'mitochondrion': (200000, 2000000)  # 200kb-2Mb
        }
        
        min_len, max_len = length_ranges.get(organelle, (50000, 5000000))
        valid_length = min_len <= length <= max_len
        
        # Exclude partial sequences
        is_complete = 'complete' in title or 'whole' in title
        is_partial = 'partial' in title or 'fragment' in title
        
        return has_keyword and valid_length and (is_complete or not is_partial)
    
    def _rank_candidates_taxonomically(self, target_species: str, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Rank candidates based on taxonomic similarity."""
        try:
            target_words = set(target_species.lower().split())
            
            for candidate in candidates:
                organism = candidate.get('organism', '').lower()
                organism_words = set(organism.split())

                # Calculate word overlap score
                overlap = len(target_words.intersection(organism_words))
                total_words = len(target_words.union(organism_words))

                if total_words > 0:
                    similarity_score = overlap / total_words
                else:
                    similarity_score = 0.0

                candidate['taxonomic_similarity_score'] = similarity_score

            # Add literature-based scores when multiple candidates are present
            if len(candidates) > 1:
                candidates = self._score_references_literature(target_species, candidates)

            # Sort by combined score
            ranked_candidates = sorted(
                candidates,
                key=lambda x: (
                    x.get('taxonomic_similarity_score', 0),
                    x.get('literature_score', 0),
                ),
                reverse=True,
            )

            return ranked_candidates

        except Exception as e:
            self.logger.warning(f"Taxonomic ranking failed: {str(e)}")
            return candidates

    def _score_references_literature(self, target_species: str, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Score candidate references based on literature co-occurrence."""
        for cand in candidates:
            organism = cand.get('organism', '')
            try:
                query = f'"{target_species}"[Title/Abstract] AND "{organism}"[Title/Abstract]'
                handle = Entrez.esearch(db="pubmed", term=query, retmax=0)
                record = Entrez.read(handle)
                handle.close()
                count = int(record.get("Count", 0))
            except Exception as e:
                self.logger.debug(f"Literature search failed for {organism}: {e}")
                count = 0
            cand['literature_score'] = count

        max_count = max((c['literature_score'] for c in candidates), default=0)
        if max_count > 0:
            for c in candidates:
                c['literature_score'] /= max_count

        return candidates
    
    def _retrieve_reference_genomes(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Retrieve reference genome sequences."""
        references = []
        
        for candidate in candidates:
            try:
                accession = candidate['accession']
                cache_file = self.cache_dir / f"{accession}.fasta"
                
                # Check cache first
                if cache_file.exists():
                    self.logger.debug(f"Loading {accession} from cache")
                    sequences = parse_fasta(str(cache_file))
                else:
                    self.logger.debug(f"Downloading {accession} from NCBI")
                    
                    # Download from NCBI
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=accession,
                        rettype="fasta",
                        retmode="text"
                    )
                    
                    # Save to cache
                    with open(cache_file, 'w') as f:
                        f.write(fetch_handle.read())
                    fetch_handle.close()
                    
                    # Parse sequences
                    sequences = parse_fasta(str(cache_file))
                    
                    # Rate limiting
                    time.sleep(0.5)
                
                if sequences:
                    reference = candidate.copy()
                    reference['sequence_file'] = str(cache_file)
                    reference['sequence_length'] = len(sequences[0].seq)
                    reference['sequence_id'] = sequences[0].id
                    references.append(reference)
                
            except Exception as e:
                self.logger.warning(f"Failed to retrieve {candidate['accession']}: {str(e)}")
                continue
        
        return references
    
    def _validate_references(self, references: List[Dict[str, Any]], organelle: str) -> List[Dict[str, Any]]:
        """Validate reference genomes."""
        validated = []
        
        for ref in references:
            try:
                # Check file exists and is readable
                if not os.path.exists(ref['sequence_file']):
                    continue
                
                # Parse and validate sequence
                sequences = parse_fasta(ref['sequence_file'])
                if not sequences:
                    continue
                
                seq_record = sequences[0]
                seq_length = len(seq_record.seq)
                
                # Validate length
                length_ranges = {
                    'chloroplast': (120000, 200000),
                    'mitochondrion': (200000, 2000000)
                }
                
                min_len, max_len = length_ranges.get(organelle, (50000, 5000000))
                if not (min_len <= seq_length <= max_len):
                    continue
                
                # Update reference info
                ref['validated'] = True
                ref['actual_length'] = seq_length
                validated.append(ref)
                
            except Exception as e:
                self.logger.warning(f"Validation failed for {ref['accession']}: {str(e)}")
                continue
        
        return validated
