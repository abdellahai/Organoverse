"""
Graph Neural Network for species similarity prediction.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import List, Dict, Optional
from pathlib import Path
import networkx as nx

from ..utils.exceptions import ModelError


class GNNSimilarityPredictor(nn.Module):
    """
    Graph Neural Network for predicting species similarity.
    
    This model builds phylogenetic graphs and predicts similarity scores
    between species based on taxonomic relationships.
    """
    
    def __init__(self, model_path: Optional[str] = None, device: str = 'cpu'):
        """
        Initialize GNN similarity predictor.
        
        Args:
            model_path: Path to pre-trained model
            device: Device for inference ('cpu' or 'cuda')
        """
        super(GNNSimilarityPredictor, self).__init__()
        
        self.device = torch.device(device)
        self.hidden_dim = 128
        self.num_layers = 3
        
        # Simple GNN implementation (in practice, would use PyTorch Geometric)
        self.node_embedding = nn.Embedding(10000, self.hidden_dim)  # Vocabulary for species
        
        self.gnn_layers = nn.ModuleList([
            nn.Linear(self.hidden_dim, self.hidden_dim) for _ in range(self.num_layers)
        ])
        
        self.similarity_head = nn.Sequential(
            nn.Linear(self.hidden_dim * 2, self.hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(self.hidden_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
        
        # Species name to index mapping (would be loaded from training data)
        self.species_vocab = {}
        self.vocab_size = 0
        
        self.to(self.device)
        
        # Load pre-trained model if provided
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
    
    def build_phylogenetic_graph(self, species_list: List[str]) -> nx.Graph:
        """
        Build phylogenetic graph from species list.
        
        Args:
            species_list: List of species names
        
        Returns:
            NetworkX graph representing phylogenetic relationships
        """
        graph = nx.Graph()
        
        # Add nodes for each species
        for species in species_list:
            graph.add_node(species)
        
        # Add edges based on taxonomic similarity (simplified)
        for i, species1 in enumerate(species_list):
            for j, species2 in enumerate(species_list[i+1:], i+1):
                similarity = self._calculate_taxonomic_similarity(species1, species2)
                if similarity > 0.1:  # Only add edges for similar species
                    graph.add_edge(species1, species2, weight=similarity)
        
        return graph
    
    def _calculate_taxonomic_similarity(self, species1: str, species2: str) -> float:
        """Calculate taxonomic similarity between two species."""
        # Simple word-based similarity
        words1 = set(species1.lower().split())
        words2 = set(species2.lower().split())
        
        if not words1 or not words2:
            return 0.0
        
        intersection = len(words1.intersection(words2))
        union = len(words1.union(words2))
        
        return intersection / union if union > 0 else 0.0
    
    def predict_similarity_scores(self, query_species: str, candidate_species: List[str]) -> np.ndarray:
        """
        Predict similarity scores between query species and candidates.
        
        Args:
            query_species: Query species name
            candidate_species: List of candidate species names
        
        Returns:
            Array of similarity scores
        """
        try:
            # For now, use simple taxonomic similarity as placeholder
            # In a real implementation, this would use the trained GNN
            scores = []
            
            for candidate in candidate_species:
                score = self._calculate_taxonomic_similarity(query_species, candidate)
                scores.append(score)
            
            return np.array(scores)
            
        except Exception as e:
            raise ModelError(f"Similarity prediction failed: {str(e)}")
    
    def save_model(self, model_path: str):
        """Save model to file."""
        try:
            torch.save({
                'model_state_dict': self.state_dict(),
                'species_vocab': self.species_vocab,
                'vocab_size': self.vocab_size,
                'hidden_dim': self.hidden_dim,
                'num_layers': self.num_layers
            }, model_path)
            
        except Exception as e:
            raise ModelError(f"Failed to save model: {str(e)}")
    
    def load_model(self, model_path: str):
        """Load model from file."""
        try:
            checkpoint = torch.load(model_path, map_location=self.device)
            self.load_state_dict(checkpoint['model_state_dict'])
            
            if 'species_vocab' in checkpoint:
                self.species_vocab = checkpoint['species_vocab']
            if 'vocab_size' in checkpoint:
                self.vocab_size = checkpoint['vocab_size']
                
        except Exception as e:
            raise ModelError(f"Failed to load model: {str(e)}")


def create_dummy_model(output_path: str):
    """Create a dummy pre-trained model for testing purposes."""
    model = GNNSimilarityPredictor()
    model.save_model(output_path)
    return model