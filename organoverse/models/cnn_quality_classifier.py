"""
CNN-based quality classifier for sequencing reads.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from typing import List, Tuple, Optional
from pathlib import Path

from ..utils.exceptions import ModelError


class CNNQualityClassifier(nn.Module):
    """
    Convolutional Neural Network for read quality classification.
    
    This model classifies sequencing reads based on quality and
    identifies organellar vs nuclear sequences.
    """
    
    def __init__(self, model_path: Optional[str] = None, device: str = 'cpu'):
        """
        Initialize CNN quality classifier.
        
        Args:
            model_path: Path to pre-trained model
            device: Device for inference ('cpu' or 'cuda')
        """
        super(CNNQualityClassifier, self).__init__()
        
        self.device = torch.device(device)
        self.sequence_length = 150  # Standard read length
        self.vocab_size = 5  # A, T, G, C, N
        
        # Define CNN architecture
        self.embedding = nn.Embedding(self.vocab_size, 32)
        
        self.conv1 = nn.Conv1d(32, 64, kernel_size=7, padding=3)
        self.conv2 = nn.Conv1d(64, 128, kernel_size=5, padding=2)
        self.conv3 = nn.Conv1d(128, 256, kernel_size=3, padding=1)
        
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.dropout = nn.Dropout(0.5)
        
        # Quality prediction head
        self.quality_head = nn.Sequential(
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
        
        # Organellar classification head
        self.organellar_head = nn.Sequential(
            nn.Linear(256, 128),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Sigmoid()
        )
        
        self.to(self.device)
        
        # Load pre-trained model if provided
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
    
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass through the network.
        
        Args:
            x: Input tensor of shape (batch_size, sequence_length)
        
        Returns:
            Tuple of (quality_scores, organellar_predictions)
        """
        # Embedding
        x = self.embedding(x)  # (batch_size, seq_len, embed_dim)
        x = x.transpose(1, 2)  # (batch_size, embed_dim, seq_len)
        
        # Convolutional layers
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = F.relu(self.conv3(x))
        
        # Global pooling
        x = self.pool(x).squeeze(-1)  # (batch_size, 256)
        x = self.dropout(x)
        
        # Prediction heads
        quality_scores = self.quality_head(x)
        organellar_predictions = self.organellar_head(x)
        
        return quality_scores, organellar_predictions
    
    def sequence_to_tensor(self, sequences: List[str]) -> torch.Tensor:
        """
        Convert DNA sequences to tensor format.
        
        Args:
            sequences: List of DNA sequences
        
        Returns:
            Tensor of encoded sequences
        """
        # Nucleotide to index mapping
        nucleotide_map = {'A': 0, 'T': 1, 'G': 2, 'C': 3, 'N': 4}
        
        encoded_sequences = []
        
        for seq in sequences:
            # Pad or truncate sequence to fixed length
            seq = seq.upper()[:self.sequence_length]
            seq = seq.ljust(self.sequence_length, 'N')
            
            # Convert to indices
            encoded = [nucleotide_map.get(base, 4) for base in seq]
            encoded_sequences.append(encoded)
        
        return torch.tensor(encoded_sequences, dtype=torch.long, device=self.device)
    
    def predict_quality_scores(self, sequences: List[str]) -> np.ndarray:
        """
        Predict quality scores for sequences.
        
        Args:
            sequences: List of DNA sequences
        
        Returns:
            Array of quality scores (0-1)
        """
        try:
            self.eval()
            
            with torch.no_grad():
                # Convert sequences to tensor
                input_tensor = self.sequence_to_tensor(sequences)
                
                # Forward pass
                quality_scores, _ = self.forward(input_tensor)
                
                # Convert to numpy
                scores = quality_scores.cpu().numpy().flatten()
                
            return scores
            
        except Exception as e:
            raise ModelError(f"Quality prediction failed: {str(e)}")
    
    def classify_organellar_reads(self, sequences: List[str]) -> List[bool]:
        """
        Classify sequences as organellar or nuclear.
        
        Args:
            sequences: List of DNA sequences
        
        Returns:
            List of boolean predictions (True = organellar)
        """
        try:
            self.eval()
            
            with torch.no_grad():
                # Convert sequences to tensor
                input_tensor = self.sequence_to_tensor(sequences)
                
                # Forward pass
                _, organellar_predictions = self.forward(input_tensor)
                
                # Convert to boolean predictions (threshold at 0.5)
                predictions = (organellar_predictions.cpu().numpy().flatten() > 0.5).tolist()
                
            return predictions
            
        except Exception as e:
            raise ModelError(f"Organellar classification failed: {str(e)}")
    
    def save_model(self, model_path: str):
        """
        Save model to file.
        
        Args:
            model_path: Path to save model
        """
        try:
            torch.save({
                'model_state_dict': self.state_dict(),
                'sequence_length': self.sequence_length,
                'vocab_size': self.vocab_size
            }, model_path)
            
        except Exception as e:
            raise ModelError(f"Failed to save model: {str(e)}")
    
    def load_model(self, model_path: str):
        """
        Load model from file.
        
        Args:
            model_path: Path to model file
        """
        try:
            checkpoint = torch.load(model_path, map_location=self.device)
            self.load_state_dict(checkpoint['model_state_dict'])
            
            # Update model parameters if saved
            if 'sequence_length' in checkpoint:
                self.sequence_length = checkpoint['sequence_length']
            if 'vocab_size' in checkpoint:
                self.vocab_size = checkpoint['vocab_size']
                
        except Exception as e:
            raise ModelError(f"Failed to load model: {str(e)}")


def create_dummy_model(output_path: str):
    """
    Create a dummy pre-trained model for testing purposes.
    
    Args:
        output_path: Path to save the dummy model
    """
    model = CNNQualityClassifier()
    
    # Initialize with random weights (in practice, this would be trained)
    model.save_model(output_path)
    
    return model