"""
LSTM-based gap filler for assembly completion.
"""

import torch
import torch.nn as nn
import numpy as np
from typing import List, Tuple, Optional
from pathlib import Path

from ..utils.exceptions import ModelError


class LSTMGapFiller(nn.Module):
    """
    LSTM-based model for filling gaps in genome assemblies.
    
    This model predicts sequences to fill gaps based on left and right context.
    """
    
    def __init__(self, model_path: Optional[str] = None, device: str = 'cpu'):
        """
        Initialize LSTM gap filler.
        
        Args:
            model_path: Path to pre-trained model
            device: Device for inference ('cpu' or 'cuda')
        """
        super(LSTMGapFiller, self).__init__()
        
        self.device = torch.device(device)
        self.vocab_size = 5  # A, T, G, C, N
        self.hidden_size = 256
        self.num_layers = 2
        self.context_length = 100  # Length of context on each side
        
        # Embedding layer
        self.embedding = nn.Embedding(self.vocab_size, 64)
        
        # Bidirectional LSTM for context encoding
        self.context_lstm = nn.LSTM(
            input_size=64,
            hidden_size=self.hidden_size,
            num_layers=self.num_layers,
            bidirectional=True,
            batch_first=True,
            dropout=0.3
        )
        
        # LSTM for sequence generation
        self.generator_lstm = nn.LSTM(
            input_size=64,
            hidden_size=self.hidden_size * 2,  # Match bidirectional output
            num_layers=self.num_layers,
            batch_first=True,
            dropout=0.3
        )
        
        # Output layer
        self.output_layer = nn.Linear(self.hidden_size * 2, self.vocab_size)
        self.dropout = nn.Dropout(0.3)
        
        self.to(self.device)
        
        # Load pre-trained model if provided
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
    
    def forward(self, left_context: torch.Tensor, right_context: torch.Tensor, 
                target_length: int) -> torch.Tensor:
        """
        Forward pass through the network.
        
        Args:
            left_context: Left context sequences (batch_size, context_length)
            right_context: Right context sequences (batch_size, context_length)
            target_length: Length of sequence to generate
        
        Returns:
            Generated sequences (batch_size, target_length, vocab_size)
        """
        batch_size = left_context.size(0)
        
        # Embed contexts
        left_embedded = self.embedding(left_context)
        right_embedded = self.embedding(right_context)
        
        # Encode contexts with bidirectional LSTM
        left_encoded, _ = self.context_lstm(left_embedded)
        right_encoded, _ = self.context_lstm(right_embedded)
        
        # Use last hidden state from left context and first from right context
        left_final = left_encoded[:, -1, :]  # (batch_size, hidden_size*2)
        right_final = right_encoded[:, 0, :]  # (batch_size, hidden_size*2)
        
        # Combine contexts
        combined_context = (left_final + right_final) / 2
        
        # Initialize generator LSTM hidden state
        h_0 = combined_context.unsqueeze(0).repeat(self.num_layers, 1, 1)
        c_0 = torch.zeros_like(h_0)
        
        # Generate sequence
        outputs = []
        hidden = (h_0, c_0)
        
        # Start with a special start token (index 0 = 'A')
        current_input = torch.zeros(batch_size, 1, dtype=torch.long, device=self.device)
        
        for _ in range(target_length):
            # Embed current input
            embedded_input = self.embedding(current_input)
            
            # LSTM step
            lstm_out, hidden = self.generator_lstm(embedded_input, hidden)
            
            # Apply dropout and output layer
            output = self.dropout(lstm_out)
            output = self.output_layer(output)  # (batch_size, 1, vocab_size)
            
            outputs.append(output)
            
            # Use predicted token as next input (teacher forcing during training)
            current_input = torch.argmax(output, dim=-1)
        
        # Concatenate outputs
        generated_sequence = torch.cat(outputs, dim=1)  # (batch_size, target_length, vocab_size)
        
        return generated_sequence
    
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
            # Pad or truncate sequence to context length
            seq = seq.upper()[:self.context_length]
            seq = seq.ljust(self.context_length, 'N')
            
            # Convert to indices
            encoded = [nucleotide_map.get(base, 4) for base in seq]
            encoded_sequences.append(encoded)
        
        return torch.tensor(encoded_sequences, dtype=torch.long, device=self.device)
    
    def tensor_to_sequence(self, tensor: torch.Tensor) -> List[str]:
        """
        Convert tensor back to DNA sequences.
        
        Args:
            tensor: Tensor of encoded sequences
        
        Returns:
            List of DNA sequences
        """
        # Index to nucleotide mapping
        index_map = {0: 'A', 1: 'T', 2: 'G', 3: 'C', 4: 'N'}
        
        sequences = []
        
        # Convert tensor to numpy
        if tensor.dim() == 3:  # (batch_size, seq_len, vocab_size)
            tensor = torch.argmax(tensor, dim=-1)  # (batch_size, seq_len)
        
        tensor_np = tensor.cpu().numpy()
        
        for seq_indices in tensor_np:
            sequence = ''.join([index_map[idx] for idx in seq_indices])
            sequences.append(sequence)
        
        return sequences
    
    def predict_gap_sequences(self, left_contexts: List[str], right_contexts: List[str], 
                            gap_size: int) -> List[str]:
        """
        Predict sequences to fill gaps.
        
        Args:
            left_contexts: Left context sequences
            right_contexts: Right context sequences
            gap_size: Size of gap to fill
        
        Returns:
            List of predicted gap sequences
        """
        try:
            self.eval()
            
            with torch.no_grad():
                # Convert contexts to tensors
                left_tensor = self.sequence_to_tensor(left_contexts)
                right_tensor = self.sequence_to_tensor(right_contexts)
                
                # Generate sequences
                generated_tensor = self.forward(left_tensor, right_tensor, gap_size)
                
                # Convert back to sequences
                gap_sequences = self.tensor_to_sequence(generated_tensor)
                
            return gap_sequences
            
        except Exception as e:
            raise ModelError(f"Gap sequence prediction failed: {str(e)}")
    
    def score_gap_predictions(self, predictions: List[str], reads_data: Optional[List] = None) -> List[float]:
        """
        Score gap predictions based on various criteria.
        
        Args:
            predictions: List of predicted gap sequences
            reads_data: Optional read data for validation
        
        Returns:
            List of confidence scores (0-1)
        """
        try:
            scores = []
            
            for pred in predictions:
                score = 0.5  # Base score
                
                # Score based on sequence complexity
                if len(set(pred)) > 1:  # Not all same nucleotide
                    score += 0.2
                
                # Score based on GC content (should be reasonable)
                gc_content = (pred.count('G') + pred.count('C')) / len(pred) if len(pred) > 0 else 0
                if 0.3 <= gc_content <= 0.7:  # Reasonable GC content
                    score += 0.2
                
                # Avoid long homopolymer runs
                max_run = self._max_homopolymer_run(pred)
                if max_run < len(pred) * 0.5:  # Less than 50% homopolymer
                    score += 0.1
                
                scores.append(min(score, 1.0))
            
            return scores
            
        except Exception as e:
            raise ModelError(f"Gap prediction scoring failed: {str(e)}")
    
    def _max_homopolymer_run(self, sequence: str) -> int:
        """Calculate maximum homopolymer run length."""
        if not sequence:
            return 0
        
        max_run = 1
        current_run = 1
        
        for i in range(1, len(sequence)):
            if sequence[i] == sequence[i-1]:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
        
        return max_run
    
    def save_model(self, model_path: str):
        """
        Save model to file.
        
        Args:
            model_path: Path to save model
        """
        try:
            torch.save({
                'model_state_dict': self.state_dict(),
                'vocab_size': self.vocab_size,
                'hidden_size': self.hidden_size,
                'num_layers': self.num_layers,
                'context_length': self.context_length
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
            if 'vocab_size' in checkpoint:
                self.vocab_size = checkpoint['vocab_size']
            if 'hidden_size' in checkpoint:
                self.hidden_size = checkpoint['hidden_size']
            if 'num_layers' in checkpoint:
                self.num_layers = checkpoint['num_layers']
            if 'context_length' in checkpoint:
                self.context_length = checkpoint['context_length']
                
        except Exception as e:
            raise ModelError(f"Failed to load model: {str(e)}")


def create_dummy_model(output_path: str):
    """
    Create a dummy pre-trained model for testing purposes.
    
    Args:
        output_path: Path to save the dummy model
    """
    model = LSTMGapFiller()
    
    # Initialize with random weights (in practice, this would be trained)
    model.save_model(output_path)
    
    return model