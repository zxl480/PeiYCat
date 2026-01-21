import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops
from model import ReactionModel

# Amino acid dictionary for indexing
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
AA_TO_INDEX = {aa: i for i, aa in enumerate(AMINO_ACIDS)}

# One-hot encode the enzyme sequence
def one_hot_encode_sequence(sequence, n_gram=1):
    # Use the global AA_TO_INDEX dictionary for encoding
    # Create n-grams from the sequence
    n_grams = [sequence[i:i+n_gram] for i in range(len(sequence)-n_gram+1)]

    # Encode the n-grams into one-hot format
    one_hot_encoded = torch.zeros((len(n_grams), n_gram * 20))  # Shape (n, n_gram * 20)
    for i, n_gram_seq in enumerate(n_grams):
        for j, aa in enumerate(n_gram_seq):
            if aa in AA_TO_INDEX:
                one_hot_encoded[i, j * 20 + AA_TO_INDEX[aa]] = 1  # Set corresponding position to 1

    return one_hot_encoded

