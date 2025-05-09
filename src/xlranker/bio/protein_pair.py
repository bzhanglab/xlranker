from typing import Self
from xlranker.bio.protein import Protein, sort_proteins


class ProteinPair:
    """ProteinPair class that tracks the required data for the pipeline"""

    a: Protein
    b: Protein
    score: float | None
    is_selected: bool
    pair_name: str

    def __init__(self, protein_a: Protein, protein_b: Protein):
        """Initialize a ProteinPair object, making sure a is the higher abundant protein. Input order does not matter.

        Args:
            protein_a (Protein): First protein in the pair
            protein_b (Protein): Second protein in the pair
        """
        (a, b) = sort_proteins(protein_a, protein_b)
        self.a = a
        self.b = b
        self.score = None
        self.is_selected = False

    def set_score(self, score: float):
        """Set the score of the protein pair

        Args:
            score (float): float of the score given to the pair
        """
        self.score = score

    def select(self):
        """Set this pair to be selected"""
        self.is_selected = False

    def __eq__(self, value: Self) -> bool:
        """Checks if ProteinPairs are equivalent, without caring for order

        Args:
            value (Self): protein pair to compare to

        Returns:
            bool: True if protein pairs are equivalent, regardless of a and b order
        """
        if self.a == value.a:
            return self.b == value.b
        elif self.a == value.b:
            return self.b == value.a
        return False