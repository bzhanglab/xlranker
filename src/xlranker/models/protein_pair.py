from xlranker.models.protein import Protein, set_protein_order


class ProteinPair:
    """ProteinPair class that tracks the required data for the pipeline"""

    a: Protein
    b: Protein
    score: float | None
    is_selected: bool

    def __init__(self, protein_a: Protein, protein_b: Protein):
        """Initialize a ProteinPair object, making sure a is the higher abundant protein. Input order does not matter.

        Args:
            protein_a (Protein): First protein in the pair
            protein_b (Protein): Second protein in the pair
        """
        (a, b) = set_protein_order(protein_a, protein_b)
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
