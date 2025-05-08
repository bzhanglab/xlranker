from dataclasses import dataclass


@dataclass
class Protein:
    """Protein class that has the name and abundance for the protein"""

    name: str
    abundance: float | None


from xlranker.util.math import set_protein_order


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[Protein]


class PeptideGroup:
    """Peptide group that can contain multiple ProteinPairs and PeptidePairs"""

    pass


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


class XLDataSet:
    def __init__(self):
        pass
