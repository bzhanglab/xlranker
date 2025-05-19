from xlranker.bio import GroupedEntity
from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein, sort_proteins


class ProteinPair(GroupedEntity):
    """ProteinPair class that tracks the required data for the pipeline"""

    a: Protein
    b: Protein
    score: float | None
    is_selected: bool
    pair_name: str
    linked_peptide_pairs: list["PeptidePair"]

    def __init__(self, protein_a: Protein, protein_b: Protein):
        """Initialize a ProteinPair object, making sure a is the higher abundant protein. Input order does not matter.

        Args:
            protein_a (Protein): First protein in the pair
            protein_b (Protein): Second protein in the pair
        """
        super().__init__()
        (a, b) = sort_proteins(protein_a, protein_b)
        self.a = a
        self.b = b
        self.score = None
        self.is_selected = False
        self.linked_peptide_pairs = []
        self.pair_name = f"{self.a.name}-{self.b.name}"

    def link_pair(self, pair: "PeptidePair") -> None:
        self.linked_peptide_pairs.append(pair)

    def set_score(self, score: float):
        """Set the score of the protein pair

        Args:
            score (float): float of the score given to the pair
        """
        self.score = score

    def select(self):
        """Set this pair to be selected"""
        self.is_selected = True

    def __eq__(self, value) -> bool:
        """Checks if ProteinPairs are equivalent, without caring for order

        Args:
            value (Self): protein pair to compare to

        Returns:
            bool: True if protein pairs are equivalent, regardless of a and b order
        """
        if value.__class__ != self.__class__:
            return False
        if self.a == value.a:
            return self.b == value.b
        elif self.a == value.b:
            return self.b == value.a
        return False


class PeptidePair(GroupedEntity):
    """Peptide group that can contain multiple ProteinPairs and PeptidePairs"""

    a: Peptide
    b: Peptide
    protein_pairs: list[ProteinPair]

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide):
        super().__init__()
        self.a = peptide_a
        self.b = peptide_b
