from xlranker.bio import PrioritizationStatus
from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein, sort_proteins
from xlranker.util import get_pair_id


class GroupedEntity:
    group_id: int
    in_group: bool
    status: PrioritizationStatus
    connections: set[str]

    def __init__(self):
        self.in_group = False
        self.group_id = -1
        self.status = PrioritizationStatus.NOT_ANALYZED
        self.connections = set()

    def set_group(self, group_id: int) -> None:
        self.in_group = True
        self.group_id = group_id

    def get_group(self) -> int:
        return self.group_id

    def set_status(self, status: PrioritizationStatus) -> None:
        self.status = status

    def add_connection(self, entity: str) -> None:
        self.connections.add(entity)

    def remove_connections(self, entities: set) -> None:
        self.connections.difference_update(entities)

    def n_connections(self) -> int:
        return len(self.connections)

    def overlap(self, entities: set[str]) -> int:
        return len(self.connections.intersection(entities))

    def same_connectivity(self, grouped_entity: "GroupedEntity") -> bool:
        return (
            len(self.connections.symmetric_difference(grouped_entity.connections)) == 0
        )

    def connectivity_id(self) -> str:
        """Returns a unique, order-independent id for the set of connections."""
        return "|".join(sorted(self.connections))


class ProteinPair(GroupedEntity):
    """ProteinPair class that tracks the required data for the pipeline"""

    a: Protein
    b: Protein
    score: float | None
    is_selected: bool
    pair_id: str

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
        self.pair_id = get_pair_id(a, b)

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

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide):
        super().__init__()
        self.a = peptide_a
        self.b = peptide_b
