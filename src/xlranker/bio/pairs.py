from enum import Enum, auto
from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein, sort_proteins
from xlranker.util import get_pair_id, safe_a_greater_or_equal_to_b


class PrioritizationStatus(Enum):
    """Prioritization status for a protein pair"""

    NOT_ANALYZED = auto()
    "No analysis performed yet"

    # Parsimony-based statuses

    PARSIMONY_NOT_SELECTED = auto()
    "Another entity was selected in group or cannot be selected"
    PARSIMONY_PRIMARY_SELECTED = auto()
    "Selected as the primary representative for group."
    PARSIMONY_SECONDARY_SELECTED = auto()
    "Selected as a secondary representative for group. Only possible for intra pairs."
    PARSIMONY_AMBIGUOUS = auto()
    "No clear candidate from parsimony analysis. Needs ML model."

    # Machine Learning-based statuses

    ML_NOT_SELECTED = auto()
    "Other candidate had higher score in group"
    ML_PRIMARY_SELECTED = auto()
    "Highest ML score in group or primary selection"
    ML_SECONDARY_SELECTED = auto()
    "High confidence ML score in group or secondary selection"


class GroupedEntity:
    group_id: int
    subgroup_id: int
    in_group: bool
    status: PrioritizationStatus
    connections: set[str]

    def __init__(self) -> None:
        self.in_group = False
        self.group_id = -1
        self.subgroup_id = 0
        self.status = PrioritizationStatus.NOT_ANALYZED
        self.connections = set()

    def set_group(self, group_id: int) -> None:
        self.in_group = True
        self.group_id = group_id

    def set_subgroup(self, subgroup_id: int) -> None:
        self.subgroup_id = subgroup_id

    def get_group_string(self) -> str:
        return f"{self.group_id}.{self.subgroup_id}"

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
    score: float
    is_selected: bool
    pair_id: str
    is_intra: bool

    def __init__(self, protein_a: Protein, protein_b: Protein) -> None:
        """Initialize a ProteinPair object, making sure a is the higher abundant protein. Input order does not matter.

        Args:
            protein_a (Protein): First protein in the pair
            protein_b (Protein): Second protein in the pair
        """
        super().__init__()
        (a, b) = sort_proteins(protein_a, protein_b)
        self.a = a
        self.b = b
        self.score = -1
        self.is_selected = False
        self.pair_id = get_pair_id(a, b)
        self.is_intra = a == b

    def set_score(self, score: float) -> None:
        """Set the score of the protein pair.

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

    def abundance_dict(self) -> dict[str, str | float | None]:
        """convert ProteinPair into dictionary of abundances, making abundances ending in a being the larger value

        Returns:
            dict[str, float | None]: dictionary where keys are the abundance name and the values being the abundance value
        """
        ret_val: dict[str, str | float | None] = {"pair": self.pair_id}
        for abundance_key in self.a.abundances:
            a = self.a.abundances[abundance_key]
            b = self.b.abundances[abundance_key]
            if safe_a_greater_or_equal_to_b(a, b):
                ret_val[f"{abundance_key}_a"] = a
                ret_val[f"{abundance_key}_b"] = b
            else:  # make sure a is the larger value
                ret_val[f"{abundance_key}_a"] = b
                ret_val[f"{abundance_key}_b"] = a
        return ret_val

    def to_tsv(self) -> str:
        """converts object into a TSV string

        Returns:
            str: TSV representation of the protein pair, including id and status
        """
        return f"{self.pair_id}\t{self.status}\t{self.get_group_string()}"

    def __hash__(self) -> int:
        return hash(self.pair_id)


class PeptidePair(GroupedEntity):
    """Peptide group that can contain multiple ProteinPairs and PeptidePairs."""

    a: Peptide
    b: Peptide
    pair_id: str

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide) -> None:
        super().__init__()
        self.a = peptide_a
        self.b = peptide_b
        self.pair_id = get_pair_id(peptide_a, peptide_b)

    def __hash__(self) -> int:
        return hash(self.pair_id)
