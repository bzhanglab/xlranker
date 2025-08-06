"""Protein and Peptide pairs."""

from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein, sort_proteins
from xlranker.status import PrioritizationStatus, ReportStatus
from xlranker.util import get_pair_id, safe_a_greater_or_equal_to_b


class GroupedEntity:
    """Entity that are part of a group with connections to other entities.

    Attributes:
        group_id (int): ID of the large connected component entity is a part of
        subgroup_id (int): ID of the unique subgroup entity is a part of
        in_group (bool): True if entity has been assigned a group
        prioritization_status (PrioritizationStatus): status of the entity during prioritization
        connections (set[str]): List of other grouped entities this entity is connected to.
                                Values are connection IDs.

    """

    group_id: int
    subgroup_id: int
    in_group: bool
    prioritization_status: PrioritizationStatus
    connections: set[str]

    def __init__(self) -> None:
        """Initialize a grouped entity."""
        self.in_group = False
        self.group_id = -1
        self.subgroup_id = 0
        self.prioritization_status = PrioritizationStatus.NOT_ANALYZED
        self.connections = set()

    def set_group(self, group_id: int) -> None:
        """Set the group this entity is a part of.

        Args:
            group_id (int): ID of the group to set this entity to.

        """
        self.in_group = True
        self.group_id = group_id

    def set_subgroup(self, subgroup_id: int) -> None:
        """Set the subgroup ID for this entity.

        Args:
            subgroup_id (int): ID of the subgroup

        """
        self.subgroup_id = subgroup_id

    def get_group_string(self) -> str:
        """Return a string representation of the group this entity is in.

        Returns:
            str: String of in the format of group_id.subgroup_id

        """
        return f"{self.group_id}.{self.subgroup_id}"

    def get_group(self) -> int:
        """Get the group id of this entity.

        Returns:
            int: group ID this entity is a part of

        """
        return self.group_id

    def set_prioritization_status(self, status: PrioritizationStatus) -> None:
        """Set the prioritization status of this entity.

        Args:
            status (PrioritizationStatus): prioritization status to assign to this entity.

        """
        self.prioritization_status = status

    def add_connection(self, entity: str) -> None:
        """Add a connection to a different entity.

        Args:
            entity (str): String representation of the entity to add

        """
        self.connections.add(entity)

    def remove_connections(self, entities: set) -> None:
        """Remove multiple connections from this entity.

        Args:
            entities (set): set of entities to remove connections to

        """
        self.connections.difference_update(entities)

    def n_connections(self) -> int:
        """Get number of connections to this entity.

        Returns:
            int: number of entity connected to this entity

        """
        return len(self.connections)

    def overlap(self, entities: set[str]) -> int:
        """Overlap of connections to a set of other entities.

        Args:
            entities (set[str]): list of entities to compare to

        Returns:
            int: number of entities in input set also connected to this entity

        """
        return len(self.connections.intersection(entities))

    def same_connectivity(self, grouped_entity: "GroupedEntity") -> bool:
        """Determine if another GroupedEntity has the same connection as this entity.

        Args:
            grouped_entity (GroupedEntity): entity to compare to

        Returns:
            bool: True if this entity and input entity have identical connections

        """
        return (
            len(self.connections.symmetric_difference(grouped_entity.connections)) == 0
        )

    def connectivity_id(self) -> str:
        """Returns a unique, ordered id for the set of connections."""
        return "|".join(sorted(self.connections))


class ProteinPair(GroupedEntity):
    """ProteinPair class that tracks the required data for the pipeline.

    Args:
        protein_a (Protein): First protein in the pair
        protein_b (Protein): Second protein in the pair

    Attributes:
        a (Protein): Protein A
        b (Protein): Protein B
        score (float): Prioritization score
        is_selected (bool): True if this pair is selected as a representative pair
        pair_id (str): str representation of this protein pair
        is_intra (bool): True if protein a and b are the same
        report_status (ReportStatus): status of what report group this pair belongs to

    """

    a: Protein
    b: Protein
    score: float
    is_selected: bool
    pair_id: str
    is_intra: bool
    report_status: ReportStatus

    def __init__(self, protein_a: Protein, protein_b: Protein) -> None:
        """Initialize the protein pair.

        Order of proteins does not matter.

        Args:
            protein_a (Protein): the first protein object
            protein_b (Protein): the second protein object

        """
        super().__init__()
        (a, b) = sort_proteins(protein_a, protein_b)
        self.a = a
        self.b = b
        self.score = -1
        self.is_selected = False
        self.pair_id = get_pair_id(a, b)
        self.is_intra = a == b
        self.report_status = ReportStatus.NONE

    def set_score(self, score: float) -> None:
        """Set the score of the protein pair.

        Args:
            score (float): float of the score given to the pair

        """
        self.score = score

    def set_report_status(self, status: ReportStatus) -> None:
        """Set the report status of the protein pair.

        Args:
            status (ReportStatus): ReportStatus enum value to set the pair to

        """
        self.report_status = status

    def select(self):
        """Set this pair to be selected."""
        self.is_selected = True

    def __eq__(self, value: "object | ProteinPair") -> bool:
        """Checks if ProteinPairs are equivalent, without caring for order.

        Args:
            value (object | ProteinPair): protein pair to compare to

        Returns:
            bool: True if protein pairs are equivalent, regardless of a and b order

        """
        if value.__class__ != self.__class__:
            return False
        if not isinstance(value, ProteinPair):
            return False
        if self.a == value.a:
            return self.b == value.b
        elif self.a == value.b:
            return self.b == value.a
        return False

    def abundance_dict(self) -> dict[str, str | float | None]:
        """Convert ProteinPair into dictionary of abundances, making abundances ending in a being the larger value.

        Returns:
            dict[str, str | float | None]: dictionary where keys are the abundance name and the values being the abundance value

        """
        ret_val: dict[str, str | float | None] = {"pair": self.pair_id}
        for abundance_key in self.a.abundances:
            a = self.a.abundances[abundance_key]
            b = self.b.abundances[abundance_key]
            ret_val[f"{abundance_key}_a"] = a
            ret_val[f"{abundance_key}_b"] = b
        return ret_val

    def to_tsv(self) -> str:
        """Converts object into a TSV string.

        Returns:
            str: TSV representation of the protein pair, including id and status

        """
        return f"{self.pair_id}\t{self.report_status}\t{self.prioritization_status}\t{self.get_group_string()}"

    def __hash__(self) -> int:
        """Generate a hash representing this protein pair.

        Returns:
            int: hash generated from the pair_id of this protein pair.

        """
        return hash(self.pair_id)


class PeptidePair(GroupedEntity):
    """Pair of two peptide sequences. Order of a and b does not matter.

    Args:
        peptide_a (Peptide): Peptide A
        peptide_b (Peptide): Peptide B

    Attributes:
        a (Peptide): Peptide a
        b (Peptide): Peptide b
        pair_id (str): String representation of this peptide pair.

    """

    a: Peptide
    b: Peptide
    pair_id: str

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide) -> None:
        """Peptide sequence pairs. Used for parsimonious selection.

        Input order does not matter.

        Args:
            peptide_a (Peptide): first peptide object
            peptide_b (Peptide): second peptide object

        """
        super().__init__()
        self.a = peptide_a
        self.b = peptide_b
        self.pair_id = get_pair_id(peptide_a, peptide_b)

    def __hash__(self) -> int:
        """Get hash of this PeptidePair.

        Returns:
            int: hash generated from pair id

        """
        return hash(self.pair_id)
