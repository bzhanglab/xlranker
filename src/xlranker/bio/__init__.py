from .protein import Protein
from .peptide import Peptide
from .protein_pair import ProteinPair
from .peptide_pair import PeptidePair

from enum import Enum, auto

__all__ = ["Protein", "Peptide", "ProteinPair", "PeptidePair", "PeptidePair"]


class PrioritizationStatus(Enum):
    NOT_ANALYZED = auto()  # No analysis performed yet

    # Parsimony-based statuses
    PARSIMONY_NOT_SELECTED = (
        auto()
    )  # Another entity was selected, or cannot be selected
    PARSIMONY_SELECTED = auto()  # Selected as the representative
    PARSIMONY_AMBIGUOUS = auto()  # No clear candidate

    # Machine Learning-based statuses
    ML_NOT_SELECTED = auto()  # Lower ML score
    ML_SELECTED = auto()  # Highest ML score


class GroupedEntity:
    group_id: int
    in_group: bool
    status: PrioritizationStatus

    def __init__(self):
        self.in_group = False
        self.group_id = -1
        self.status = PrioritizationStatus.NOT_ANALYZED

    def set_group(self, group_id: int) -> None:
        self.in_group = True
        self.group_id = group_id

    def get_group(self) -> int:
        return self.group_id

    def set_status(self, status: PrioritizationStatus) -> None:
        self.status = status
