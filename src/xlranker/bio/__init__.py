from .protein import Protein
from .peptide import Peptide
from .protein_pair import ProteinPair
from .peptide_pair import PeptidePair

__all__ = ["Protein", "Peptide", "ProteinPair", "PeptidePair", "PeptidePair"]


class GroupedEntity:
    group_id: int
    in_group: int

    def __init__(self):
        self.in_group = False
        self.group_id = -1

    def set_group(self, group_id: int) -> None:
        self.in_group = True
        self.group_id = group_id

    def get_group(self) -> int:
        if self.in_group:
            return self.group_id
        return -1
