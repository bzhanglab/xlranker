from xlranker.bio import PeptidePair, ProteinPair
from xlranker.lib import XLDataSet
import logging

logger = logging.getLogger(__name__)


class ParsimonyGroup:
    protein_pairs: list[ProteinPair]
    peptide_pairs: list[PeptidePair]


class ParsimonySelector:
    data_set: XLDataSet

    def __init__(self, data_set: XLDataSet):
        self.data_set = data_set

    def prioritize_protein_pair(self, protein_pair: ProteinPair, group_id: int) -> None:
        if protein_pair.in_group:
            if protein_pair.group_id != group_id:
                raise ValueError(
                    f"Protein Pair {protein_pair} has incorrect group id. Expected {group_id}. Got {protein_pair.group_id}."
                )
            return
        protein_pair.set_group(group_id)
        for peptide_pair in protein_pair.linked_peptide_pairs:
            self.prioritize_peptide_pair(peptide_pair, group_id)

    def prioritize_peptide_pair(self, peptide_pair: PeptidePair, group_id: int) -> None:
        if peptide_pair.in_group:
            if peptide_pair.group_id != group_id:
                raise ValueError(
                    f"Protein Pair {peptide_pair} has incorrect group id. Expected {group_id}. Got {peptide_pair.group_id}."
                )
            return
        peptide_pair.set_group(group_id)
        for protein_pair in peptide_pair.protein_pairs:
            self.prioritize_protein_pair(protein_pair, group_id)

    def create_groups(self) -> None:
        next_group_id = 1
        for pair in self.data_set.network:
            if pair.in_group:
                continue
            self.prioritize_peptide_pair(pair, next_group_id)
            next_group_id += 1
