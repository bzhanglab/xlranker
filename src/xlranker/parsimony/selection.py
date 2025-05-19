from xlranker.bio import PrioritizationStatus
from xlranker.bio.pairs import PeptidePair
from xlranker.bio.pairs import ProteinPair
from xlranker.lib import XLDataSet
import networkx as nx
import logging

from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class ParsimonyGroup:
    protein_pairs: list[ProteinPair]
    peptide_pairs: list[PeptidePair]


class ParsimonySelector:
    data_set: XLDataSet
    groups: dict[int, list[ProteinPair]]
    can_prioritize: bool
    network: nx.Graph

    def __init__(self, data_set: XLDataSet):
        self.data_set = data_set
        self.groups = {}
        self.can_prioritize = False
        self.network = None

    def assign_protein_pair(self, protein_pair: ProteinPair, group_id: int) -> None:
        if protein_pair.in_group:
            if protein_pair.group_id != group_id:
                raise ValueError(
                    f"Protein Pair {protein_pair} has incorrect group id. Expected {group_id}. Got {protein_pair.group_id}."
                )
            return
        protein_pair.set_group(group_id)
        if group_id not in self.groups:
            self.groups[group_id] = []
        self.groups[group_id].append(protein_pair)
        for peptide_pair in protein_pair.linked_peptide_pairs:
            self.assign_peptide_pair(peptide_pair, group_id)

    def assign_peptide_pair(self, peptide_pair: PeptidePair, group_id: int) -> None:
        if peptide_pair.in_group:
            if peptide_pair.group_id != group_id:
                raise ValueError(
                    f"Protein Pair {peptide_pair} has incorrect group id. Expected {group_id}. Got {peptide_pair.group_id}."
                )
            return
        peptide_pair.set_group(group_id)
        for protein_pair in peptide_pair.protein_pairs:
            self.assign_protein_pair(protein_pair, group_id)

    def create_groups(self) -> None:
        next_group_id = 1
        for pair in self.data_set.network:
            if pair.in_group:
                continue
            self.assign_peptide_pair(pair, next_group_id)
            next_group_id += 1
        self.can_prioritize = True

    def build_network(self) -> None:
        g = nx.Graph()
        g.add_nodes_from(self.data_set.network)
        g.add_nodes_from(self.data_set.network)

    def prioritize(self) -> None:
        if not self.can_prioritize:
            logger.warning(
                "Parsimony group creation not performed before prioritization. Running now."
            )
            self.create_groups()
        for group in self.groups:
            best_pairs: list[ProteinPair] = []
            max_connections = 0
            for protein_pair in self.groups[group]:
                conn_len = len(protein_pair.linked_peptide_pairs)
                if max_connections == conn_len:
                    best_pairs.append(protein_pair)
                elif max_connections < conn_len:
                    best_pairs = [protein_pair]
                    max_connections = conn_len
            best_status = (
                PrioritizationStatus.PARSIMONY_SELECTED
                if len(best_pairs) == 1
                else PrioritizationStatus.PARSIMONY_AMBIGUOUS
            )
            for protein_pair in self.groups[group]:
                if protein_pair in best_pairs:
                    protein_pair.set_status(best_status)
                else:
                    protein_pair.set_status(PrioritizationStatus.PARSIMONY_NOT_SELECTED)

    def run(self) -> None:
        self.create_groups()
        self.prioritize()
