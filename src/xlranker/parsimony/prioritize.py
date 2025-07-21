import logging
import random
from dataclasses import dataclass

from xlranker.bio.pairs import PeptidePair, ProteinPair
from xlranker.lib import XLDataSet
from xlranker.status import PrioritizationStatus, ReportStatus
from xlranker.util import get_pair_id

logger = logging.getLogger(__name__)


def select_random(
    data_set: XLDataSet,
) -> None:  # FIXME: This needs to be updated to handle the selection process
    """Resolve ambiguous groups by selected a random pair.

    This is normally done by the machine learning model. However, if training data is not available, use this function to resolve the remaining ambiguity.

    Args:
        data_set (XLDataSet): data set to resolve ambiguity

    """
    ambiguity: dict[str, list[ProteinPair]] = {}
    for pair in data_set.protein_pairs.values():
        if pair.prioritization_status != PrioritizationStatus.PARSIMONY_AMBIGUOUS:
            continue  # No ambiguity
        conn_id = pair.connectivity_id()
        if conn_id not in ambiguity:
            ambiguity[conn_id] = []
        ambiguity[conn_id].append(pair)
    for conn_id in ambiguity:
        selected_location = random.randrange(len(ambiguity[conn_id]))
        for i in range(len(ambiguity[conn_id])):
            if selected_location == i:
                ambiguity[conn_id][i].set_prioritization_status(
                    PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
                )
            else:
                ambiguity[conn_id][i].set_prioritization_status(
                    PrioritizationStatus.PARSIMONY_NOT_SELECTED
                )


@dataclass
class ParsimonyGroup:
    protein_pairs: list[ProteinPair]
    peptide_pairs: list[PeptidePair]


class ParsimonySelector:
    data_set: XLDataSet
    protein_groups: dict[int, list[ProteinPair]]
    peptide_groups: dict[int, list[PeptidePair]]
    can_prioritize: bool

    def __init__(self, data_set: XLDataSet):
        """Initialize the ParsimonySelector object

        Args:
            data_set (XLDataSet): cross-linking dataset
        """
        self.data_set = data_set
        self.protein_groups = {}
        self.peptide_groups = {}
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
        if group_id not in self.protein_groups:
            self.protein_groups[group_id] = []
        self.protein_groups[group_id].append(protein_pair)
        for peptide_pair_id in protein_pair.connections:
            self.assign_peptide_pair(
                self.data_set.peptide_pairs[peptide_pair_id], group_id
            )

    def assign_peptide_pair(self, peptide_pair: PeptidePair, group_id: int) -> None:
        if peptide_pair.in_group:
            if peptide_pair.group_id != group_id:
                raise ValueError(
                    f"Protein Pair {peptide_pair} has incorrect group id. Expected {group_id}. Got {peptide_pair.group_id}."
                )
            return
        peptide_pair.set_group(group_id)
        if group_id not in self.peptide_groups:
            self.peptide_groups[group_id] = []
        self.peptide_groups[group_id].append(peptide_pair)
        for protein_pair_id in peptide_pair.connections:
            self.assign_protein_pair(
                self.data_set.protein_pairs[protein_pair_id], group_id
            )

    def create_groups(self) -> None:
        next_group_id = 1
        for pair in self.data_set.peptide_pairs.values():
            if pair.in_group or len(pair.connections) == 0:
                continue
            self.assign_peptide_pair(pair, next_group_id)
            next_group_id += 1
        self.can_prioritize = True

    def prioritize_group(self, group_id: int) -> None:
        peptide_names = set(
            [get_pair_id(pep.a, pep.b) for pep in self.peptide_groups[group_id]]
        )
        proteins = set(self.protein_groups[group_id])
        protein_pair_groups: dict[str, list[ProteinPair]] = {}
        for protein_pair in proteins:
            conn_id = protein_pair.connectivity_id()
            if conn_id not in protein_pair_groups:
                protein_pair_groups[conn_id] = []
            protein_pair_groups[conn_id].append(protein_pair)
        while len(peptide_names) > 0:
            max_connections = 0
            best_pairs: set[str] = (
                set()
            )  # set so there is no bias towards larger groups
            for conn_id in protein_pair_groups:
                n_conn = protein_pair_groups[conn_id][0].overlap(peptide_names)
                if max_connections < n_conn:
                    max_connections = n_conn
                    best_pairs = set()
                    best_pairs.add(conn_id)
                elif max_connections == n_conn:
                    best_pairs.add(conn_id)
            selected_index = random.randint(
                0, len(best_pairs) - 1
            )  # select one random index to move forward
            best_pair_group = protein_pair_groups[list(best_pairs)[selected_index]]
            peptide_names.difference_update(best_pair_group[0].connections)
            intra_pairs: list[ProteinPair] = []
            # for pair in best_pair_group
            status = (
                PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
                if len(best_pair_group) == 1
                else PrioritizationStatus.PARSIMONY_AMBIGUOUS
            )
            for pair in best_pair_group:
                if pair.is_intra and len(best_pair_group) > 1:
                    intra_pairs.append(pair)
                else:
                    pair.set_prioritization_status(status)
                    if len(best_pair_group) == 1:
                        pair.set_score(1.01)
                        pair.set_report_status(
                            ReportStatus.CONSERVATIVE
                        )  # Unambiguous pairs are reported as CONSERVATIVE
            if len(intra_pairs) > 0:
                intra_pairs.sort(
                    key=lambda pair: (
                        -pair.a.abundance()  # type: ignore
                        if pair.a.abundance() is not None
                        else float("-inf"),
                        pair.a.name,
                    )
                )
                intra_pairs[0].set_prioritization_status(
                    PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
                )
                intra_pairs[0].set_score(1.0 + 0.01 * len(intra_pairs))
                intra_pairs[0].set_report_status(
                    ReportStatus.MINIMAL
                )  # Selected intra pairs are reported at least as MINIMAL
                for i in range(1, len(intra_pairs)):
                    intra_pairs[i].set_prioritization_status(
                        PrioritizationStatus.PARSIMONY_NOT_SELECTED  # Not selected by default
                    )
                    intra_pairs[i].set_report_status(
                        ReportStatus.EXPANDED
                    )  # Secondary intra pairs are reported as EXPANDED
                    intra_pairs[i].set_score(1.0 + 0.01 * (len(intra_pairs) - i))
        for pair_group in protein_pair_groups.values():
            for protein_pair in pair_group:
                if (
                    protein_pair.prioritization_status
                    == PrioritizationStatus.NOT_ANALYZED
                ):  # if not analyzed then pair is not selected
                    protein_pair.set_prioritization_status(
                        PrioritizationStatus.PARSIMONY_NOT_SELECTED
                    )
                    protein_pair.set_report_status(ReportStatus.ALL)

    def prioritize(self) -> None:
        if not self.can_prioritize:
            logger.warning(
                "Parsimony group creation not performed before prioritization. Running now."
            )
            self.create_groups()
        for group in self.protein_groups:
            self.prioritize_group(group)

    def run(self) -> None:
        self.create_groups()
        self.prioritize()
