"""Methods and classes for the final selection of protein pairs after Parsimony and ML processing."""

from abc import ABC, abstractmethod

from xlranker.bio.pairs import ProteinPair
from xlranker.status import PrioritizationStatus, ReportStatus


def filter_for_undecided_pairs(protein_pairs: list[ProteinPair]) -> list[ProteinPair]:
    """Only get pairs that don't have any selected status.

    Args:
        protein_pairs (list[ProteinPair]): list of protein pairs to filter.

    Returns:
        list[ProteinPair]: list of all protein pairs still needing a selection status.

    """
    return [
        pair
        for pair in protein_pairs
        if not (
            pair.prioritization_status == PrioritizationStatus.PARSIMONY_NOT_SELECTED
            and pair.score == -1.0
        )
        and pair.prioritization_status
        != PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
    ]


def assign_not_selected_status(protein_pair: ProteinPair) -> None:
    """Assign the correct not selected status to protein pair.

    Args:
        protein_pair (ProteinPair): protein pair needing status assignment.

    """
    if protein_pair.score > 1.0:
        protein_pair.set_prioritization_status(
            PrioritizationStatus.PARSIMONY_NOT_SELECTED
        )
    else:
        protein_pair.set_prioritization_status(PrioritizationStatus.ML_NOT_SELECTED)


def assign_secondary_selected_status(protein_pair: ProteinPair) -> None:
    """Assign the correct secondary selected status to protein pair.

    Args:
        protein_pair (ProteinPair): protein pair needing status assignment.

    """
    if protein_pair.score >= 1.01:
        protein_pair.set_prioritization_status(
            PrioritizationStatus.PARSIMONY_SECONDARY_SELECTED
        )
    else:
        protein_pair.set_prioritization_status(
            PrioritizationStatus.ML_SECONDARY_SELECTED
        )
    protein_pair.set_report_status(
        ReportStatus.EXPANDED
    )  # Secondary selections are reported as EXPANDED


def assign_primary_selected_status(protein_pair: ProteinPair) -> None:
    """Assign the correct primary selected status to protein pair.

    Args:
        protein_pair (ProteinPair): protein pair needing status assignment.

    """
    if protein_pair.score > 1.0:
        protein_pair.set_prioritization_status(
            PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
        )  # No change to report status for parsimony primary selections, done in parsimony step
    else:
        protein_pair.set_prioritization_status(PrioritizationStatus.ML_PRIMARY_SELECTED)
        protein_pair.set_report_status(
            ReportStatus.MINIMAL
        )  # Primary ML selections are reported as MINIMAL


class PairSelector(ABC):
    """Abstract class describing methods for a protein pair selector."""

    @abstractmethod
    def __init__(self) -> None:
        """Initialize the pair selector."""
        super().__init__()

    @abstractmethod
    def process(self, protein_pairs: list[ProteinPair]) -> None:
        """Process and assign all protein pairs a status.

        Args:
            protein_pairs (list[ProteinPair]): protein pairs to process.

        """
        pass

    def assign_subgroups_and_get_best(
        self, protein_pairs: list[ProteinPair]
    ) -> dict[str, float]:
        """Assign subgroups to pairs and get the best score for each subgroup.

        Returns:
            dict[str, float]: dict where key is the connectivity ID (str)
                              and the values are the highest score (float)

        """
        best_score: dict[str, float] = {}
        subgroups: dict[str, int] = {}
        subgroup_id = 1
        for pair in protein_pairs:
            conn_id = pair.connectivity_id()
            if conn_id not in best_score:
                best_score[conn_id] = pair.score
                subgroups[conn_id] = subgroup_id
                subgroup_id += 1
            elif best_score[conn_id] < pair.score:
                best_score[conn_id] = pair.score
            pair.set_subgroup(subgroups[conn_id])
        return best_score


class BestSelector(PairSelector):
    """PairSelector that only keeps the best score. Optionally can allow secondary selections for ties.

    Attributes:
        with_secondary (bool): if True, assign secondary selections to tied protein pairs.

    """

    with_secondary: bool

    def __init__(self, with_secondary: bool = False) -> None:
        """Select the pair with the highest score.

        For pairs tied for best score, pair alphabetically first is selected.
        If with_secondary is True, assign pairs with best score but not alphabetically
        first as secondary selections.

        Args:
            with_secondary (bool, optional): Flag to keep tied pairs. Defaults to False.

        """
        super().__init__()
        self.with_secondary = with_secondary

    def process(self, protein_pairs: list[ProteinPair]) -> None:
        """Process and assign all protein pairs a status.

        Args:
            protein_pairs (list[ProteinPair]): protein pairs to process.

        """
        best_score = self.assign_subgroups_and_get_best(protein_pairs)
        best_pair: dict[str, ProteinPair] = {}
        replaced_status = (
            PrioritizationStatus.ML_SECONDARY_SELECTED
            if self.with_secondary
            else PrioritizationStatus.ML_NOT_SELECTED
        )  # get status for scores with best score but not alphabetically first
        for pair in protein_pairs:
            if (
                pair.score == best_score[pair.connectivity_id()]
            ):  # NOTE: Multiple pairs with best score only possible if all pairs are inter pairs
                if pair.connectivity_id() not in best_pair:
                    assign_primary_selected_status(pair)
                    best_pair[pair.connectivity_id()] = pair
                elif (
                    pair.pair_id < best_pair[pair.connectivity_id()].pair_id
                ):  # alphabetically sort, only triggered by inter pairs
                    best_pair[
                        pair.connectivity_id()
                    ].prioritization_status = replaced_status  # replace previous best
                    assign_primary_selected_status(pair)
                    best_pair[pair.connectivity_id()] = pair
            else:
                assign_not_selected_status(pair)


class ThresholdSelector(PairSelector):
    """PairSelector that selects all pairs that pass a scoring threshold.

    Attributes:
        threshold (float): minimum scoring threshold for pair to be selected.
        top_n (int | None): if not None, only keep top_n pairs in a group above threshold.
    """

    threshold: float
    top_n: int | None

    def __init__(self, threshold: float, top_n: int | None = None) -> None:
        """Initialize the ThresholdSelector.

        Args:
            threshold (float): minimum scoring threshold for pair to be selected.
            top_n (int | None, optional): if not None, only keep top_n pairs in a group above threshold. Defaults to None.

        """
        super().__init__()
        self.threshold = threshold
        self.top_n = top_n

    def process(self, protein_pairs: list[ProteinPair]) -> None:
        """Process and assign all protein pairs a status.

        Args:
            protein_pairs (list[ProteinPair]): protein pairs to process.

        """
        best_score = self.assign_subgroups_and_get_best(protein_pairs)
        best_pair: dict[str, ProteinPair] = {}
        subgroups: dict[int, list[ProteinPair]] = {}
        protein_pairs = filter_for_undecided_pairs(protein_pairs)
        for pair in protein_pairs:
            conn_id = pair.connectivity_id()
            if pair.score == best_score[conn_id]:
                if conn_id not in best_pair:
                    assign_primary_selected_status(pair)
                    best_pair[conn_id] = pair
                elif (
                    pair.pair_id < best_pair[conn_id].pair_id
                ):  # alphabetically sort, only possible if all inter
                    best_pair[
                        conn_id
                    ].prioritization_status = PrioritizationStatus.ML_NOT_SELECTED
                    best_pair[conn_id].set_report_status(ReportStatus.ALL)
                    pair.prioritization_status = (
                        PrioritizationStatus.ML_PRIMARY_SELECTED
                    )
                    pair.set_report_status(ReportStatus.MINIMAL)
                    best_pair[conn_id] = pair
            else:
                assign_not_selected_status(pair)
        for pair in protein_pairs:
            subgroup = pair.subgroup_id
            conn_id = pair.connectivity_id()
            if subgroup not in subgroups:
                subgroups[subgroup] = []
            if (
                pair.score > self.threshold  # greater than threshold
                and pair.prioritization_status
                != PrioritizationStatus.ML_PRIMARY_SELECTED
                and pair.prioritization_status
                != PrioritizationStatus.PARSIMONY_PRIMARY_SELECTED
            ):  # Check if within
                subgroups[subgroup].append(pair)
        for subgroup in subgroups:
            group_list = subgroups[subgroup]
            if self.top_n is None:  # Select all
                for pair in group_list:
                    assign_secondary_selected_status(pair)
            else:
                if len(group_list) < self.top_n:  # no need to sort
                    for pair in group_list:
                        assign_secondary_selected_status(pair)
                else:
                    group_list.sort(
                        key=lambda pair: (-pair.score, pair.pair_id)
                    )  # -pair.score makes it so higher scores come first
                    for i in range(self.top_n - 1):
                        assign_secondary_selected_status(group_list[i])
