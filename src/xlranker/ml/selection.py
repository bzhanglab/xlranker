from abc import ABC, abstractmethod

from xlranker.bio.pairs import PrioritizationStatus, ProteinPair


class PairSelector(ABC):
    @abstractmethod
    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def process(self, protein_pairs: list[ProteinPair]) -> None:
        pass

    def assign_subgroups_and_get_best(
        self, protein_pairs: list[ProteinPair]
    ) -> dict[str, float]:
        """assign subgroups to pairs and get the best score for each subgroup

        Returns:
            dict[str, float]: dict where key is the connectivity ID (str) and the values are the highest score (float)
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
    with_secondary: bool

    def __init__(self, with_secondary: bool = False) -> None:
        """Selects the pair with the highest score. For pairs tied for best score, pair alphabetically first is selected. If with_secondary is True, assign pairs with best score but not alphabetically first as secondary selections.

        Args:
            with_secondary (bool, optional): Flag to keep tied pairs that are not the primary selection. Defaults to False.
        """
        super().__init__()
        self.with_secondary = with_secondary

    def process(self, protein_pairs: list[ProteinPair]) -> None:
        best_score = self.assign_subgroups_and_get_best(protein_pairs)
        best_pair: dict[str, ProteinPair] = {}
        replaced_status = (
            PrioritizationStatus.ML_SECONDARY_SELECTED
            if self.with_secondary
            else PrioritizationStatus.ML_NOT_SELECTED
        )  # get status for scores with best score but not alphabetically first
        for pair in protein_pairs:
            if pair.score == best_score[pair.connectivity_id()]:
                if pair.connectivity_id() not in best_pair:
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
                elif (
                    pair.pair_id < best_pair[pair.connectivity_id()].pair_id
                ):  # alphabetically sort
                    best_pair[
                        pair.connectivity_id()
                    ].status = replaced_status  # replace previous best
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
            else:
                pair.status = PrioritizationStatus.ML_NOT_SELECTED


class ThresholdSelector(PairSelector):
    threshold: float

    def __init__(self, threshold: float) -> None:
        super().__init__()
        self.threshold = threshold

    def process(self, protein_pairs: list[ProteinPair]) -> None:
        best_score = self.assign_subgroups_and_get_best(protein_pairs)
        best_pair: dict[str, ProteinPair] = {}
        for pair in protein_pairs:
            if pair.score == best_score[pair.connectivity_id()]:
                if pair.connectivity_id() not in best_pair:
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
                elif (
                    pair.pair_id < best_pair[pair.connectivity_id()].pair_id
                ):  # alphabetically sort
                    best_pair[pair.connectivity_id()].status = (
                        PrioritizationStatus.ML_SECONDARY_SELECTED
                        if pair.score >= self.threshold
                        else PrioritizationStatus.ML_NOT_SELECTED
                    )  # replace previous best, but keep as secondary if above threshold
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
            elif pair.score >= self.threshold:  # select if above threshold
                pair.status = PrioritizationStatus.ML_SECONDARY_SELECTED
            else:
                pair.status = PrioritizationStatus.ML_NOT_SELECTED


class WithinBestScoreSelector(PairSelector):
    top_n: int
    within: float

    def __init__(self, top_n: int, within: float) -> None:
        super().__init__()
        if within > 1 or within < 0:  # TODO: Decide if necessary
            raise ValueError(
                "within must be between 0 and 1 for WithinBestScoreSelector!"
            )
        self.top_n = top_n
        self.within = within

    def process(self, protein_pairs: list[ProteinPair]) -> None:
        best_score = self.assign_subgroups_and_get_best(protein_pairs)
        best_pair: dict[str, ProteinPair] = {}
        subgroups: dict[int, list[ProteinPair]] = {}
        for pair in protein_pairs:
            if pair.score == best_score[pair.connectivity_id()]:
                if pair.connectivity_id() not in best_pair:
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
                elif (
                    pair.pair_id < best_pair[pair.connectivity_id()].pair_id
                ):  # alphabetically sort
                    best_pair[
                        pair.connectivity_id()
                    ].status = (
                        PrioritizationStatus.ML_NOT_SELECTED
                    )  # replace previous best
                    pair.status = PrioritizationStatus.ML_PRIMARY_SELECTED
                    best_pair[pair.connectivity_id()] = pair
            else:
                pair.status = PrioritizationStatus.ML_NOT_SELECTED
        for pair in protein_pairs:
            subgroup = pair.subgroup_id
            conn_id = pair.connectivity_id()
            if subgroup not in subgroups:
                subgroups[subgroup] = []
            if (
                best_score[conn_id] * (1.0 - self.within) <= pair.score
                and pair.status != PrioritizationStatus.ML_PRIMARY_SELECTED
            ):  # Check if within
                subgroups[subgroup].append(pair)
        for subgroup in subgroups:
            group_list = subgroups[subgroup]
            if len(group_list) <= self.top_n:  # no need to sort
                for pair in group_list:
                    pair.status = PrioritizationStatus.ML_SECONDARY_SELECTED
            else:
                group_list.sort(
                    key=lambda pair: (-pair.score, pair.pair_id)
                )  # -pair.score makes it so higher scores come first
                for i in range(self.top_n):
                    group_list[i].status = PrioritizationStatus.ML_SECONDARY_SELECTED
