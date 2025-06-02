"""
Model Process:

1. Identify Positive Dataset
    - All representative pairs from parsimonious selection
2. Generate Negative Dataset
    - Random protein pairs that are not selected pairs
"""

import logging
import polars as pl
from xlranker.bio.pairs import ProteinPair, PrioritizationStatus
from xlranker.lib import XLDataSet

logger = logging.getLogger(__name__)


class PrioritizationModel:
    positives: list[ProteinPair]
    dataset: XLDataSet

    def __init__(self, dataset: XLDataSet):
        self.dataset = dataset
        self.positives = []

    def get_positives(self):
        for protein_pair in self.dataset.protein_pairs.values():
            if protein_pair.status == PrioritizationStatus.PARSIMONY_SELECTED:
                self.positives.append(protein_pair)

    def get_negatives(self, n: int) -> list[ProteinPair]:
        negatives: list[ProteinPair] = []
        n_prot = len(self.dataset.proteins.values())
        if n > (n_prot * (n_prot - 1)) // 2 - len(self.positives):
            logger.warning(
                f"n value for get_negatives ({n}) is too large. Setting to maximum value: {(n_prot * (n_prot - 1)) // 2 - len(self.positives)}"
            )
            n = (n_prot * (n_prot - 1)) // 2 - len(self.positives)
        return negatives

    def construct_df(self, negative_pairs: list[ProteinPair]) -> pl.DataFrame:
        """generate a Polars DataFrame from the positive pairs and a list of negative ProteinPair

        Args:
            negative_pairs (list[ProteinPair]): the list of negative pairs to add to DataFrame

        Returns:
            pl.DataFrame: DataFrame where the first column is 'pair', followed by abundances. Last column is 'label'
        """
        df_array: list[dict[str, str | int | float | None]] = []
        headers = ["pair"]  # headers in the correct order
        for pair in self.positives:
            pair_dict = pair.abundance_dict()
            pair_dict["label"] = 1
        headers.extend(
            [k for k in self.positives[0].abundance_dict().keys() if k != "pair"]
        )
        headers.append("label")
        for pair in negative_pairs:
            pair_dict = pair.abundance_dict()
            pair_dict["label"] = 0
        return pl.DataFrame(df_array).select(headers)
