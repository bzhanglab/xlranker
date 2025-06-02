"""
Model Process:

1. Identify Positive Dataset
    - All representative pairs from parsimonious selection
2. Generate Negative Dataset
    - Random protein pairs that are not candidate pairs
"""

import logging
import random
import polars as pl
from xlranker.bio.pairs import ProteinPair, PrioritizationStatus
from xlranker.bio.protein import Protein
from xlranker.lib import XLDataSet
from xlranker import config

logger = logging.getLogger(__name__)


class ModelConfig:
    runs: int
    folds: int

    def __init__(self, runs: int = 10, folds: int = 5):
        self.runs = runs
        self.folds = folds

    def validate(self) -> bool:
        attrs = {
            "runs": (int, lambda x: x >= 1),
            "folds": (int, lambda x: x >= 1),
        }
        for attr, (typ, cond) in attrs.items():
            value = getattr(self, attr, None)
            if not isinstance(value, typ):
                return False
            if cond and not cond(value):
                return False
        return True

    def as_dict(self) -> dict:
        # Utility to get all config as a dict
        return {k: getattr(self, k) for k in self.__dict__}


class PrioritizationModel:
    positives: list[ProteinPair]
    dataset: XLDataSet
    existing_pairs: set[tuple[Protein, Protein]]
    model_config: ModelConfig
    na_count: int

    def __init__(self, dataset: XLDataSet, model_config: ModelConfig = ModelConfig()):
        self.dataset = dataset
        self.positives = []
        self.existing_pairs: set[tuple[Protein, Protein]] = set(
            (p.a, p.b) for p in self.dataset.protein_pairs.values()
        )
        self.model_config = model_config

    def get_positives(self):
        for protein_pair in self.dataset.protein_pairs.values():
            if protein_pair.status == PrioritizationStatus.PARSIMONY_SELECTED:
                self.positives.append(protein_pair)

    def get_negatives(self, n: int) -> list[ProteinPair]:
        """get a list of negative protein pairs

        Args:
            n (int): the number of pairs to generate

        Returns:
            list[ProteinPair]: list of negative protein pairs
        """
        negatives: list[ProteinPair] = []
        n_prot = len(self.dataset.proteins.values())
        if n > (n_prot * (n_prot - 1)) // 2 - len(self.positives):
            msg = f"n value for get_negatives ({n}) is too large. Setting to maximum value: {(n_prot * (n_prot - 1)) // 2 - len(self.positives)}"
            if config.fragile:
                logger.error(msg)
                raise ValueError(
                    "get_negatives(n: int) n value is too large and fragile is True"
                )
            logger.warning(msg)
            n = (n_prot * (n_prot - 1)) // 2 - len(self.positives)
        protein_ids = list(self.dataset.proteins.keys())

        generated = set()

        while len(negatives) < n:
            a, b = random.sample(protein_ids, 2)
            pair_key = tuple(sorted([a, b]))
            if pair_key in self.existing_pairs or pair_key in generated:
                continue
            negatives.append(
                ProteinPair(self.dataset.proteins[a], self.dataset.proteins[b])
            )
            generated.add(pair_key)
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

    def get_predictions(self, df: pl.DataFrame):
        pass

    def run(self):
        pass
