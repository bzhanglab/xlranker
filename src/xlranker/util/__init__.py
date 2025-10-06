"""Pipeline utility functions and classes."""

import random

import numpy as np
import polars as pl

from xlranker.config import config
from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein


def set_seed(seed: int) -> None:
    """Set seed to provide consistent results between runs.

    Args:
        seed (int): number to initialize random number generators with

    """
    random.seed(seed)
    np.random.seed(int(random.random() * 1000000))


def get_pair_id(a: Protein | Peptide, b: Protein | Peptide) -> str:
    """Get a string representation of the pair. Input order independent.

    Order is determine alphabetically.

    Args:
        a (Protein | Peptide): entity a
        b (Protein | Peptide): entity b

    Returns:
        str: pair representation with entities separated by config.advanced.pair_separator (default '+').

    """
    name_a = ""
    name_b = ""
    if isinstance(a, Protein):
        name_a = a.name
    else:
        name_a = a.sequence
    if isinstance(b, Protein):
        name_b = b.name
    else:
        name_b = b.sequence
    if name_a < name_b:
        return f"{name_a}{config.advanced.pair_separator}{name_b}"
    return f"{name_b}{config.advanced.pair_separator}{name_a}"


def safe_a_greater_or_equal_to_b(a: float | None, b: float | None) -> bool:
    """Returns True if a is greater or equal to b, with checks for None.

    None is treated as missing value. Any float is greater than None. If both are None, return True.

    Args:
        a (float | None): a value
        b (float | None): b value

    Returns:
        bool: True if a is greater or equal to b. If both are None, return True. Any float is greater than None.

    """
    if a is None:
        return b is None  # if a is None, then if b is not None, b is greater
    else:
        if b is None:
            return True  # Non-None is always greater than None
        return a >= b  # both are not None, so compare normally


def get_abundance(
    omic_df: pl.DataFrame, analyte: str, use_median: bool = False
) -> float | None:
    """Get the mean or median abundance of an analyte from an omics dataset.

    Args:
        omic_df (pl.DataFrame): Polars dataframe containing the omics data, with the first column being the index.
        analyte (str): analyte that should have an exact match in omic_df.
        use_median (bool): Aggregate samples by median instead of mean. Defaults to False.

    Returns:
        float | None: abundance value or None if not found.

    """
    # Assume first column is the index/search space
    index_col = omic_df.columns[0]
    # Filter rows where index_col matches analyte
    filtered = omic_df.filter(pl.col(index_col) == analyte)
    if filtered.is_empty():
        return None
    # Get numeric columns (excluding index)
    value_cols = [col for col in omic_df.columns if col != index_col]
    if not value_cols:
        return None
    # Compute mean across all value columns for the analyte row(s)
    all_vals = (
        filtered.select([pl.col(col).mean() for col in value_cols]).to_numpy().flatten()
    )
    if all_vals.size == 0:
        return None
    if use_median:  # use median?
        return float(np.median(all_vals))
    return float(all_vals.mean())


def validate_species(species: str = config.species) -> None:
    """Validate that the species is supported.

    Params:
        species (str): Species to verify. Defaults to species set in config.

    Raises:
        ValueError: Raised if species is not supported by XLRanker

    """
    supported_species = ["hsapiens", "mmusculus", "other"]
    if species not in supported_species:
        raise ValueError(
            f"Species {species} is not supported. Supported species are: {', '.join(supported_species)}"
        )
