import polars as pl

from xlranker.bio.peptide import Peptide
from xlranker.bio.protein import Protein


def get_pair_id(a: Protein | Peptide, b: Protein | Peptide) -> str:
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
        return f"{name_a}+{name_b}"
    return f"{name_b}+{name_a}"


def get_abundance(omic_df: pl.DataFrame, analyte: str) -> float | None:
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
    mean_val = (
        filtered.select([pl.col(col).mean() for col in value_cols]).to_numpy().flatten()
    )
    if mean_val.size == 0:
        return None
    return float(mean_val.mean())
