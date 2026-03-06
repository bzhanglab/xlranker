"""Creates a homolog database mapping gene symbols from different species to its human homolog."""  # noqa: E501

import gzip
import pickle

import polars as pl


def build_index(df: pl.LazyFrame) -> dict[int, str]:
    """Builds an index mapping DB Class Key to human gene symbol.

    Args:
        df (pl.LazyFrame): Dataframe containing homolog data

    Returns:
        dict[int, str]: index mapping DB Class Key to human gene symbol

    """
    full_df = df.collect()
    return dict(zip(full_df["RAT_GENE_SYMBOL"], full_df["HUMAN_ORTHOLOG_SYMBOL"]))


def create_homolog_db(
    df: pl.LazyFrame, index: dict[int, str], output_file: str
) -> None:
    """Output mapping table of gene symbols from a species to its human homolog.

    File format is pkl.gz.

    Args:
        df (pl.LazyFrame): Dataframe containing homolog data
        index (dict[int, str]): Index mapping DB Class Key to human gene symbol
        output_file (str): Output file path

    """
    full_df = (
        df.with_columns(
            pl.col("RAT_GENE_SYMBOL")
            .replace_strict(index, return_dtype=pl.Utf8, default="NO_HOMOLOG")
            .alias("Homolog"),
            pl.col("RAT_GENE_SYMBOL").alias("Symbol"),
        )
        .unique(subset=["Symbol", "Homolog"])
        .collect()
    )
    with gzip.open(output_file, "wb") as f:
        pickle.dump(
            dict(
                zip(
                    full_df["HUMAN_ORTHOLOG_SYMBOL"],
                    full_df["Homolog"],
                )
            ),
            f,
        )


if __name__ == "__main__":
    df = pl.scan_csv("data/RGD_ORTHOLOGS_NCBI.txt", separator="\t", comment_prefix="#")
    if df is None:
        raise ValueError("Failed to load data")
    index = build_index(df)

    mk_dirs = [
        "output/rnorvegicus",
    ]

    for d in mk_dirs:
        import os

        os.makedirs(d, exist_ok=True)

    create_homolog_db(df, index, "output/rnorvegicus/homologs.pkl.gz")
