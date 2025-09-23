"""Creates a homolog database mapping gene symbols from different species to its human homolog."""

import polars as pl
import pickle
import gzip


def build_index(df: pl.LazyFrame) -> dict[int, str]:
    """Builds an index mapping DB Class Key to human gene symbol.

    Args:
        df (pl.LazyFrame): Dataframe containing homolog data

    Returns:
        dict[int, str]: index mapping DB Class Key to human gene symbol

    """
    df = df.filter(pl.col("Common Organism Name").eq(pl.lit("human"))).collect()
    return dict(zip(df["DB Class Key"], df["Symbol"]))


def create_homolog_db(
    df: pl.LazyFrame, species: str, index: dict[int, str], output_file: str
) -> None:
    """Outputs a pkl.gz file mapping gene symbols from a given species to its human homolog.

    Args:
        df (pl.LazyFrame): Dataframe containing homolog data
        species (str): Species to filter for (e.g. "mouse, laboratory")
        index (dict[int, str]): Index mapping DB Class Key to human gene symbol
        output_file (str): Output file path

    """
    df = (
        df.filter(pl.col("Common Organism Name").eq(pl.lit(species)))
        .with_columns(
            pl.col("DB Class Key")
            .replace_strict(index, return_dtype=pl.Utf8, default="NO_HOMOLOG")
            .alias("Homolog")
        )
        .collect()
    )
    with gzip.open(output_file, "wb") as f:
        pickle.dump(
            dict(
                zip(
                    df["Symbol"],
                    df["Homolog"],
                )
            ),
            f,
        )


if __name__ == "__main__":
    df = pl.scan_csv("data/HOM_ALLOrganism.rpt", separator="\t")
    if df is None:
        raise ValueError("Failed to load data")
    index = build_index(df)

    mk_dirs = [
        "output/mmusculus",
    ]

    for d in mk_dirs:
        import os

        os.makedirs(d, exist_ok=True)

    create_homolog_db(
        df, "mouse, laboratory", index, "output/mmusculus/homologs.pkl.gz"
    )
