"""Creates a Polars dataframe with P1 and P2 columns where each row is a set of proteins with known PPI.

Input file is a TSV file with two columns constituting a known PPI
"""

import polars as pl


def process_ppi(input_file: str, output_file: str = "ppi.parquet"):
    """Process a TSV file containing known PPIs into a Polars DataFrame.

    Args:
        input_file (str): Path to the TSV file containing known PPIs.
        output_file (str): Path to the output Parquet file to save the processed DataFrame.

    Returns:
        pl.DataFrame: A Polars DataFrame with two columns, P1 and P2, where each row represents a known PPI.
    """

    df = pl.read_csv(
        input_file,
        separator="\t",
        has_header=False,
        new_columns=["P1", "P2"],
    )

    # make sure first column is lower value
    df = df.with_columns(
        [
            pl.when(pl.col("P1") > pl.col("P2"))
            .then(pl.col("P2"))
            .otherwise(pl.col("P1"))
            .alias("P1"),
            pl.when(pl.col("P1") > pl.col("P2"))
            .then(pl.col("P1"))
            .otherwise(pl.col("P2"))
            .alias("P2"),
        ]
    )

    print(df.head())
    df.write_parquet(output_file)


dirs_to_create = ["output", "output/hsapiens", "output/mmusculus"]
for d in dirs_to_create:
    import os

    if not os.path.exists(d):
        os.mkdir(d)

process_ppi("data/hsapiens/ppi.tsv", "output/hsapiens/ppi.parquet")
process_ppi("data/mmusculus/ppi.tsv", "output/mmusculus/ppi.parquet")
