# Input file is a TSV file with two columns constituting a known PPI

import polars as pl

input_file = "ppi.tsv"

df = pl.read_csv(
    input_file,
    separator="\t",
    has_header=False,
    new_columns=["P1", "P2"],
)

# consistent case

df = df.with_columns(
    [
        pl.col("P1").str.to_uppercase().alias("P1"),
        pl.col("P2").str.to_uppercase().alias("P2"),
    ]
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
df.write_parquet("ppi.parquet")
