from importlib.resources import files
import polars as pl
import gzip
import pickle


def load_default_ppi() -> pl.DataFrame:
    """load default pre-generated table of known PPIs from parquet file into polars DataFrame.

    Returns:
        pl.DataFrame: Two column database with column names of P1 and P2 where P1 and P2 have a known PPI.
    """
    ppi_path = files("xlranker.data") / "ppi.parquet"
    return pl.read_parquet(str(ppi_path))


def load_gmts() -> list[list[set[str]]]:
    with gzip.open(str(files("xlranker.data") / "gmt.pkl.gz"), "rb") as r:
        return pickle.load(r)
