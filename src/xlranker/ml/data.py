from importlib.resources import files
import gzip
import polars as pl


def load_default_ppi() -> pl.DataFrame:
    ppi_path = files("xlranker.data") / "ppi.parquet"
    return pl.read_parquet(ppi_path)
