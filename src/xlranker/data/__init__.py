"""Retrieve datasets required for the pipeline."""

import gzip
import lzma
import pickle
import tarfile
import tempfile
from importlib.resources import files

import polars as pl


def load_default_ppi() -> pl.DataFrame:
    """Load default pre-generated table of known PPIs from parquet file into polars DataFrame. Human only.

    Returns:
        pl.DataFrame: Two column database with column names of P1 and P2 where P1 and P2 have a known PPI.

    """
    ppi_path = files("xlranker.data") / "ppi.parquet"
    return pl.read_parquet(str(ppi_path))


def load_gmts() -> list[list[set[str]]]:
    """Load gmt collection. Used to determine negative sets for ML step.

    Contains Gene Ontology Biological Process and Reactome.

    Returns:
        list[list[set[str]]]: list of gmts, which are collections of sets.

    """
    with gzip.open(str(files("xlranker.data") / "gmt.pkl.gz"), "rb") as r:
        return pickle.load(r)


def get_default_fasta() -> str:
    """Extract and write a default UNIPROT FASTA database from May 2022.

    Raises:
        FileNotFoundError: Raised if FASTA file has issues with extraction.

    Returns:
        str: path to temporary file containing FASTA file.

    """
    fasta_path = str(files("xlranker.data") / "uniprot_5_22.fa.tar.xz")
    with lzma.open(fasta_path) as r:
        with tarfile.open(fileobj=r) as tar:
            fa_file = next(
                (m for m in tar.getmembers() if m.name.endswith(".fa")), None
            )
            if not fa_file:
                raise FileNotFoundError(
                    "No .fa file found in the tar archive. Please report issue."
                )
            temp_dir = tempfile.mkdtemp()
            tar.extract(fa_file, path=temp_dir)
            temp_fa_path = f"{temp_dir}/{fa_file.name}"
            return temp_fa_path
