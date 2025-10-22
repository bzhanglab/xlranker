"""Retrieve datasets required for the pipeline."""

# TODO: Need to add mouse data files and ignore species if not exist.

import gzip
import lzma
import pickle
import tarfile
import tempfile
from importlib.resources import files

import polars as pl
from xlranker.config import config
from xlranker.util import validate_species


def load_localization_data(
    species: str = config.species,
) -> dict[str, dict[str, set[str]]]:
    """Load protein localization data from COMPARTMENTS and Human Protein Atlas.

    Returns:
         dict[str, dict[str, set[str]]]: Key is the data set (compartments or HumanProteinAtlas) and the values are dicts containing the annotated locations for a gene.

    """
    validate_species()
    try:
        coloc_path = files(f"xlranker.data.species.{config.species}") / "coloc.pkl.gz"
        with gzip.open(str(coloc_path), "rb") as r:
            return pickle.load(r)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"No localization data found for species '{species}'. Please report issue or provide your own localization data."
        ) from e


def load_default_ppi(species: str = config.species) -> pl.DataFrame:
    """Load default pre-generated table of known PPIs from parquet file into polars DataFrame. Human only.

    Returns:
        pl.DataFrame: Two column database with column names of P1 and P2 where P1 and P2 have a known PPI.

    """
    validate_species()

    try:
        ppi_path = files(f"xlranker.data.species.{species}") / "ppi.parquet"
        return pl.read_parquet(str(ppi_path))
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"No PPI database found for species '{species}'. Please report issue or provide your own PPI database."
        ) from e


def load_gmts(species=config.species) -> list[list[set[str]]]:
    """Load gmt collection. Used to determine negative sets for ML step.

    Contains Gene Ontology Biological Process and Reactome.

    Returns:
        list[list[set[str]]]: list of gmts, which are collections of sets.

    """
    validate_species()

    with gzip.open(
        str(files(f"xlranker.data.species.{species}") / "gmt.pkl.gz"), "rb"
    ) as r:
        return pickle.load(r)


def load_homologs() -> dict[str, str] | None:
    """Load homolog database for species.

    Returns:
        dict[str, str]: key is species gene symbol, value is homolog Human gene symbol
    """
    validate_species()

    if config.species == "hsapiens":
        return None
    homolog_path = files(f"xlranker.data.species.{config.species}") / "homologs.pkl.gz"
    with gzip.open(str(homolog_path), "rb") as r:
        return pickle.load(r)


def get_default_fasta() -> str:
    """Extract and write a default UNIPROT FASTA database from May 2022.

    Raises:
        FileNotFoundError: Raised if FASTA file has issues with extraction.

    Returns:
        str: path to temporary file containing FASTA file.

    """
    validate_species()

    fasta_path = str(
        files(f"xlranker.data.species.{config.species}") / "uniprot_5_22.fa.tar.xz"
    )
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
