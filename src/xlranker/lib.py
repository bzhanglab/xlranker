import sys
from typing import Self
import logging
import polars as pl

from xlranker.util.mapping import PeptideMapper
from xlranker.util.readers import read_data_folder, read_network_file

from .bio import Protein, PeptideGroup

logger = logging.getLogger(__name__)


def setup_logging(
    verbose: bool = False, log_file: str = None, silent_all: bool = False
):
    if silent_all:
        # Remove all handlers and disable logging
        logging.getLogger().handlers.clear()
        logging.disable(logging.CRITICAL + 1)
        return
    level = logging.DEBUG if verbose else logging.INFO

    # Create root logger
    logger = logging.getLogger()
    logger.setLevel(level)

    # Console handler (stderr)
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(level)
    console_formatter = logging.Formatter("[%(levelname)s] %(message)s")
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    # Optional file handler
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_formatter = logging.Formatter(
            "%(asctime)s [%(levelname)s] %(name)s: %(message)s"
        )
        file_handler.setFormatter(file_formatter)
        logger.addHandler(file_handler)


class XLDataSet:
    """XLRanker cross-linking dataset object.

    Attributes:
        peptides (dict[str, Peptide]): dictionary of Peptide objects
    """

    network: list[PeptideGroup]
    omic_data: list[pl.DataFrame]

    def __init__(self, network: list[PeptideGroup], omic_data: list[pl.DataFrame]):
        self.network = network
        self.omic_data = omic_data

    def get_all_proteins(self) -> list[Protein]:
        all_proteins = []
        for peptide in self.peptides.values():
            all_proteins.extend(peptide.mapped_proteins)
        return all_proteins

    @classmethod
    def load_from_network(
        cls: Self,
        network_path: str,
        omics_data_folder: str,
        custom_mapping_path: str | None = None,
        is_fasta: bool = True,
        split_by: str = "|",
        split_index: int = 6,
    ) -> Self:
        network = read_network_file(network_path)
        omic_data: list[pl.DataFrame] = read_data_folder(omics_data_folder)
        peptide_sequences = set()
        for group in network:
            peptide_sequences.add(group.a.sequence)
            peptide_sequences.add(group.b.sequence)
        mapper = PeptideMapper(
            mapping_table_path=custom_mapping_path,
            split_by=split_by,
            split_index=split_index,
            is_fasta=is_fasta,
        )
        mapping_results = mapper.map_sequences(list(peptide_sequences))
        for group in network:
            group.a.mapped_proteins = mapping_results[group.a.sequence]
            group.b.mapped_proteins = mapping_results[group.b.sequence]
        return cls(network, omic_data)
