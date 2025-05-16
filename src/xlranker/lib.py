import sys
import logging
import polars as pl

from xlranker.util.mapping import PeptideMapper
from xlranker.util.readers import read_data_folder, read_network_file

from .bio import Protein, PeptidePair, ProteinPair

logger = logging.getLogger(__name__)


def setup_logging(
    verbose: bool = False, log_file: str | None = None, silent_all: bool = False
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

    network: list[PeptidePair]
    omic_data: list[pl.DataFrame]
    proteins: dict[str, Protein]

    def __init__(self, network: list[PeptidePair], omic_data: list[pl.DataFrame]):
        self.network = network
        self.omic_data = omic_data

    def build_protein_pairs(self, remove_intra: bool = True) -> None:
        """build protein pairs of the XLDataSet network

        Args:
            remove_intra (bool, optional): if true, only creates protein pairs between different proteins. Defaults to True.
        """
        self.proteins = {}
        for pair in self.network:
            protein_pairs: list[ProteinPair] = []
            for protein_a_name in pair.a.mapped_proteins:
                # Check if protein already encountered
                if protein_a_name not in self.proteins:
                    protein_a = Protein(protein_a_name)
                    self.proteins[protein_a_name] = protein_a
                else:
                    protein_a = self.proteins[protein_a_name]
                for protein_b_name in pair.b.mapped_proteins:
                    if protein_b_name not in self.proteins:
                        protein_b = Protein(protein_b_name)
                        self.proteins[protein_b_name] = protein_b
                    else:
                        protein_b = self.proteins[protein_b_name]
                    if (
                        remove_intra and protein_a != protein_b
                    ):  # Check if is intra linkage
                        new_pair = ProteinPair(protein_a, protein_b)
                        protein_pairs.append(new_pair)
                    else:  # Add pair regardless of protein identity
                        new_pair = ProteinPair(protein_a, protein_b)
                        protein_pairs.append(new_pair)
            pair.protein_pairs = protein_pairs

    @classmethod
    def load_from_network(
        cls,
        network_path: str,
        omics_data_folder: str,
        custom_mapping_path: str | None = None,
        is_fasta: bool = True,
        split_by: str | None = "|",
        split_index: int | None = 6,
    ) -> "XLDataSet":
        split_by = "|" if split_by is None else split_by
        split_index = 6 if split_index is None else split_index
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
