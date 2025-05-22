import sys
import logging
import polars as pl

from xlranker.util import get_abundance, get_pair_id
from xlranker.util.mapping import PeptideMapper
from xlranker.util.readers import read_data_folder, read_network_file

from .bio import Protein
from .bio.pairs import PeptidePair
from .bio.pairs import ProteinPair

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

    network: dict[str, PeptidePair]
    omic_data: dict[str, pl.DataFrame]
    proteins: dict[str, Protein]
    protein_pairs: dict[str, ProteinPair]

    def __init__(
        self, network: dict[str, PeptidePair], omic_data: dict[str, pl.DataFrame]
    ):
        self.network = network
        self.omic_data = omic_data
        self.protein_pairs = {}
        self.proteins = {}

    def build_proteins(self, remove_intra: bool = True) -> None:
        """build protein pairs of the XLDataSet network

        Args:
            remove_intra (bool, optional): if true, only creates protein pairs between different proteins. Defaults to True.
        """
        all_proteins: set[str] = set()
        for peptide_pairs in self.network.values():
            all_proteins = all_proteins.union(set(peptide_pairs.a.mapped_proteins))
            all_proteins = all_proteins.union(set(peptide_pairs.b.mapped_proteins))
        for protein in all_proteins:
            abundances = {}
            for omic_file in self.omic_data:
                abundances[omic_file] = get_abundance(
                    self.omic_data[omic_file], protein
                )
            self.proteins[protein] = Protein(protein, abundances)
        for peptide_pair in self.network.values():
            peptide_pair_id = get_pair_id(peptide_pair.a, peptide_pair.b)
            for protein_a_name in peptide_pairs.a.mapped_proteins:
                protein_a = self.proteins[protein_a_name]
                for protein_b_name in peptide_pair.b.mapped_proteins:
                    protein_b = self.proteins[protein_b_name]
                    if not remove_intra or (
                        remove_intra and protein_a != protein_b
                    ):  # Check if is intra linkage
                        protein_pair_id = get_pair_id(protein_a, protein_b)
                        if protein_pair_id not in self.protein_pairs:
                            new_pair = ProteinPair(protein_a, protein_b)
                            self.protein_pairs[protein_pair_id] = new_pair
                            peptide_pair.add_connection(protein_pair_id)
                            new_pair.add_connection(peptide_pair_id)
                        else:
                            self.protein_pairs[protein_pair_id].add_connection(
                                peptide_pair_id
                            )
                            peptide_pair.add_connection(protein_pair_id)

    @classmethod
    def load_from_network(
        cls,
        network_path: str,
        omics_data_folder: str,
        custom_mapping_path: str | None = None,
        is_fasta: bool = True,
        split_by: str | None = "|",
        split_index: int | None = 3,
    ) -> "XLDataSet":
        split_by = "|" if split_by is None else split_by
        split_index = 6 if split_index is None else split_index
        network = read_network_file(network_path)
        omic_data: dict[str, pl.DataFrame] = read_data_folder(omics_data_folder)
        peptide_sequences = set()
        for group in network.values():
            peptide_sequences.add(group.a.sequence)
            peptide_sequences.add(group.b.sequence)
        mapper = PeptideMapper(
            mapping_table_path=custom_mapping_path,
            split_by=split_by,
            split_index=split_index,
            is_fasta=is_fasta,
        )
        mapping_results = mapper.map_sequences(list(peptide_sequences))
        for group in network.values():
            group.a.mapped_proteins = mapping_results[group.a.sequence]
            group.b.mapped_proteins = mapping_results[group.b.sequence]
        return cls(network, omic_data)
