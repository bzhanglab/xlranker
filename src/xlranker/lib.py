from dataclasses import dataclass
from typing import Self
import logging

from xlranker.util.mapping import PeptideMapper
from xlranker.util.readers import read_network_file

from .bio import Protein, Peptide, ProteinPair, PeptideGroup

logger = logging.getLogger(__name__)


class XLDataSet:
    """Crosslinking dataset object.

    Attributes:
        peptides (dict[str, Peptide]): dictionary of Peptide objects
    """

    peptides: dict[str, Peptide]

    def __init__(self):
        pass

    def get_all_proteins(self) -> list[Protein]:
        all_proteins = []
        for peptide in self.peptides.values():
            all_proteins.extend(peptide.mapped_proteins)
        return all_proteins

    def load_from_network(
        network_path: str,
        custom_mapping_path: str | None = None,
        is_fasta: bool = True,
        split_by: str = "|",
        split_index: int = 6,
    ) -> Self:
        network = read_network_file(network_path)
        peptide_sequences = set()
        for group in network:
            peptide_sequences.add(group.a.sequence)
            peptide_sequences.add(group.b.sequence)
        mapper = PeptideMapper(
            mapping_table_path=custom_mapping_path,
            split_by=split_by,
            split_index=split_index,
        )
        mapping_results = mapper.map_sequences(list(peptide_sequences))
        pairs = []
