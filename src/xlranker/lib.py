from dataclasses import dataclass
from typing import Self

from .bio import Protein, Peptide, ProteinPair, PeptideGroup


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
        all_sequences = None
