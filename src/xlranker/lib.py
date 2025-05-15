from dataclasses import dataclass
from typing import Self
import logging

from xlranker.util.mapping import PeptideMapper

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
        all_peptides = set()
        try:
            with open(network_path) as r:
                data = r.read().split("\n")
            for line in data:
                if "\t" in line:
                    vals = line.split("\t")
                    all_peptides.add(vals[0])
                    all_peptides.add(vals[1])
        except IndexError:
            logger.error(
                "Index out of bound. Make sure network is in the correct format."
            )
            raise IndexError()
        except FileNotFoundError:
            logger.error(f"File not found: {network_path}")
            raise FileNotFoundError
        all_peptides = list(all_peptides)
        mapper = PeptideMapper(
            mapping_table_path=custom_mapping_path,
            split_by=split_by,
            split_index=split_index,
        )
        mapping_results = mapper.map_sequences(all_peptides)
        pairs = []
