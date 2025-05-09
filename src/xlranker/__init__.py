from dataclasses import dataclass

from .models import Protein, Peptide, ProteinPair, PeptideGroup


class XLDataSet:
    peptides: dict[str, Peptide]

    def __init__(self):
        pass

    def get_all_proteins(self) -> list[Protein]:
        all_proteins = []
        for peptide in self.peptides.values():
            all_proteins.extend(peptide.mapped_proteins)
        return all_proteins
