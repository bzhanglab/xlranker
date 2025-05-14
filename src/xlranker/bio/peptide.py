from dataclasses import dataclass

from xlranker.bio.protein import Protein


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[Protein]

    def __init__(self, sequence: str, mapped_proteins: list[Protein] = []):
        self.sequence = sequence
        self.mapped_proteins = mapped_proteins
