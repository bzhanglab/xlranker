from dataclasses import dataclass

from xlranker.bio.protein import Protein


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[Protein]
