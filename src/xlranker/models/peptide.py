from dataclasses import dataclass

from xlranker.models.protein import Protein


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[Protein]
