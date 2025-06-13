from dataclasses import dataclass
from typing import override


@dataclass
class Peptide:
    """Peptide sequence object

    Attributes:
        sequence (str): Peptide sequence from peptide network
        mapped_proteins (list[str]): list of all proteins mapping to sequence

    """

    sequence: str
    mapped_proteins: list[str]

    def __init__(self, sequence: str, mapped_proteins: list[str] = []):
        self.sequence = sequence
        self.mapped_proteins = mapped_proteins

    @override
    def __str__(self) -> str:
        return self.sequence
