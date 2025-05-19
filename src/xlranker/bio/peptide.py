from dataclasses import dataclass
from typing import override


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[str]

    def __init__(self, sequence: str, mapped_proteins: list[str] = []):
        self.sequence = sequence
        self.mapped_proteins = mapped_proteins

    @override
    def __str__(self) -> str:
        return self.sequence
