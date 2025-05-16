from dataclasses import dataclass


@dataclass
class Peptide:
    sequence: str
    mapped_proteins: list[str]

    def __init__(self, sequence: str, mapped_proteins: list[str] = []):
        self.sequence = sequence
        self.mapped_proteins = mapped_proteins

    def __str__(self):
        return self.sequence
