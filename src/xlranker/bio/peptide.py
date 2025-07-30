from dataclasses import dataclass


@dataclass
class Peptide:
    """Peptide sequence object.

    Args:
        sequence (str): peptide amino acid sequence
        mapped_proteins (list[str]): list of proteins the peptide sequence maps to

    Attributes:
        sequence (str): Peptide sequence from peptide network
        mapped_proteins (list[str]): list of all proteins mapping to sequence

    """

    sequence: str
    mapped_proteins: list[str]

    def __init__(self, sequence: str, mapped_proteins: list[str] = []):
        self.sequence = sequence
        self.mapped_proteins = mapped_proteins

    def __str__(self) -> str:
        return self.sequence
