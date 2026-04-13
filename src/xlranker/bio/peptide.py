"""Peptide sequence objects and functions."""

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
    linkage: int | None

    def __init__(
        self, sequence: str, linkage: int | None = None, mapped_proteins: list[str] = []
    ):
        """Initialize the Peptide object with a sequence and mapped proteins.

        Args:
            sequence (str): amino acid sequence
            linkage (int | None): Optional 1-based index of the linkage
                on the amino acid sequence. Defaults to None.
            mapped_proteins (list[str], optional): list of protein names this sequence
                maps to. Defaults to [].

        """
        self.sequence = sequence
        if linkage is not None and linkage <= len(self.sequence):
            self.linkage = linkage
        self.mapped_proteins = mapped_proteins

    def __str__(self) -> str:
        """Get the string representation of this peptide.

        If linkage is present, it is added to string and separated by ':'.

        Returns:
            str: the amino acid sequence of this peptide.

        """
        if self.linkage is not None:
            return f"{self.sequence}:{self.linkage}"
        return self.sequence
