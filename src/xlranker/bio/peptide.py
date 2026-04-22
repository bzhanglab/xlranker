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
        protein_linkages (dict[str, str]): mapping of protein name to linkage location

    """

    sequence: str
    mapped_proteins: list[str]
    linkage: int | None
    protein_linkages: dict[str, str]
    protein_locations: dict[str, tuple[int, int]]

    def __init__(
        self,
        sequence: str,
        linkage: int | None = None,
        mapped_proteins: list[str] | None = None,
        protein_linkages: dict[str, str] | None = None,
        protein_locations: dict[str, tuple[int, int]] | None = None,
    ):
        """Initialize the Peptide object with a sequence and mapped proteins.

        Args:
            sequence (str): amino acid sequence
            linkage (int | None): Optional 1-based index of the linkage
                on the amino acid sequence. Defaults to None.
            mapped_proteins (list[str], optional): list of protein names this sequence
                maps to. Defaults to [].
            protein_linkages (dict[str, str], optional): mapping of protein name to
                linkage location. Defaults to {}.
            protein_locations (dict[str, tuple[int, int]], optional): mapping of protein
                name to start and end index. Defaults to {}.

        """
        self.sequence = sequence
        self.linkage = None
        if linkage is not None and linkage <= len(self.sequence):
            self.linkage = linkage
        self.mapped_proteins = [] if mapped_proteins is None else list(mapped_proteins)
        self.protein_linkages = (
            {} if protein_linkages is None else dict(protein_linkages)
        )
        self.protein_locations = (
            {} if protein_locations is None else dict(protein_locations)
        )

    def __str__(self) -> str:
        """Get the string representation of this peptide.

        If linkage is present, it is added to string and separated by ':'.

        Returns:
            str: the amino acid sequence of this peptide.

        """
        if self.linkage is not None:
            return f"{self.sequence}:{self.linkage}"
        return self.sequence
