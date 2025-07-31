"""Protein objects and functions."""

from abc import ABC, abstractmethod


class ProteinNameExtractor(ABC):
    """Abstract class describing methods for a protein name extractor from isoform name."""

    @abstractmethod
    def __init__(self) -> None:
        """Initialize the protein name extractor."""
        super().__init__()

    @abstractmethod
    def extract(self, isoform_name: str) -> str:
        """Method to extract protein from isoform.

        Args:
            isoform_name (str): name of the isoform that needs extraction

        Returns:
            str: string of the protein name without any isoform information

        """
        pass


class NoExtractor(ProteinNameExtractor):
    """Protein name extractor that performs no operations.

    Used when no extraction needed.

    """

    def __init__(self) -> None:
        """Initialize a functionless name extractor."""
        super().__init__()

    def extract(self, isoform_name: str) -> str:
        """Returns the input isoform_name without any modifications.

        Args:
            isoform_name (str): name of the isoform

        Returns:
            str: same string as input

        """
        return isoform_name


class SplitExtractor(ProteinNameExtractor):
    """Protein name extractor that extracts a section from a isoform name using a splitting char and index.

    Args:
        split_by (str): string or char to split the isoform name by
        split_index (int): index to extract the protein name from after splitting the isoform name

    Attributes:
        split_by (str): string or char to split the isoform name by
        split_index (int): index to extract the protein name from after splitting the isoform name

    """

    split_by: str
    split_index: int

    def __init__(self, split_by: str, split_index: int) -> None:
        """Initialize the split extractor.

        Args:
            split_by (str): string or char to split the isoform name by
            split_index (int): index to extract the protein name from after splitting the isoform name

        """
        super().__init__()
        self.split_by = split_by
        self.split_index = split_index

    def extract(self, isoform_name: str) -> str:
        """Extract a section from isoform name using split_by and split_index.

        Args:
            isoform_name (str): isoform name needing extraction

        Returns:
            str: string of the protein name

        """
        try:
            return isoform_name.split(self.split_by)[self.split_index]
        except IndexError:
            return isoform_name


class Protein:
    """Protein class that has the name and abundance for the protein.

    Args:
        name (str): Name of the protein
        protein_name (str): name of the protein specific to the isoform level.
        abundances (dict[str, float | None]): Abundance values of the protein where keys is the name of the data and the value is the abundance
        main_column (str | None): column/key in abundance dictionary that represents the main abundance for the protein. Used for sorting the proteins. If none, uses first key in abundances dict.


    Attributes:
        name (str): Name of the protein
        abundances (dict[str, float | None]): Abundance values of the protein where keys is the name of the data and the value is the abundance
        main_column (str): column/key in abundance dictionary that represents the main abundance for the protein. Used for sorting the proteins.
        protein_name (str): name of the protein specific to the isoform level.

    """

    name: str
    abundances: dict[str, float | None]
    main_column: str
    protein_name: str

    def __init__(
        self,
        name: str,
        protein_name: str,
        abundances: dict[str, float | None] = {},
        main_column: str | None = None,
    ):
        """Protein class that has the name and abundance for the protein.

        Args:
            name (str): Name of the protein
            protein_name (str): name of the protein specific to the isoform level.
            abundances (dict[str, float | None]): Abundance values of the protein where keys is the name of the data and the value is the abundance
            main_column (str | None): column/key in abundance dictionary that represents the main abundance for the protein. Used for sorting the proteins. If none, uses first key in abundances dict.

        """
        self.name = name
        self.protein_name = protein_name
        self.abundances = abundances
        if main_column is None:
            self.main_column = next(iter(abundances))
        else:
            self.main_column = main_column

    def __eq__(self, value: object) -> bool:
        """Determine if object is equal to another.

        Args:
            value (object): object to compare to

        Returns:
            bool: If input is not a Protein object, returns False. Returns true if protein_name is the same for both proteins.

        """
        if isinstance(value, Protein):
            return value.protein_name == self.protein_name
        return False

    def __hash__(self) -> int:
        """Get hash representation of this object.

        Returns:
            int: hash using the protein name of this protein

        """
        return hash(self.protein_name)

    def abundance(self) -> float | None:
        """Get the representative abundance value for this protein.

        Uses the main_column attribute to get a value from the abundances dictionary.

        Returns:
            float | None: Abundance value if available. If main_column is invalid, returns None.

        """
        return self.abundances.get(self.main_column, None)


def sort_proteins(a: Protein, b: Protein) -> tuple[Protein, Protein]:
    """Takes into two Proteins and returns them so the first protein is higher abundant. Handles missing values.

    In the case of missing values for both proteins or equal values, the input order is maintained.

    Args:
        a (Protein): first protein
        b (Protein): second protein

    Returns:
        tuple[Protein, Protein]: protein tuple where the first protein is the higher abundant protein

    """
    a_abundance = a.abundance()
    b_abundance = b.abundance()
    if a_abundance is None:
        if b_abundance is None:
            return (a, b)
        return (b, a)
    if b_abundance is None:
        return (a, b)
    if b_abundance <= a_abundance:
        return (a, b)
    return (b, a)
