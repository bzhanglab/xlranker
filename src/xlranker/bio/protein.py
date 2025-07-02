from abc import ABC, abstractmethod


class ProteinNameExtractor(ABC):
    @abstractmethod
    def __init__(self) -> None:
        super().__init__()

    @abstractmethod
    def extract(self, isoform_name: str) -> str:
        pass


class NoExtractor(ProteinNameExtractor):
    def __init__(self) -> None:
        super().__init__()

    def extract(self, isoform_name: str) -> str:
        return isoform_name


class SplitExtractor(ProteinNameExtractor):
    split_by: str
    split_index: int

    def __init__(self, split_by: str, split_index: int) -> None:
        super().__init__()
        self.split_by = split_by
        self.split_index = split_index

    def extract(self, isoform_name: str) -> str:
        try:
            return isoform_name.split(self.split_by)[self.split_index]
        except IndexError:
            return isoform_name


class Protein:
    """Protein class that has the name and abundance for the protein

    Attributes:
        name (str): Name of the protein
        abundance (float | None): Abundance value of the protein

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
        self.name = name
        self.protein_name = protein_name
        self.abundances = abundances
        if main_column is None:
            self.main_column = next(iter(abundances))
        else:
            self.main_column = main_column

    def __eq__(self, value):
        return value.name == self.name

    def __hash__(self) -> int:
        return hash(self.name)

    def abundance(self) -> float | None:
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
    if a.abundance() is None:
        if b.abundance() is None:
            return (a, b)
        return (b, a)
    if b.abundance() is None:
        return (a, b)
    if b.abundance() <= a.abundance():
        return (a, b)
    return (b, a)
