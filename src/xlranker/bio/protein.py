from dataclasses import dataclass


@dataclass
class Protein:
    """Protein class that has the name and abundance for the protein

    Attributes:
        name (str): Name of the protein
        abundance (float | None): Abundance value of the protein
    """

    name: str
    abundance: float | None

    def __eq__(self, value: "Protein"):
        return value.name == self.name

    def __hash__(self):
        hash(self.name)


def sort_proteins(a: Protein, b: Protein) -> tuple[Protein, Protein]:
    """Takes into two Proteins and returns them so the first protein is higher abundant. Handles missing values.

    In the case of missing values for both proteins or equal values, the input order is maintained.

    Args:
        a (Protein): first protein
        b (Protein): second protein

    Returns:
        tuple[Protein, Protein]: protein tuple where the first protein is the higher abundant protein
    """
    if a.abundance is None:
        if b.abundance is None:
            return (a, b)
        return (b, a)
    if b.abundance is None:
        return (a, b)
    if b.abundance <= a.abundance:
        return (a, b)
    return (b, a)
