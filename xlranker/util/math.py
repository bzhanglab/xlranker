from xlranker.xlranker import Protein


def set_protein_order(a: Protein, b: Protein) -> tuple[Protein, Protein]:
    if a.abundance is None:
        if b.abundance is None:
            return (a, b)
        return (b, a)
    if b.abundance is None:
        return (a, b)
    if b.abundance < a.abundance:
        return (a, b)
    return (b, a)
