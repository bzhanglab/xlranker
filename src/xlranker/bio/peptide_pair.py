from xlranker.bio import Peptide, ProteinPair, GroupedEntity


class PeptidePair(GroupedEntity):
    """Peptide group that can contain multiple ProteinPairs and PeptidePairs"""

    a: Peptide
    b: Peptide
    protein_pairs: list[ProteinPair]

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide):
        super().__init__()
        self.a = peptide_a
        self.b = peptide_b
