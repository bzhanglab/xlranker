from xlranker.bio.peptide import Peptide


class PeptideGroup:
    """Peptide group that can contain multiple ProteinPairs and PeptidePairs"""

    a: Peptide
    b: Peptide

    def __init__(self, peptide_a: Peptide, peptide_b: Peptide):
        self.a = peptide_a
        self.b = peptide_b
