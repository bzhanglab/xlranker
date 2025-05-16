from xlranker.bio import PeptidePair, ProteinPair
from xlranker.lib import XLDataSet
import logging

logger = logging.getLogger(__name__)


class ParsimonyGroup:
    protein_pairs: list[ProteinPair]
    peptide_pairs: list[PeptidePair]


class ParsimonySelector:
    data_set: XLDataSet

    def __init__(self, data_set: XLDataSet):
        self.data_set = data_set

    def prioritize_pair(self, peptide_pair: PeptidePair) -> None:
        peptide_pair.protein_pairs

    def prioritize(self) -> None:
        for pair in self.data_set.network:
            pass
