"""
Model Process:

1. Identify Positive Dataset
    - All representative pairs from parsimonious selection
2. Generate Negative Dataset
    - Random protein pairs that are not selected pairs
"""

from xlranker.lib import XLDataSet


class PrioritizationModel:
    dataset: XLDataSet

    def __init__(self, dataset: XLDataSet):
        pass

    def get_positives(self):
        for peptide in self.dataset.peptides.values():
            pass
