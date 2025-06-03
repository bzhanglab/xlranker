"""Pipeline helper functions"""

from xlranker.lib import XLDataSet
from xlranker.ml.models import PrioritizationModel
from xlranker.parsimony.selection import ParsimonySelector


def run_full_pipeline(
    data_set: XLDataSet,
) -> XLDataSet:
    data_set.build_proteins()  # TODO: Determine if this should be done when loaded/initialized
    parsimony = ParsimonySelector(data_set)
    parsimony.run()
    model = PrioritizationModel(data_set)
    model.run_model()
    model.get_selections()
    return data_set
