"""Pipeline helper functions"""

from xlranker.lib import XLDataSet
from xlranker.ml.models import PrioritizationModel
from xlranker.parsimony.selection import ParsimonySelector


def run_full_pipeline(
    data_set: XLDataSet,
) -> XLDataSet:
    """Run the full XLRanker pipeline

    Args:
        data_set (XLDataSet): Cross-linking dataset that needs prioritization

    Returns:
        XLDataSet: XLDataSet with full prioritization

    """
    data_set.build_proteins()  # TODO: Determine if this should be done when loaded/initialized
    parsimony = ParsimonySelector(data_set)
    parsimony.run()
    model = PrioritizationModel(data_set)
    model.run_model()
    model.get_selections()
    return data_set


def parsimony_only(data_set: XLDataSet) -> XLDataSet:
    """Run the XLRanker pipeline with only the parsimonious selection step.

    This will likely result in many PARSIMONY_AMBIGUOUS protein pairs. You can have them removed, or run the xlranker.parsimony.selection.random(data_set) function to randomly select pairs from each unresolved group.

    Args:
        data_set (XLDataSet): Cross-linking dataset that needs prioritization

    Returns:
        XLDataSet: XLDataSet with only parsimonious selection performed.

    """
    data_set.build_proteins()  # TODO: Determine if this should be done when loaded/initialized
    parsimony = ParsimonySelector(data_set)
    parsimony.run()
    return data_set
