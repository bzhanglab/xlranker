"""Pipeline helper functions"""

from xlranker.lib import XLDataSet, get_final_network
from xlranker.ml.models import PrioritizationModel
from xlranker.parsimony.prioritize import ParsimonySelector, select_random
from xlranker.selection import ThresholdSelector


def run_full_pipeline(data_set: XLDataSet, threshold: float = 0.5) -> XLDataSet:
    """Run the full XLRanker pipeline.

    Args:
        data_set (XLDataSet): Cross-linking dataset that needs prioritization
        threshold (float): Score threshold for the expanded report

    Returns:
        XLDataSet: XLDataSet with full prioritization

    """
    data_set.build_proteins()  # TODO: Determine if this should be done when loaded/initialized
    parsimony = ParsimonySelector(data_set)
    parsimony.run()
    model = PrioritizationModel(data_set)
    model.run_model()
    get_final_network(data_set, ThresholdSelector(threshold))
    return data_set


def parsimony_only(data_set: XLDataSet, full_prioritization: bool = False) -> XLDataSet:
    """Run the XLRanker pipeline with only the parsimonious selection step.

    This will likely result in many PARSIMONY_AMBIGUOUS protein pairs. To avoid ambiguity, you can set full_prioritization to True. This will select one random pair as the representative pair for ambiguous groups.

    Args:
        data_set (XLDataSet): Cross-linking dataset that needs prioritization
        full_prioritization (bool): Default to False. If True, randomly select representative pairs for ambiguous groups.

    Returns:
        XLDataSet: XLDataSet with only parsimonious selection performed.

    """
    data_set.build_proteins()  # TODO: Determine if this should be done when loaded/initialized
    parsimony = ParsimonySelector(data_set)
    parsimony.run()
    if full_prioritization:
        select_random(data_set)
    return data_set
