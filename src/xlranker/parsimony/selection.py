from xlranker.bio import PeptidePair
import logging

logger = logging.getLogger(__name__)


class ParsimonySelector:
    def __init__(
        self, peptide_groups: list[PeptidePair], mapping_table: dict[str, list[str]]
    ):
        pass
