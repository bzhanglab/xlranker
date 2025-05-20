from .protein import Protein
from .peptide import Peptide

from enum import Enum, auto

__all__ = ["Protein", "Peptide"]


class PrioritizationStatus(Enum):
    NOT_ANALYZED = auto()  # No analysis performed yet

    # Parsimony-based statuses
    PARSIMONY_NOT_SELECTED = (
        auto()
    )  # Another entity was selected, or cannot be selected
    PARSIMONY_SELECTED = auto()  # Selected as the representative
    PARSIMONY_AMBIGUOUS = auto()  # No clear candidate

    # Machine Learning-based statuses
    ML_NOT_SELECTED = auto()  # Lower ML score
    ML_SELECTED = auto()  # Highest ML score
