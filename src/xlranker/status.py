from enum import Enum, auto


class ReportStatus(Enum):
    """Pair reporting status"""

    CONSERVATIVE = (
        auto(),
        "High confidence pairs, parsimonious unambiguous",
    )
    MINIMAL = (
        auto(),
        "Medium confidence pairs, parsimonious ambiguous included, all peptides represented",
    )
    EXPANDED = (
        auto(),
        "All pairs, including non-parsimonious pairs, high scoring ML pairs",
    )
    ALL = auto(), "All pairs, regardless of status"
    NONE = auto(), "No assigned status"


class PrioritizationStatus(Enum):
    """Prioritization status for a protein pair"""

    NOT_ANALYZED = auto()
    "No analysis performed yet"

    # Parsimony-based statuses

    PARSIMONY_NOT_SELECTED = auto()
    "Another entity was selected in group or cannot be selected"
    PARSIMONY_PRIMARY_SELECTED = auto()
    "Selected as the primary representative for group."
    PARSIMONY_SECONDARY_SELECTED = auto()
    "Selected as a secondary representative for group. Only possible for intra pairs."
    PARSIMONY_AMBIGUOUS = auto()
    "No clear candidate from parsimony analysis. Needs ML model."

    # Machine Learning-based statuses

    ML_NOT_SELECTED = auto()
    "Other candidate had higher score in group"
    ML_PRIMARY_SELECTED = auto()
    "Highest ML score in group or primary selection"
    ML_SECONDARY_SELECTED = auto()
    "High confidence ML score in group or secondary selection"
