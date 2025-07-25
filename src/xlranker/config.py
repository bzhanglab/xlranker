import json
from dataclasses import dataclass, field
from typing import Any

DEFAULT_CONFIG = {
    "seed": None,
    "mapping_table": None,
    "is_fasta": True,
    "fasta_type": "UNIPROT",
    "only_human": True,
    "intra_in_training": False,
    "advanced": {
        "intra_in_training": False,  # allow intra in training data
    },
}


@dataclass
class AdvancedConfig:
    """Advanced config options for XLRanker

    Attributes:
        intra_in_training (bool): Default to False. If True, intra pairs are included in the positive set for model training. # TODO: May remove this option in future versions

    """

    intra_in_training: bool = False  # allow intra in training data


@dataclass
class MappingConfig:
    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    custom_table: str | None = None
    is_fasta: bool = True
    split_by: str | None = None
    split_index: int | None = None


@dataclass
class Config:
    """Config for XLRanker

    Attributes:
        fragile (bool): Default to False. If True, throw error on any warning
        detailed (bool): Default to False. If True, perform more analysis about dataset
        reduce_fasta (bool): Default to True. If True, when a gene has multiple sequences, only accept longest sequence as the canonical sequence.
        intra_in_training (bool): Default to False. If True, intra pairs are included in the positive set for model training.
        output (str): Default to "xlranker_output/". Directory where output files are saved.
        additional_null_values (list[str]): Default to []. Additional null values to consider when reading data files.

    """

    fragile: bool = False  # Break on any warning
    detailed: bool = False  # Show more detailed information about dataset and analysis
    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    human_only: bool = True  # Is all data human only?
    output: str = "xlranker_output/"  # output directory
    additional_null_values: list[str] = field(
        default_factory=list
    )  # additional null values to consider when reading data files
    advanced: AdvancedConfig = field(
        default_factory=AdvancedConfig
    )  # advanced config options
    mapping = field(default_factory=MappingConfig)


config = Config()


def set_config_from_dict(config_dict: dict[str, Any]) -> None:
    """set config from a dict object

    Args:
        config_dict (dict[str, Any]): dictionary with config settings

    """

    for key in config_dict:
        setattr(config, key, config_dict[key])


def load_from_json(json_file: str) -> None:
    """set config to settings in JSON file

    Args:
        json_file (str): path to JSON file

    """
    with open(json_file) as r:
        json_obj = json.load(r)
    for key in json_obj:
        setattr(config, key, json_obj[key])
