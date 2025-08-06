"""Config related classes and methods. Contains the global config object."""

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
    "primary_column": None,
    "advanced": {
        "intra_in_training": False,  # allow intra in training data
    },
}


@dataclass
class AdvancedConfig:
    """Advanced config options for XLRanker.

    Attributes:
        intra_in_training (bool): Default to False. If True, intra pairs are included in the positive set for model training. # TODO: May remove this option in future versions

    """

    intra_in_training: bool = False  # allow intra in training data


@dataclass
class MappingConfig:
    """Mapping configuration object.

    Attributes:
        reduce_fasta (bool): If True, only keep longest sequence for duplicated protein entries
        custom_table (str | None): Path to custom table for peptide mapping. If None use default FASTA file.
        is_fasta (bool): True if custom table is in FASTA format
        split_by (str | None): string to split FASTA file for gene symbol extraction
        split_index (int | None): 0-based index of section containing gene symbol after string splitting
        fasta_type (str | None): UNIPROT or GENCODE fasta type. If None, will default to UNIPROT

    """

    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    custom_table: str | None = None
    is_fasta: bool = True
    split_by: str | None = None
    split_index: int | None = None
    fasta_type: str | None = "UNIPROT"


@dataclass
class Config:
    """Config for XLRanker.

    Attributes:
        fragile (bool): Default to False. If True, throw error on any warning
        detailed (bool): Default to False. If True, perform more analysis about dataset
        reduce_fasta (bool): Default to True. If True, when a gene has multiple sequences, only accept longest sequence as the canonical sequence.
        human_only (bool): True if all data is human only.
        output (str): Default to "xlranker_output/". Directory where output files are saved.
        additional_null_values (list[str]): Default to []. Additional null values to consider when reading data files.
        advanced (AdvancedConfig): Advanced configuration options
        primary_column (str | None): Column name of which omic file should be the representative. If None, default to the first file alphabetically.
        mapping (MappingConfig): Configuration related to peptide sequence mapping.

    """

    fragile: bool = False  # Break on any warning
    detailed: bool = False  # Show more detailed information about dataset and analysis
    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    human_only: bool = True  # Is all data human only?
    output: str = "xlranker_output/"  # output directory
    primary_column: str | None = (
        None  # Which omic file should be the representative omic set for ordering
    )
    additional_null_values: list[str] = field(
        default_factory=list
    )  # additional null values to consider when reading data files
    advanced: AdvancedConfig = field(
        default_factory=AdvancedConfig
    )  # advanced config options
    mapping: MappingConfig = field(default_factory=MappingConfig)


config = Config()


def set_config_from_dict(config_dict: dict[str, Any]) -> None:
    """Set config from a dict object.

    Args:
        config_dict (dict[str, Any]): dictionary with config settings

    """
    for key in config_dict:
        if key == "advanced":
            for sub_key in config_dict[key]:
                setattr(config.advanced, sub_key, config_dict[key][sub_key])
        elif key == "mapping":
            for sub_key in config_dict[key]:
                setattr(config.mapping, sub_key, config_dict[key][sub_key])
        else:
            setattr(config, key, config_dict[key])


def load_from_json(json_file: str) -> None:
    """Set config to settings in JSON file.

    Args:
        json_file (str): path to JSON file

    """
    with open(json_file) as r:
        json_obj = json.load(r)
    set_config_from_dict(json_obj)


def config_to_dict(config_obj: Config) -> dict[str, Any] | list[dict[str, Any]]:
    """Convert Config object to a dictionary.

    Args:
        config_obj (config): config to convert to dictionary like object.

    Returns:
        dict[str, Any] | list[dict[str, Any]]: JSON/YAML serializable object representing the input config

    """

    def dataclass_to_dict(obj) -> dict[str, Any] | list[dict[str, Any]]:
        if hasattr(obj, "__dataclass_fields__"):
            result = {}
            for field_name in obj.__dataclass_fields__:
                value = getattr(obj, field_name)
                result[field_name] = dataclass_to_dict(value)
            return result
        elif isinstance(obj, list):
            return [dataclass_to_dict(item) for item in obj]
        else:
            return obj

    return dataclass_to_dict(config_obj)
