from dataclasses import dataclass
import json
from typing import Any


@dataclass
class Config:
    """config for XLRanker

    Attributes:
        fragile (bool): Default to False. If True, throw error on any warning
        detailed (bool): Default to False. If True, perform more analysis about dataset
        reduce_fasta (bool): Default to True. If True, when a gene has multiple sequences, only accept longest sequence as the canonical sequence.
    """

    fragile: bool = False  # Break on any warning
    detailed: bool = False  # Show more detailed information about dataset and analysis
    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    human_only: bool = True  # Is all data human only?


config = Config()


def set_config_from_dict(config_dict: dict[str, Any]) -> None:
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
