from dataclasses import dataclass
import json


@dataclass
class Config:
    fragile: bool = False  # Break on any warning
    detailed: bool = False  # Show more detailed information about dataset and analysis
    reduce_fasta = False  # Reduce FASTA file by only keeping the largest sequence


config = Config()


def load_from_json(json_file: str) -> None:
    """set config to settings in JSON file

    Args:
        json_file (str): path to JSON file
    """
    with open(json_file) as r:
        json_obj = json.load(r)
    for key in json_obj:
        setattr(config, key, json_obj[key])
