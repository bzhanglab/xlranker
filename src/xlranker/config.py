"""Config related classes and methods. Contains the global config object."""

import json
from dataclasses import field
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict


class AdvancedConfig(BaseModel):
    """Advanced config options for XLRanker.

    Attributes:
        intra_in_training (bool): Default to False. If True, intra pairs are included in
            the positive set for model training.
            # TODO: May remove this option in future versions
        pair_separator (str): Default to "+". String separating pairs of peptides
            or proteins (i.e. ABC1+DEF2)
        save_model_files (bool): Default to False. Save model files required for SHAP
            value evaluation and further model evaluation
        binary_compartments (bool): Default to False. Use binary compartments
            (1 if any localization, 0 if none) rather than counts of compartments

    """

    model_config = ConfigDict(validate_assignment=True)
    intra_in_training: bool = False
    pair_separator: str = "+"
    save_model_files: bool = False
    binary_compartments: bool = False


class MappingConfig(BaseModel):
    """Mapping configuration object.

    Attributes:
        reduce_fasta (bool): If True, only keep longest sequence for duplicated protein
            entries
        custom_table (str | None): Path to custom table for peptide mapping.
            If None use default FASTA file.
        is_fasta (bool): True if custom table is in FASTA format
        split_by (str | None): string to split FASTA file for gene symbol extraction
        split_index (int | None): 0-based index of section containing gene symbol
            after string splitting
        fasta_type (str | None): UNIPROT or GENCODE fasta type. If None, will default to
            UNIPROT

    """

    model_config = ConfigDict(validate_assignment=True)
    reduce_fasta: bool = False  # Reduce FASTA file by only keeping the largest sequence
    custom_table: str | None = None
    is_fasta: bool = True
    split_by: str | None = None
    split_index: int | None = None
    fasta_type: str | None = "UNIPROT"


class Config(BaseModel):
    """Config for XLRanker.

    Attributes:
        network_path (str): Default to network.tsv. Path to the peptide network TSV file
        omic_data_folder (str): Default to omic_data/. Path to the folder containing the
            omics data in TSV format for the ML step
        fragile (bool): Default to False. If True, throw error on any warning
        detailed (bool): Default to False. If True, perform more analysis about dataset
        reduce_fasta (bool): Default to True.
            If True, when a gene has multiple sequences, only accept longest sequence
            as the canonical sequence.
        output (str): Default to "xlranker_output/".
            Directory where output files are saved.
        additional_null_values (list[str]): Default to []. Additional null values to
            consider when reading data files.
        advanced (AdvancedConfig): Advanced configuration options
        primary_column (str | None): Column name of which omic file should be the
            representative. If None, default to the first file alphabetically.
        mapping (MappingConfig): Configuration related to peptide sequence mapping.
        species (Literal["hsapiens", "mmusculus", "rnorvegicus", "other"]): Species name for mapping.
            Default to human.
        use_homologs (bool): If True, map to human homologs for non-human species.
        seed (int | None): Random state to initialize code. If None, seed is randomly
            generated. Defaults to None.
        threshold (float): Default to 0.5. Score threshold to accept pairs during
            ML selection.

    """

    network_path: str = "network.tsv"
    omic_data_folder: str = "omic_data/"
    omic_data_files: list[str] | None = None
    model_config = ConfigDict(validate_assignment=True)
    fragile: bool = False  # Break on any warning
    detailed: bool = False  # Show more detailed information about dataset and analysis
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
    species: Literal["hsapiens", "mmusculus", "rnorvegicus", "other"] = (
        "hsapiens"  # species name for mapping. Default to human.
    )
    use_homologs: bool = False  # If True, map to human homologs for non-human species
    seed: int | None = (
        None  # If None, set_seed will create random seed from 0 to 1000000
    )
    threshold: float = 0.5

    def is_safe(self) -> tuple[bool, str]:
        """Check if config is set in a valid way.

        Returns:
            (bool, str): True if all options set are compatible and error message.
        """
        is_safe = True
        error_msg = ""
        if self.network_path is None:  # input network must be provided
            is_safe = False
            error_msg += "Network path is none\n"
        if (
            not self.mapping.is_fasta and self.mapping.custom_table is None
        ):  # If not FASTA, custom table must be provided
            is_safe = False
            error_msg += "If not using FASTA, custom table must be provided\n"
        return (is_safe, error_msg)

    def check_is_safe(self) -> None:
        """Check if safe and raise error if error found."""
        (was_safe, error_msg) = self.is_safe()
        if not was_safe:
            message_string = (
                "Error Message" if error_msg.count("\n") == 1 else "Error Messages"
            )
            raise ValueError(
                f"Config not correctly set.\n{message_string}:\n{error_msg}"
            )


config = Config()


def update_model_in_place(model: BaseModel, updates: dict[str, Any]):
    """Recursively update a Pydantic model in-place from a dict."""
    for key, value in updates.items():
        if hasattr(model, key):
            current = getattr(model, key)
            if isinstance(current, BaseModel) and isinstance(value, dict):
                # recurse into nested BaseModel
                update_model_in_place(current, value)
            else:
                setattr(model, key, value)


def set_config_from_dict(config_dict: dict[str, Any]) -> None:
    """Set config from a dict object.

    Args:
        config_dict (dict[str, Any]): dictionary with config settings

    """
    update_model_in_place(config, config_dict)
    config.model_validate(config)


def load_from_json(json_file: str) -> None:
    """Set config to settings in JSON file.

    Args:
        json_file (str): path to JSON file

    """
    with open(json_file) as r:
        json_obj = json.load(r)
    set_config_from_dict(json_obj)


def config_to_dict(config_obj: Config) -> dict[str, Any]:
    """Convert Config object to a dictionary.

    Args:
        config_obj (config): config to convert to dictionary like object.

    Returns:
        dict[str, Any] | list[dict[str, Any]]: JSON/YAML serializable object
            representing the input config

    """
    return config_obj.model_dump()


DEFAULT_CONFIG = config_to_dict(Config())
