"""Command Line Interface for XLRanker."""

import json
import logging
import os
from pathlib import Path
from typing import Annotated, Any

import cyclopts
import questionary
import yaml
from questionary import Choice

from xlranker.config import (
    DEFAULT_CONFIG,
    Config,
    config_to_dict,
    set_config_from_dict,
)
from xlranker.config import config as xlr_config
from xlranker.lib import XLDataSet, setup_logging
from xlranker.pipeline import run_full_pipeline
from xlranker.util import set_seed
from xlranker.util.readers import base_name

app = cyclopts.App()
logger = logging.getLogger(__name__)


def load_config(path: str) -> None:
    """Load a JSON or YAML config to dictionary from path.

    Args:
        path (str): path to JSON or YAML file

    Raises:
        ValueError: raised if config does not end in .json, .yaml, or .yml

    """
    input_dict = {}
    if path.lower().endswith(".json"):
        input_dict = json.load(open(path))
    elif path.lower().endswith(".yaml") or path.lower().endswith(".yml"):
        input_dict = yaml.safe_load(open(path))
    else:
        raise ValueError("Unsupported config file format.")
    set_config_from_dict(input_dict)


def save_config(path: str, config_obj: dict[str, Any]) -> None:
    """Save config to path in either JSON or YAML format.

    Args:
        path (str): path to write config to
        config_obj (dict[str, Any]): config to save

    Raises:
        ValueError: raised if path does not end in .json, .yaml, or .yml

    """
    path = path.lower()
    if path.endswith(".json"):
        return json.dump(config_obj, open(path, "w"))
    elif path.endswith(".yaml") or path.endswith(".yml"):
        return yaml.dump(config_obj, open(path, "w"))
    raise ValueError("Unsupported config file format.")


def is_folder(path_to_validate: str) -> bool | str:
    """Check if the input path is a folder for form verification.

    Args:
        path_to_validate (str): path to check

    Returns:
        bool | str: returns True if is_folder, else provide error message.

    """
    return (
        True
        if not os.path.isfile(path_to_validate)
        else "Input a directory, not an existing file."
    )


@app.command()
def init(
    default: bool = False,
    output: Annotated[str | None, cyclopts.Parameter(name=["--output", "-o"])] = None,
) -> None:  # noqa: DOC105
    """Initialize a config file, with an optional interactive form.

    If no default flag provided, config is created through a interactive form.

    Args:
        default (bool, optional): Create a simple default config. Defaults to False.
        output (Annotated[str | None, cyclopts.Parameter], optional): Output config file
            Can either be JSON or YAML format. Defaults to None.

    Raises:
        ValueError: raises ValueError if output is not set when default is passed

    """
    if default:
        if output is None:
            output = "xlranker_config.yaml"
        save_config(output, DEFAULT_CONFIG)
        print(f"Saved config to {output}.")
        return

    network = questionary.path(
        "Path to your peptide sequence network:",
    ).ask()
    omic_data = questionary.path(
        "Path to omic data folder:",
        only_directories=True,
        validate=is_folder,
    ).ask()
    globs = [str(p) for p in list(Path(omic_data).glob("*"))]
    primary_column = None
    if len(globs) > 0:
        selected_file = questionary.select(
            "Which file is the primary abundance value?", choices=globs
        ).ask()
        primary_column = base_name(str(selected_file))

    mapping_table = questionary.select(
        "What mapping table will you use?",
        choices=[
            "Custom FASTA database",
            "TSV Table",
            "Default: Human UNIPROT from May 2025",
        ],
    ).ask()
    fasta_type = None
    match mapping_table:
        case "Custom FASTA database":
            is_fasta = True
            fasta_type = questionary.select(
                "Type of FASTA file:", choices=["GENCODE", "UNIPROT"]
            ).ask()
            mapping_table_path = questionary.path(
                "Path to fasta file:",
                validate=lambda x: True
                if x.lower().endswith(".fa") or x.lower().endswith(".fasta")
                else "Please input a FASTA file (.fasta or .fa)",
            ).ask()
            if fasta_type == "GENCODE":
                print(
                    "\nGENCODE additional configuration to read the FASTA file.\nSee https://bzhanglab.github.io/xlranker/latest/usage/input_data/fasta/ for more information.\n"  # noqa: E501
                )
        case "TSV Table":
            is_fasta = False
            mapping_table_path = questionary.path("Path to TSV file:").ask()
        case _:
            is_fasta = True
            mapping_table_path = None

    species = questionary.select(
        "What species is your data from?",
        choices=[
            Choice("Human (hsapiens)", "hsapiens"),
            Choice("Mouse (mmusculus)", "mmusculus"),
        ],
    ).ask()
    if species != "hsapiens":
        use_homologs = questionary.confirm(
            "Do you want to use human homologs? (Default is Yes)", True
        ).ask()
    else:
        use_homologs = False
    new_config = Config()

    # Set variables

    new_config.network_path = network
    new_config.omic_data_folder = omic_data
    new_config.species = species
    new_config.use_homologs = use_homologs
    new_config.primary_column = primary_column

    new_config.mapping.fasta_type = fasta_type
    new_config.mapping.custom_table = mapping_table_path
    new_config.mapping.is_fasta = is_fasta

    good_output = output is not None and any(
        output.endswith(s) for s in [".json", ".yml", ".yaml"]
    )
    if not good_output and output is not None:
        print("Output path did not end in .yaml, .yml, or .json. Asking for output")
    if good_output:
        output_path = output  # type: ignore
    else:
        output_path: str = questionary.path(
            "Output file for config (JSON or YAML format):",
            validate=lambda x: True
            if x.lower().endswith(".json")
            or x.lower().endswith(".yaml")
            or x.lower().endswith(".yml")
            else "File must end with .json, .yaml, or .yml",
        ).ask()

    save_config(output_path, config_to_dict(new_config))  # type: ignore


@app.command()
def start(
    config_file: Annotated[
        str | None, cyclopts.Parameter(name=["--config", "-c"])
    ] = None,
    network: Annotated[str | None, cyclopts.Parameter(name=["--network", "-n"])] = None,
    data_folder: Annotated[
        str | None, cyclopts.Parameter(name=["--data-folder", "-d"])
    ] = None,
    seed: Annotated[int | None, cyclopts.Parameter(name=["--seed", "-s"])] = None,
    verbose: Annotated[bool, cyclopts.Parameter(name=["--verbose", "-v"])] = False,
):
    """Run the full prioritization pipeline.

    Configuration is loaded with the following priority:
    1. Command-line arguments (e.g., --seed)
    2. Values from the --config file (if provided)
    3. Default values

    It is strongly recommended to use a config file for detailed settings
    like peptide mapping, logging, and advanced options.

    Examples:
    `xlranker start -c config.yml`
    `xlranker start -c config.yml --seed 42 --verbose`
    `xlranker start --network net.tsv --data-folder data/ --seed 42`

    Args:
        config_file (str | None): path to the config file.
            If none provided, use the default configuration.
        network (str | None): path to the network file,
            and override path set in config if not None.
        data_folder (str | None): path the omic data folder
            and override path set in config if not None.
        seed (int | None): seed for random generators. If none, use random seed.
        verbose: if set, log includes more messages.

    """
    if config_file:
        try:
            load_config(config_file)
        except Exception as e:
            print(f"Error: Failed to load config file {config_file}. {e}")
            raise

    if network is not None:
        xlr_config.network_path = network

    if data_folder is not None:
        xlr_config.omic_data_folder = data_folder

    if seed is not None:
        xlr_config.seed = seed

    if xlr_config.network_path is None:
        raise ValueError("network_path not provided via --network or in config file!")
    if xlr_config.omic_data_folder is None:
        raise ValueError(
            "omic_data_folder not provided via --data-folder or in config file!"
        )
    if not xlr_config.mapping.is_fasta and xlr_config.mapping.custom_table is None:
        raise ValueError(
            "Mapping table must be provided in config (mapping.custom_table) if is_fasta is False."  # noqa: E501
        )
    xlr_config.check_is_safe()

    # Set up and run the pipeline
    setup_logging(
        verbose=verbose,
    )

    set_seed(xlr_config.seed)

    data_set = XLDataSet.load_from_config()

    _ = run_full_pipeline(data_set, threshold=xlr_config.threshold)


def cli():
    """Start the CLI."""
    app()
