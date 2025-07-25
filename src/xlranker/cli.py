import json
import logging
import os
import random
from typing import Annotated, Any

import cyclopts
import questionary
import yaml

from xlranker.config import DEFAULT_CONFIG
from xlranker.lib import XLDataSet, setup_logging
from xlranker.pipeline import run_full_pipeline
from xlranker.util import set_seed
from xlranker.util.mapping import PeptideMapper

app = cyclopts.App()
logger = logging.getLogger(__name__)


def load_config(path: str) -> dict[str, Any]:
    if path.lower().endswith(".json"):
        return json.load(open(path))
    elif path.lower().endswith(".yaml") or path.lower().endswith(".yml"):
        return yaml.safe_load(open(path))
    else:
        raise ValueError("Unsupported config file format.")


def save_config(path: str, config_obj: dict[str, Any]) -> None:
    path = path.lower()
    if path.endswith(".json"):
        return json.dump(config_obj, open(path, "w"))
    elif path.endswith(".yaml") or path.endswith(".yml"):
        return yaml.dump(config_obj, open(path, "w"))
    raise ValueError("Unsupported config file format.")


def is_folder(path_to_validate: str) -> bool | str:
    return (
        True
        if not os.path.isfile(path_to_validate)
        else "Input a directory, not an existing file."
    )


@app.command()
def init(
    default: bool = False,
    output: Annotated[str | None, cyclopts.Parameter(name=["--output", "-o"])] = None,
) -> None:
    """Initialize a config file. If no default flag provided, config is created through a interactive form.

    Args:
        default (bool, optional): Create a simple default config. Defaults to False.
        output (str | None, optional): Output config file. Can either be JSON or YAML format. Defaults to None.

    Raises:
        ValueError: raises ValueError if output is not set when default is passed

    """

    if default:
        if output is None:
            raise ValueError("Output must be specified if using default!")
        save_config(output, DEFAULT_CONFIG)
        return

    network = questionary.path(
        "Path to your peptide sequence network:",
    ).ask()
    omic_data = questionary.path(
        "Path to omic data folder:",
        only_directories=True,
        validate=is_folder,
    ).ask()
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
        case "TSV Table":
            is_fasta = False
            mapping_table_path = questionary.path("Path to TSV file:").ask()
        case _:
            is_fasta = True
            mapping_table_path = None
    if mapping_table_path is not None:
        only_human = questionary.confirm("Is your data only human data?").ask()
    else:
        only_human = True
    output_config = {
        "network": network,
        "data_folder": omic_data,
        "is_fasta": is_fasta,
        "only_human": only_human,
    }
    if mapping_table_path is not None:
        output_config["mapping_table"] = mapping_table_path
        if is_fasta:
            output_config["fasta_type"] = fasta_type
    output_path = questionary.path(
        "Output file for config (JSON or YAML format):",
        validate=lambda x: True
        if x.lower().endswith(".json")
        or x.lower().endswith(".yaml")
        or x.lower().endswith(".yml")
        else "File must end with .json, .yaml, or .yml",
    ).ask()
    save_config(output_path, output_config)


@app.command()
def test_fasta(
    fasta_file: str,
    split: str,
    gs_index: int,
    verbose: Annotated[bool, cyclopts.Parameter(name=["--verbose", "-v"])] = False,
):
    setup_logging(verbose=verbose)
    mapper = PeptideMapper(
        mapping_table_path=fasta_file, split_by=split, split_index=gs_index
    )
    sequences = ["QKTPK", "MGSGKK"]
    mapping_res = mapper.map_sequences(sequences)
    nl_char = "\n"
    for seq in sequences:
        print(f"Sequence: {seq}")
        print(f"Results:\n{nl_char.join(mapping_res.peptide_to_protein[seq])}\n")
    print("Verify results are in gene symbol!")


@app.command()
def start(
    network: Annotated[str | None, cyclopts.Parameter(name=["--network", "-n"])] = None,
    data_folder: Annotated[
        str | None, cyclopts.Parameter(name=["--data-folder", "-d"])
    ] = None,
    config_file: Annotated[
        str | None, cyclopts.Parameter(name=["--config", "-c"])
    ] = None,
    seed: Annotated[int | None, cyclopts.Parameter(name=["--seed", "-s"])] = None,
    verbose: Annotated[bool, cyclopts.Parameter(name=["--verbose", "-v"])] = False,
    log_file: Annotated[
        str | None, cyclopts.Parameter(name=["--log-file", "-l"])
    ] = None,
    mapping_table: Annotated[
        str | None, cyclopts.Parameter(name=["--mapping-table", "-m"])
    ] = None,
    split: Annotated[str | None, cyclopts.Parameter(name=["--split"])] = None,
    gs_index: Annotated[int | None, cyclopts.Parameter(name=["--gs-index"])] = None,
    is_fasta: Annotated[bool, cyclopts.Parameter(name=["--is-fasta"])] = False,
):
    """Run the full prioritization pipeline

    Requires input file to be in the format specified in the project documentation.

    Examples:

    `xlranker start network.tsv omic_data_folder/ -s 42`

    Args:
        network (Annotated[str, cyclopts.Parameter], optional): path to TSV file containing peptide network.
        data_folder (Annotated[str, cyclopts.Parameter], optional): folder containing the omics data for the model prediction.
        config_file (Annotated[ str  |  None, cyclopts.Parameter], optional): if set, read and load options from config file. Can be in JSON or YAML format.
        seed (Annotated[int  |  None, cyclopts.Parameter], optional): seed for machine learning pipeline. If not set, seed is randomly selected.
        verbose (Annotated[bool, cyclopts.Parameter], optional): enable verbose logging.
        log_file (Annotated[ str  |  None, cyclopts.Parameter], optional): if set, saves logging to path
        mapping_table (Annotated[ str  |  None, cyclopts.Parameter], optional): path to custom mapping table for peptide sequences
        split (Annotated[ str  |  None, cyclopts.Parameter], optional): character used for splitting the FASTA file header
        gs_index (Annotated[int  |  None, cyclopts.Parameter], optional): index in the FASTA file that contains the gene symbol. Index starts at 0.
        is_fasta (Annotated[bool, cyclopts.Parameter], optional): Enable if mapping table is a FASTA file.

    """

    if config_file is not None:
        config_data = load_config(config_file)
    else:
        config_data = DEFAULT_CONFIG

    # Check if network and data_folder are set, which are required

    # Use CLI arg if provided, otherwise fall back to config
    network = network if network is not None else config_data.get("network", None)
    data_folder = (
        data_folder if data_folder is not None else config_data.get("data_folder", None)
    )
    if network is None:
        raise ValueError("network not provided in command or in config!")
    if data_folder is None:
        raise ValueError("data_folder not provided in command or in config!")
    seed = seed if seed is not None else config_data.get("seed", None)
    verbose = verbose or config_data.get("verbose", False)
    log_file = log_file or config_data.get("log_file", None)
    mapping_table = mapping_table or config_data.get("mapping_table", None)
    split = split or config_data.get("split", None)
    gs_index = gs_index if gs_index is not None else config_data.get("gs_index", None)
    is_fasta = is_fasta or config_data.get("is_fasta", False)

    setup_logging(verbose=verbose, log_file=log_file)
    if seed is None:
        seed = random.randint(0, 10000000)
        logger.info(f"Randomly generated seed: {seed}")

    set_seed(seed)

    data_set = XLDataSet.load_from_network(
        network,
        data_folder,
        custom_mapping_path=mapping_table,
        is_fasta=is_fasta,
        split_by=split,
        split_index=gs_index,
    )

    # run the full pipeline
    _ = run_full_pipeline(data_set)


def cli():
    app()
