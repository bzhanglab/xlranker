import logging
import random
import cyclopts
from xlranker.util.mapping import PeptideMapper
import xlranker.ml.data as xlr_data
from xlranker.lib import setup_logging
from typing import Annotated


app = cyclopts.App()
logger = logging.getLogger(__name__)


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
    sequences = ["SGGLSNL", "MIYD", "NGLEEKRKS"]
    mapping_res = mapper.map_sequences(sequences)
    for seq in sequences:
        print(f"Sequence: {seq}")
        print(f"Results:\n{'\n'.join(mapping_res[seq])}\n")
    print("Verify results are in gene symbol!")


@app.command()
def test(
    verbose: Annotated[bool, cyclopts.Parameter(name=["--verbose", "-v"])] = False,
):
    setup_logging(verbose)
    mapper = PeptideMapper()
    print(mapper.map_sequences(["SGGLSNL", "MIYD", "NGLEEKRKS"]))


@app.command()
def test_loading():
    df = xlr_data.load_default_ppi()
    print(df.head())


@app.command()
def start(
    network: Annotated[str, cyclopts.Parameter(name=["--network", "-n"])],
    data_folder: Annotated[str, cyclopts.Parameter(name=["--data-folder", "-d"])],
    seed: Annotated[int | None, cyclopts.Parameter(name=["--seed", "-s"])] = None,
    verbose: Annotated[bool, cyclopts.Parameter(name=["--verbose", "-v"])] = False,
    log_file: Annotated[
        str | None, cyclopts.Parameter(name=["--log-file", "-l"])
    ] = None,
):
    """Run the full prioritization pipeline

    Requires input file to be in the format specified in the project documentation.

    Examples:

    `xlranker start network.tsv omic_data_folder/ -s 42`

    Args:
        network (Annotated[str, cyclopts.Parameter, optional): path to TSV file containing peptide network.
        data_folder (Annotated[str, cyclopts.Parameter, optional): folder containing the omics data for the model prediction.
        seed (Annotated[int  |  None, cyclopts.Parameter, optional): seed for machine learning pipeline. If not set, seed is randomly selected.
        verbose (Annotated[bool, cyclopts.Parameter, optional): enable verbose logging.
        log_file (Annotated[ str  |  None, cyclopts.Parameter, optional): if set, saves logging to path
    """
    setup_logging(verbose=verbose, log_file=log_file)
    if seed is None:
        seed = random.randint(0, 10000000)
        logger.info(f"Randomly generated seed: {seed}")


def cli():
    app()
