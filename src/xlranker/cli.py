import logging
import random
import cyclopts
from xlranker.parsimony.selection import ParsimonySelector
from xlranker.util.mapping import PeptideMapper
import xlranker.ml.data as xlr_data
from xlranker.lib import XLDataSet, setup_logging
import xlranker.config
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
    sequences = ["QKTPK", "MGSGKK"]
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
    print(mapper.map_sequences(["VVEPKR", "MGSGKK"]))


@app.command()
def test_loading():
    df = xlr_data.load_default_ppi()
    print(df.head())


@app.command()
def test_prioritization(network: str, omic_folder: str):
    setup_logging()
    logger.info("Loading data")
    dataset = XLDataSet.load_from_network(network, omic_folder)
    logging.info("Building protein data")
    dataset.build_proteins()
    prioritizer = ParsimonySelector(dataset)
    logging.info("Running prioritization")
    prioritizer.run()
    logging.info("Saving results")
    tsv_data = ["pair_id\tstatus\tgroup_id"]
    for protein_pair_groups in prioritizer.protein_groups.values():
        for protein_pair in protein_pair_groups:
            tsv_data.append(protein_pair.to_tsv())
    with open("output.tsv", "w") as w:
        w.write("\n".join(tsv_data))


@app.command()
def start(
    network: Annotated[str, cyclopts.Parameter(name=["--network", "-n"])],
    data_folder: Annotated[str, cyclopts.Parameter(name=["--data-folder", "-d"])],
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
    config_file: Annotated[str | None, cyclopts.Parameter(name=["--config"])] = None,
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
        mapping_table (Annotated[ str  |  None, cyclopts.Parameter, optional): path to custom mapping table for peptide sequences
        split (Annotated[ str  |  None, cyclopts.Parameter, optional): character used for splitting the FASTA file header
        gs_index (Annotated[int  |  None, cyclopts.Parameter, optional): index in the FASTA file that contains the gene symbol. Index starts at 0.
        is_fasta (Annotated[bool, cyclopts.Parameter, optional): Enable if mapping table is a FASTA file.
        config (Annotated[str | None, cyclopts.Parameter, optional): Path of JSON file containing config parameters
    """
    if config_file is not None:
        xlranker.config.load_from_json(config_file)
    setup_logging(verbose=verbose, log_file=log_file)
    if seed is None:
        seed = random.randint(0, 10000000)
        logger.info(f"Randomly generated seed: {seed}")
    _data_set = XLDataSet.load_from_network(
        network,
        data_folder,
        custom_mapping_path=mapping_table,
        is_fasta=is_fasta,
        split_by=split,
        split_index=gs_index,
    )


def cli():
    app()
