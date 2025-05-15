import cyclopts
from xlranker.util.mapping import PeptideMapper
import xlranker.ml.data as xlr_data
from xlranker.lib import setup_logging
from typing import Annotated


app = cyclopts.App()


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
        print(f"Results:\n{"\n".join(mapping_res[seq])}\n")
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
def start(input: str):
    """Run the full prioritization pipeline

    Requires input file to be in the format specified in the project documentation.

    Args:
        input (str): input file in TSV format
    """
    pass


def cli():
    app()
