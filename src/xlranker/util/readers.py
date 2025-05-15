import logging
import polars as pl
from xlranker.bio import Peptide
from xlranker.bio import PeptideGroup
from pathlib import Path

logger = logging.getLogger(__name__)


def read_data_matrix(data_path: str, additional_null_values=[]) -> pl.DataFrame:
    """reads data matrix into a Polars DataFrame

    Format:
     - Has Header (any names allowed)
     - null/missing values: "", "NA". More can be added.

    Args:
        data_path (str): path to the data matrix
        additional_null_values (list[str]): list of str of additional values that should considered as null

    Returns:
        pl.DataFrame: Polars DataFrame of the input data
    """
    null_values = ["", "NA"]
    null_values.extend(additional_null_values)
    return pl.read_csv(
        data_path, has_header=True, separator="\t", null_values=null_values
    )


def read_data_folder(folder_path: str, additional_null_values=[]) -> list[pl.DataFrame]:
    """reads all TSV files in a folder

    Args:
        folder_path (str): path of the folder that contains files ending in .tsv
        additional_null_values (list[str]): list of str of additional values that should considered as null in the data files

    Raises:
        FileNotFoundError: raised if no TSV files are found

    Returns:
        list[pl.DataFrame]: list of all of the data files in a Polars DataFrame, as read by the read_data_matrix function
    """
    file_list = Path(folder_path).glob("*.tsv")
    if len(file_list) == 0:
        raise FileNotFoundError(f"No TSV files were found in directory: {folder_path}")
    return [
        read_data_matrix(file, additional_null_values=additional_null_values)
        for file in file_list
    ]


def read_network_file(network_path: str) -> list[PeptideGroup]:
    """reads TSV network file to a list of PeptideGroup

    Args:
        network_path (str): path to the TSV file

    Returns:
        list[PeptideGroup]: list of PeptideGroup representing the network
    """
    try:
        with open(network_path) as r:
            text = r.read().split("\n")
        new_rows = set()  # Track unique rows
        valid_rows = 0  # Keeps track of number of edges in original file
        for row in text:
            if "\t" in row:
                valid_rows += 1
                vals = row.split("\t")
                val_a = vals[0]
                val_b = vals[1]
                if val_a > val_b:  # Make sure edges are all sorted the same.
                    temp = val_a
                    val_a = val_b
                    val_b = temp
                new_rows.add(f"{val_a}\t{val_b}")
    except IndexError:
        logger.error("Index out of bound. Make sure network is in the correct format.")
        raise IndexError()
    except FileNotFoundError:
        logger.error(f"File not found: {network_path}")
        raise FileNotFoundError
    duplicate_rows = valid_rows - len(new_rows)  # Count number of duplicated rows
    if duplicate_rows > 0:  # Send warning that duplicate edges were removed.
        logger.warning(
            f"Found and removed {duplicate_rows} duplicated edge(s) in network."
        )
    network = []
    for row in new_rows:
        vals = row.split("\t")
        a = Peptide(vals[0])
        b = Peptide(vals[1])
        group = PeptideGroup(a, b)
        network.append(group)
    return network


def read_mapping_table_file(file_path: str) -> dict[str, list[str]]:
    """read mapping file where the first column is the peptide sequence and the following columns are proteins that map to that sequence

    Args:
        file_path (str): path to the tab-separated mapping table
    """
    try:
        with open(file_path, mode="r") as r:
            raw_text = r.read().split("\n")
        mapping_res: dict[str, list[str]] = dict()
        uniq_sequences: set[str] = set()
        for line in raw_text:
            if "\t" in line:
                vals = line.split("\t")
                seq = vals[0]
                if seq in uniq_sequences:
                    logging.warning(
                        f"Peptide sequence {seq} duplicated! Keeping first instance."
                    )
                else:
                    uniq_sequences.add(seq)
                    mapping_res[seq] = vals[1:]
        if len(mapping_res) == 0:
            logging.error(f"No peptide sequences found in mapping file: {file_path}")
            raise ValueError("No peptide sequence identified")
        return mapping_res
    except FileNotFoundError:
        logging.error(f"Could not find mapping table file at {file_path}!")
        raise ValueError("Could not read mapping table: File not found.")
