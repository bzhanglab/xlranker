"""Functions for reading data files, mapping files, and networks."""

import logging
import re
from pathlib import Path

import polars as pl

from xlranker.bio import Peptide
from xlranker.bio.pairs import PeptidePair
from xlranker.config import config
from xlranker.util import get_pair_id

logger = logging.getLogger(__name__)


def parse_linkage(link: str) -> int | None:
    """Parse a link from a linkage string by removing all non-numeric chars.

    Args:
        link (str): the original linkage string

    Returns:
        int | None: 1-based index of the link or None if unable to parse
    """
    cleaned_link = re.sub("[^0-9]", "", link)
    if len(cleaned_link) > 0:
        return int(cleaned_link)
    return None


def read_data_matrix(
    data_path: str, additional_null_values: list[str] = []
) -> pl.DataFrame:
    """Read data matrix into a Polars DataFrame with samples/measurements being columns.

    Format:
     - Has header (any names allowed).
     - First column must be the protein/gene followed by measurements.
     - Null/missing values: "", "NA". More can be added.

    Args:
        data_path (str): path to the data matrix
        additional_null_values (list[str]): list of str of additional values that should
            considered as null

    Returns:
        pl.DataFrame: Polars DataFrame of the input data

    """
    null_values = ["", "NA", "None"]
    null_values.extend(additional_null_values)
    with open(data_path) as f:
        header = f.readline().strip().split("\t")

    # Force all except first column to Float64
    dtype_overrides: dict[str, pl.Float64 | pl.Utf8] = {
        col: pl.Float64 for col in header[1:]
    }
    # Keep first column as string
    dtype_overrides[header[0]] = pl.Utf8

    return pl.read_csv(
        data_path,
        has_header=True,
        separator="\t",
        null_values=null_values,
        schema_overrides=dtype_overrides,
    )


def base_name(file_path: Path | str) -> str:
    """Get the base name from a path.

    Example:
        ```
        file_path = "example/test.tsv"
        base = base_name(file_path) # base = "test"
        assert base == "test"
        ```

    Args:
        file_path (Path | str): path of file to get base name of

    Returns:
        str: the base file name

    """
    return Path(file_path).stem


def read_data_folder(
    folder_path: str, omic_data_files: list[str] | None
) -> dict[str, pl.DataFrame]:
    """Reads all TSV files in a folder.

    Args:
        folder_path (str): path of the folder that contains files ending in .tsv
        omic_data_files (list[str] | None): list of files to read from the folder
            If none, use all files in the folder.
            Otherwise, use only the files in the list.

    Returns:
        list[pl.DataFrame]: list of all of the data files in a Polars DataFrame, as read
            by the read_data_matrix function

    """
    if omic_data_files is None:
        file_glob = Path(folder_path).glob("*.tsv")
        file_list: list[Path] = list(file_glob)
    else:
        file_list = [Path(file) for file in omic_data_files]
    if len(file_list) == 0:
        logger.warning(f"No TSV files were found in directory: {folder_path}")
    ret_dict = {}
    for file in file_list:
        ret_dict[base_name(file)] = read_data_matrix(
            str(file), additional_null_values=config.additional_null_values
        )
    return ret_dict


def read_network_file(network_path: str) -> tuple[dict[str, PeptidePair], bool]:
    """Reads TSV network file to a list of PeptideGroup.

    Args:
        network_path (str): path to the TSV file

    Returns:
        tuple[list[PeptideGroup], bool]: tuple with first element
            is list of PeptideGroup representing the network
            the bool is True if linkage information is present

    Raises:
        IndexError: Raised if there are not 2 columns in the network.
        FileNotFoundError: Raised if no file is found in network_path

    """
    try:
        with open(network_path) as r:
            text = r.read().split("\n")
        new_rows = set()  # Track unique rows
        valid_rows = 0  # Keeps track of number of edges in original file
        has_linkages = False
        for row in text:
            if "\t" in row:
                valid_rows += 1
                vals = row.split("\t")
                if len(vals) == 2:
                    val_a = vals[0]
                    val_b = vals[1]
                    if val_a > val_b:  # Make sure edges are all sorted the same.
                        val_a, val_b = val_b, val_a
                    new_rows.add(f"{val_a}\t{val_b}")
                elif len(vals) == 4:
                    has_linkages = True
                    val_a = vals[0]
                    link_a = vals[1]
                    val_b = vals[2]
                    link_b = vals[3]
                    if val_a > val_b:  # Make sure edges are all sorted the same.
                        val_a, val_b = val_b, val_a
                        link_a, link_b = link_b, link_a
                    new_rows.add(f"{val_a}\t{val_b}\t{link_a}\t{link_b}")
                else:
                    raise ValueError(
                        """Network must contain either 2 or 4 columns.
                        See documentation for details."""
                    )

    except IndexError as e:
        logger.error("Index out of bound. Make sure network is in the correct format.")
        raise IndexError() from e
    except FileNotFoundError as e:
        logger.error(f"File not found: {network_path}")
        raise FileNotFoundError from e
    duplicate_rows = valid_rows - len(new_rows)  # Count number of duplicated rows
    if duplicate_rows > 0:  # Send warning that duplicate edges were removed.
        logger.info(
            f"Found and removed {duplicate_rows} duplicated edge(s) in network."
        )
    network: dict[str, PeptidePair] = {}
    for row in new_rows:
        vals = row.split("\t")
        link_a = None
        link_b = None
        if has_linkages:
            link_a = parse_linkage(vals[2])
            link_b = parse_linkage(vals[3])
        a = Peptide(vals[0], linkage=link_a)
        b = Peptide(vals[1], linkage=link_b)
        group = PeptidePair(a, b)
        network[get_pair_id(a, b)] = group
    return (network, has_linkages)


def read_mapping_table_file(file_path: str) -> dict[str, list[str]]:
    """Read mapping file into dict.

    Format:
        The first column is the peptide sequence and the following columns are proteins
        that map to that sequence.

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
