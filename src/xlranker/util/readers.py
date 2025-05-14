import logging
from xlranker.bio import Peptide, Protein
from xlranker.bio import PeptideGroup

logger = logging.getLogger(__name__)


def read_network_file(network_path: str) -> list[PeptideGroup]:
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


def read_mapping_table_file(file_path: str) -> list[Peptide]:
    """read mapping file where the first column is the peptide sequence and the following columns are proteins that map to that sequence

    Args:
        file_path (str): path to the tab-separated mapping table
    """
    try:
        with open(file_path, mode="r") as r:
            raw_text = r.read().split("\n")
        peptides: list[Peptide] = list()
        uniq_sequences: set[str] = set()
        for line in raw_text:
            if "\t" in line:
                vals = line.split("\t")
                seq = vals[0]
                if seq in uniq_sequences:
                    logging.warning(
                        f"Peptide Sequence {seq} duplicated! Keeping first instance."
                    )
                else:
                    uniq_sequences.add(seq)
                    proteins: list[Protein] = []
                    for val in vals[1:]:
                        new_prot = Protein(val)
                        proteins.append(new_prot)
                    peptide = Peptide(sequence=seq, mapped_proteins=proteins)
                    peptides.append(peptide)
        if len(peptides) == 0:
            logging.error(f"No peptide sequences found in mapping file: {file_path}")
            raise ValueError("No peptide sequence identified")
        return peptides
    except FileNotFoundError:
        logging.error(f"Could not find mapping table file at {file_path}!")
        raise ValueError("Could not read mapping table: File not found.")
