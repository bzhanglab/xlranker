import logging
from xlranker.bio import Peptide, Protein

logger = logging.getLogger(__name__)


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
