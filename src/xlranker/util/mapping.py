from xlranker.data import get_gencode_fasta
from Bio import SeqIO
import logging
from enum import Enum, auto

from xlranker.util.readers import read_mapping_table_file

logger = logging.getLogger(__name__)


class FastaType(Enum):
    UNIPROT = auto()
    GENCODE = auto()


def extract_gene_symbol_uniprot(fasta_description: str) -> str:
    splits = fasta_description.split(" ")
    for split in splits:
        if "GN=" in split:  # check if gene name split
            return split[3:]  # Remove GN= from string
    return fasta_description  # return if failed


def extract_gene_symbol(fasta_description: str, fasta_type: FastaType) -> str:
    match fasta_type:
        case FastaType.UNIPROT:
            return extract_gene_symbol_uniprot(fasta_description)
        case FastaType.GENCODE:
            return ""


class PeptideMapper:
    """Peptide mapper class

    Raises:
        ValueError: Raises error if there is an issue with mapping tables
    """

    mapping_table_path: str
    split_by: str
    split_index: int
    is_fasta: bool

    def __init__(
        self,
        mapping_table_path: str | None = None,
        split_by: str = "|",
        split_index: int = 3,
        is_fasta: bool = True,
    ):
        """initializes PeptideMapper

        Args:
            mapping_table_path (str | None, optional): path to mapping table. Can be in fasta or mapping table. If none, then uses the default gencode v42 version. Defaults to None.
            split_by (str, optional): character in fasta description to split into id components. Defaults to "|".
            split_index (int, optional): index of gene symbol in fasta file. Defaults to 3.
            is_fasta (bool, optional): is input file fasta file. Defaults to True.
        """
        if mapping_table_path is None:
            logger.info("Using default gencode fasta file for peptide mapping")
            self.mapping_table_path = get_gencode_fasta()
            # Make sure variables match defaults
            split_by = "|"
            split_index = 3
            is_fasta = True
        else:
            logger.info("Using custom fasta file for peptide mapping")
            logging.debug(f"FASTA File Path: {mapping_table_path}")
            self.mapping_table_path = mapping_table_path
        self.split_by = split_by
        self.split_index = split_index
        self.is_fasta = is_fasta

    def map_sequences(self, sequences: list[str]) -> dict[str, list[str]]:
        map_res = dict()
        if self.is_fasta:  # determine which mapping function to use
            map_res = self.map_fasta(sequences)
        else:  # mapping table just needs to be read
            map_res = read_mapping_table_file(self.mapping_table_path)
        had_error = False
        for seq in sequences:  # verify all sequences have mapping information
            if seq not in map_res:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.error(f"{seq} not found in mapping table!")
                had_error = True
            elif len(map_res[seq]) == 0:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.error(f"{seq} maps to no proteins!")
                had_error = True
        if had_error:
            raise ValueError(
                "Mapping incomplete. Verify peptide sequences and mapping table is correct."
            )
        return map_res

    def map_fasta(self, sequences: list[str]) -> dict[str, list[str]]:
        matches: dict[str, set[str]] = {}
        logger.info(f"Mapping {len(sequences)} peptide sequences")
        for record in SeqIO.parse(self.mapping_table_path, "fasta"):
            for sequence in sequences:
                if sequence in record.seq:
                    if sequence not in matches:
                        matches[sequence] = set()
                    split_res = record.description.split(self.split_by)
                    if self.split_index >= len(split_res):
                        matches[sequence].add(split_res[0])
                    elif len(split_res) != 0:
                        matches[sequence].add(split_res[self.split_index].split(" ")[0])

        final_matches: dict[str, list[str]] = {}
        for key in matches:
            final_matches[key] = list(matches[key])
        return final_matches
