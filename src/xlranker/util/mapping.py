from xlranker.data import get_gencode_fasta
from Bio import SeqIO
import logging

from xlranker.util.readers import read_mapping_table_file

logger = logging.getLogger(__name__)


class PeptideMapper:
    mapping_table_path: str
    split_by: str
    split_index: int
    is_fasta: bool

    def __init__(
        self,
        mapping_table_path: str | None = None,
        split_by: str = "|",
        split_index: int = 6,
        is_fasta: bool = True,
    ):
        if mapping_table_path is None:
            logger.info("Using default gencode fasta file for peptide mapping")
            self.mapping_table_path = get_gencode_fasta()
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
        for seq in sequences:
            if seq not in map_res:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.error(f"{seq} not found in mapping table!")
                had_error = True
            elif len(map_res[seq]) == 0:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.error(f"{seq} mapped to no proteins!")
                had_error = True
        if had_error:
            raise ValueError(
                "Mapping incomplete. Verify peptide sequences and mapping table is correct."
            )
        return map_res

    def map_fasta(self, sequences: list[str]) -> dict[str, list[str]]:
        matches = {}
        logger.info(f"Mapping {len(sequences)} peptide sequences")
        for record in SeqIO.parse(self.mapping_table_path, "fasta"):
            for sequence in sequences:
                if sequence in record.seq:
                    if sequence not in matches:
                        matches[sequence] = set()
                    matches[sequence].add(
                        record.description.split(self.split_by)[
                            self.split_index
                        ].strip()
                    )
        for key in matches:
            matches[key] = list(matches[key])
        return matches
