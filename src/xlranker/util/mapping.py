"""Mapping related classes and functions."""

import logging
from enum import Enum, auto

from Bio import SeqIO

from xlranker.data import get_gencode_fasta
from xlranker.util.readers import read_mapping_table_file

logger = logging.getLogger(__name__)


class FastaType(Enum):
    """Types of Fasta files supported by XLRanker."""

    UNIPROT = auto(), "UNIPROT FASTA type"
    GENCODE = auto(), "Gencode FASTA type"


def extract_gene_symbol_uniprot(fasta_description: str) -> str:
    """Get the gene symbol from a UNIPROT style FASTA description.

    Method:
        1. Split the description by spaces
        2. Find split with GN= (Gene Name)
        3. Remove GN= from split and return

    If split with GN= not found, return the UNIPROT symbol.
        1. Using first split (when splitting by space), split again by |
        2. If there is at least 2 elements in split, return second element

    If can't get UNIPROT symbol, return original description.

    Args:
        fasta_description (str): FASTA description string

    Returns:
        str: Gene Symbol from description. If can't be extracted, try getting UNIPROT ID.
             If all fails, return original description
    """
    splits = fasta_description.split(" ")
    for split in splits:
        if "GN=" in split:  # check if gene name split
            return split[3:]  # Remove GN= from string
    splits = splits[0].split("|")
    if len(splits) >= 2:
        return splits[1]
    return fasta_description  # return if failed


def extract_gene_symbol_gencode(fasta_description: str, **kwargs) -> str:
    """Get the gene symbol from a UNIPROT style FASTA description.

    Method:
        1. Split the description by spaces
        2. Find split with GN= (Gene Name)
        3. Remove GN= from split and return

    If split with GN= not found, return the UNIPROT symbol.
        1. Using first split (when splitting by space), split again by |
        2. If there is at least 2 elements in split, return second element

    If can't get UNIPROT symbol, return original description.

    Args:
        fasta_description (str): FASTA description string
        split_by (str): Character to split description string
        split_index (str): Index (0-based) of gene symbol after splitting.
                           All characters after first space are removed.

    Returns:
        str: Gene Symbol from description. If can't be extracted, return original description

    """
    split_by = kwargs["split_by"]
    split_index = kwargs["split_index"]
    split_res = fasta_description.split(split_by)
    if split_index >= len(split_res):
        return split_res[0]  # keep first split if split_index is too large
    if len(split_res) != 0:
        return split_res[split_index].split(" ")[0]  # remove elements after space
    return fasta_description  # return if failed


def extract_gene_symbol(fasta_description: str, fasta_type: FastaType, **kwargs) -> str:
    match fasta_type:
        case FastaType.UNIPROT:
            return extract_gene_symbol_uniprot(fasta_description).upper()
        case FastaType.GENCODE:
            return extract_gene_symbol_gencode(fasta_description, **kwargs).upper()


class PeptideMapper:
    """Peptide mapper class.

    Raises:
        ValueError: Raises error if there is an issue with mapping tables

    """

    mapping_table_path: str
    split_by: str
    split_index: int
    is_fasta: bool
    fasta_type: FastaType

    def __init__(
        self,
        mapping_table_path: str | None = None,
        split_by: str = "|",
        split_index: int = 3,
        is_fasta: bool = True,
        fasta_type: FastaType = FastaType.UNIPROT,
    ) -> None:
        """Initialize PeptideMapper.

        Args:
            mapping_table_path (str | None, optional): Path to mapping table.
                                                       Can be in fasta or mapping table.
                                                       If none, then uses the default uniprot version
                                                       Defaults to None.
            split_by (str, optional): character in fasta description to split into id components.
                                      Defaults to "|".
            split_index (int, optional): index of gene symbol in fasta file. Defaults to 3.
            is_fasta (bool, optional): is input file fasta file. Defaults to True.
            fasta_type (FastaType): Type of FASTA header. Can be UNIPROT or GENCODE

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
        self.fasta_type = fasta_type

    def map_sequences(self, sequences: list[str]) -> dict[str, list[str]]:
        """Map a list of sequences to genes.

        Args:
            sequences (list[str]): list of sequences to map to genes

        Returns:
            dict[str, list[str]]: dictionary where keys are peptide sequences
                                  values are list of genes that map to that sequence

        """
        map_res = dict()
        if self.is_fasta:  # determine which mapping function to use
            map_res = self.map_fasta(sequences)
        else:  # mapping table just needs to be read
            map_res = read_mapping_table_file(self.mapping_table_path)
        no_maps = 0
        for seq in sequences:  # verify all sequences have mapping information
            if seq not in map_res:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.warning(f"{seq} not found in mapping table!")
            elif len(map_res[seq]) == 0:
                logger.debug(f"is_fasta: {self.is_fasta}")
                logger.warning(f"{seq} maps to no proteins!")
                no_maps += 1
        if no_maps != 0:
            logger.warning(f"{no_maps} sequences do not have mapped proteins")
        return map_res

    def map_fasta(self, sequences: list[str]) -> dict[str, list[str]]:
        matches: dict[str, set[str]] = {}
        for seq in sequences:
            matches[seq] = set()
        logger.info(f"Mapping {len(sequences)} peptide sequences")
        for record in SeqIO.parse(self.mapping_table_path, "fasta"):
            for sequence in sequences:
                if sequence in record.seq:
                    matches[sequence].add(
                        extract_gene_symbol(
                            record.description,
                            self.fasta_type,
                            split_by=self.split_by,
                            split_index=self.split_index,
                        )
                    )

        final_matches: dict[str, list[str]] = {}
        for key in matches:
            final_matches[key] = list(matches[key])
        return final_matches
