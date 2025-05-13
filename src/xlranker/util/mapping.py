from xlranker.data import get_gencode_fasta
from Bio import SeqIO
import logging

logger = logging.getLogger()


class PeptideMapper:
    fasta_path: str

    def __init__(self, fasta_path: str | None = None):
        if fasta_path is None:
            logger.info("Using default gencode fasta file for peptide mapping")
            self.fasta_path = get_gencode_fasta()
        else:
            logger.info("Using custom fasta file for peptide mapping")
            self.fasta_path = fasta_path

    def map_sequences(self, sequences: str) -> dict[str, list[str]]:
        matches = {}
        logger.info(f"Mapping {len(sequences)} peptide sequences")
        for record in SeqIO.parse(self.fasta_path, "fasta"):
            for sequence in sequences:
                if sequence in record.seq:
                    if sequence not in matches:
                        matches[sequence] = []
                    matches[sequence].append(record.description.split("|")[-2])
        return matches
