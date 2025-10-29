"""Report helper functions."""

from pathlib import Path

from xlranker.util import get_pair_id_from_str
from xlranker.bio.pairs import ProteinPair
from xlranker.config import config
from xlranker.lib import write_pair_to_network, XLDataSet
from xlranker.status import ReportStatus
from pydantic import BaseModel
import csv


class DetailedRow(BaseModel):
    """Row for the detailed report."""

    peptide_pair: str
    protein_pair: str
    mapping_order: str
    protein_pair_status: str
    selection_criteria: str
    pair_score: float


def make_report(
    pairs: list[ProteinPair], status: ReportStatus, output_path: Path
) -> None:
    """Write network of protein pairs that have the level of the inputted status or higher.

    Status order (highest on top):
        1. Conservative
        2. Minimal
        3. Expanded
        4. All

    Args:
        pairs (list[ProteinPair]): list of all protein pairs
        status (ReportStatus): Minimum status to accept.
        output_path (Path): output path to save pair to

    """
    valid_pairs = [pair for pair in pairs if pair.report_status <= status]
    write_pair_to_network(valid_pairs, str(output_path))


def make_all_reports(pairs: list[ProteinPair]) -> None:
    """Make reports for all status levels into the output_folder specified in the config.

    Args:
        pairs (list[ProteinPair]): list of all protein pairs.

    """
    output_folder = Path(config.output).joinpath("reports")
    output_folder.mkdir(exist_ok=True)
    make_report(pairs, ReportStatus.CONSERVATIVE, output_folder / "conservative.tsv")
    make_report(pairs, ReportStatus.MINIMAL, output_folder / "minimal.tsv")
    make_report(pairs, ReportStatus.EXPANDED, output_folder / "expanded.tsv")
    make_report(pairs, ReportStatus.ALL, output_folder / "all.tsv")


def make_mapping_tables(data_set: XLDataSet) -> None:
    """Make all mapping tables including peptide pairs to protein pairs and peptide to protein.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
    """
    mapping_peptide_pair_to_protein_pairs(data_set)
    mapping_peptide_to_protein(data_set)


def mapping_peptide_pair_to_protein_pairs(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Generate mapping table that shows peptide pairs and their corresponding mapped protein pairs.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path. Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_protein_pair.tsv")
    data = []
    for peptide_pair in sorted(data_set.peptide_pairs.keys()):
        row = [peptide_pair]
        for protein_pair in data_set.peptide_pairs[peptide_pair].connections:
            row.append(protein_pair)
        data.append("\t".join(row))
    with open(output_path, "w") as w:
        w.write("\n".join(data))


def mapping_peptide_to_protein(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Generate mapping table that shows a peptide and its corresponding mapped proteins.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path. Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_protein.tsv")
    data = []
    uniq_peptides = set()
    for peptide_pair in sorted(
        data_set.peptide_pairs.values(), key=lambda p: p.a.sequence + p.b.sequence
    ):
        if peptide_pair.a.sequence not in uniq_peptides:
            uniq_peptides.add(peptide_pair.a.sequence)
            row = [peptide_pair.a.sequence]
            for prot in peptide_pair.a.mapped_proteins:
                row.append(prot)
            data.append("\t".join(row))
        if peptide_pair.b.sequence not in uniq_peptides:
            uniq_peptides.add(peptide_pair.b.sequence)
            row = [peptide_pair.b.sequence]
            for prot in peptide_pair.b.mapped_proteins:
                row.append(prot)
            data.append("\t".join(row))
    with open(output_path, "w") as w:
        w.write("\n".join(data))


def detailed_report(data_set: XLDataSet, output_path: str | None = None) -> None:
    """Generates final detailed report.

    Args:
        data_set (XLDataSet): XL Dataset to make detailed_report of
        output_path (str): final output
    """
    rows: list[DetailedRow] = []
    for pair_id in sorted(data_set.peptide_pairs.keys()):
        peptide_pair = data_set.peptide_pairs[pair_id]
        for mapped_prot_a in sorted(peptide_pair.a.mapped_proteins):
            for mapped_prot_b in sorted(peptide_pair.b.mapped_proteins):
                prot_pair = data_set.protein_pairs.get(
                    get_pair_id_from_str(mapped_prot_a, mapped_prot_b), None
                )
                if prot_pair is None:
                    continue
                row_val = DetailedRow(
                    peptide_pair=pair_id,
                    protein_pair=prot_pair.pair_id,
                    mapping_order=f"{mapped_prot_a}{config.advanced.pair_separator}{mapped_prot_b}",
                    protein_pair_status=str(
                        prot_pair.report_status,
                    ),
                    selection_criteria=str(prot_pair.prioritization_status),
                    pair_score=float(prot_pair.score),
                )
                rows.append(row_val)
    if output_path is None:
        output_folder = Path(config.output)
        output_path = str(output_folder / "detailed_report.tsv")
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=DetailedRow.model_fields.keys(), delimiter="\t"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row.model_dump())
