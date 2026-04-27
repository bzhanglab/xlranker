"""Report helper functions."""

import csv
from pathlib import Path

from pydantic import BaseModel

from xlranker.bio.pairs import ProteinPair
from xlranker.config import config
from xlranker.lib import (
    XLDataSet,
    write_pair_to_network,
    write_pair_to_network_with_linkages,
)
from xlranker.status import ReportStatus
from xlranker.util import get_pair_id_from_str


class DetailedRow(BaseModel):
    """Row for the detailed report."""

    peptide_pair: str
    protein_pair: str
    mapping_order: str
    peptide_linkage_pair: str = ""
    protein_linkage_pair: str = ""
    protein_pair_status: str
    selection_criteria: str
    pair_score: float


def make_report(
    pairs: list[ProteinPair], status: ReportStatus, output_path: Path
) -> None:
    """Write network of protein pairs that have the level of the input status or higher.

    Status order (highest on top):
        1. Unique
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


def make_report_with_linkages(
    pairs: list[ProteinPair], status: ReportStatus, output_path: Path
) -> None:
    """Write network of protein pairs with linkages.

    Args:
        pairs (list[ProteinPair]): list of all protein pairs
        status (ReportStatus): Minimum status to accept.
        output_path (Path): output path to save pair to

    """
    valid_pairs = [pair for pair in pairs if pair.report_status <= status]
    write_pair_to_network_with_linkages(valid_pairs, str(output_path))


def make_all_reports(data_set: XLDataSet) -> None:
    """Write reports for all status levels to the output_folder specified in the config.

    Args:
        data_set (XLDataSet): data set containing protein pairs and linkage metadata.

    """
    pairs = list(data_set.protein_pairs.values())
    output_folder = Path(config.output).joinpath("reports")
    output_folder.mkdir(exist_ok=True)
    make_report(pairs, ReportStatus.UNIQUE, output_folder / "unique.tsv")
    make_report(pairs, ReportStatus.MINIMAL, output_folder / "minimal.tsv")
    make_report(pairs, ReportStatus.EXPANDED, output_folder / "expanded.tsv")
    make_report(pairs, ReportStatus.ALL, output_folder / "all.tsv")

    if data_set.has_linkages:
        make_report_with_linkages(
            pairs, ReportStatus.UNIQUE, output_folder / "unique_with_linkages.tsv"
        )
        make_report_with_linkages(
            pairs, ReportStatus.MINIMAL, output_folder / "minimal_with_linkages.tsv"
        )
        make_report_with_linkages(
            pairs, ReportStatus.EXPANDED, output_folder / "expanded_with_linkages.tsv"
        )
        make_report_with_linkages(
            pairs, ReportStatus.ALL, output_folder / "all_with_linkages.tsv"
        )


def make_mapping_tables(data_set: XLDataSet) -> None:
    """Make all mapping tables.

    Includes peptide pairs to protein pairs and peptide to protein.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
    """
    mapping_peptide_pair_to_protein_pairs(data_set)
    mapping_peptide_to_protein(data_set)
    mapping_peptide_pair_to_best_protein_pair(data_set)
    if data_set.has_linkages:
        mapping_report(data_set)
        mapping_peptide_pair_to_protein_pairs_with_linkages(data_set)
        mapping_peptide_to_protein_with_linkages(data_set)


def mapping_report(data_set: XLDataSet, output_path: str | None = None) -> None:
    """Write detailed mapping report.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "mapping_report.tsv")

    header = [
        "pair_id",
        "peptide a",
        "pep_residue_a",
        "start_a",
        "end_a",
        "peptide b",
        "pep_residue_b",
        "start_b",
        "end_b",
        "gene_a",
        "residue_a",
        "gene_b",
        "residue_b",
        "ambiguity",
        "in_unique",
        "in_minimal",
        "in_expanded",
    ]
    data = []

    for peptide_pair_id in sorted(data_set.peptide_pairs.keys()):
        peptide_pair = data_set.peptide_pairs[peptide_pair_id]

        # Calculate ambiguity
        expanded_protein_pairs = 0
        for prot_pair_id in peptide_pair.connections:
            prot_pair = data_set.protein_pairs[prot_pair_id]
            if prot_pair.report_status <= ReportStatus.EXPANDED:
                expanded_protein_pairs += 1

        ambiguity = (
            "ambiguous_protein" if expanded_protein_pairs > 1 else "no_ambiguity"
        )

        if len(peptide_pair.linkage_pairs) > 0:
            for linkage_pair in peptide_pair.sorted_linkage_pairs():
                prot_pair = data_set.protein_pairs[linkage_pair.protein_pair_id]
                in_unique = (
                    "True"
                    if prot_pair.report_status <= ReportStatus.UNIQUE
                    else "False"
                )
                in_minimal = (
                    "True"
                    if prot_pair.report_status <= ReportStatus.MINIMAL
                    else "False"
                )
                in_expanded = (
                    "True"
                    if prot_pair.report_status <= ReportStatus.EXPANDED
                    else "False"
                )
                row = [
                    str(linkage_pair),
                    linkage_pair.peptide_a,
                    str(linkage_pair.peptide_linkage_a or ""),
                    str(linkage_pair.start_a or ""),
                    str(linkage_pair.end_a or ""),
                    linkage_pair.peptide_b,
                    str(linkage_pair.peptide_linkage_b or ""),
                    str(linkage_pair.start_b or ""),
                    str(linkage_pair.end_b or ""),
                    linkage_pair.protein_a,
                    str(linkage_pair.protein_linkage_a or ""),
                    linkage_pair.protein_b,
                    str(linkage_pair.protein_linkage_b or ""),
                    ambiguity,
                    in_unique,
                    in_minimal,
                    in_expanded,
                ]
                data.append("\t".join(row))
        else:
            for mapped_prot_a in sorted(peptide_pair.a.mapped_proteins):
                for mapped_prot_b in sorted(peptide_pair.b.mapped_proteins):
                    prot_pair = data_set.protein_pairs[
                        get_pair_id_from_str(mapped_prot_a, mapped_prot_b)
                    ]
                    in_unique = (
                        "True"
                        if prot_pair.report_status <= ReportStatus.UNIQUE
                        else "False"
                    )
                    in_minimal = (
                        "True"
                        if prot_pair.report_status <= ReportStatus.MINIMAL
                        else "False"
                    )
                    in_expanded = (
                        "True"
                        if prot_pair.report_status <= ReportStatus.EXPANDED
                        else "False"
                    )
                    start_a, end_a = peptide_pair.a.protein_locations.get(
                        mapped_prot_a, ("", "")
                    )
                    start_b, end_b = peptide_pair.b.protein_locations.get(
                        mapped_prot_b, ("", "")
                    )

                    pair_id_str = f"{mapped_prot_a}{config.advanced.pair_separator}{mapped_prot_b}_{peptide_pair.a.sequence}{config.advanced.pair_separator}{peptide_pair.b.sequence}"  # noqa: E501

                    row = [
                        pair_id_str,
                        peptide_pair.a.sequence,
                        "",
                        str(start_a) if start_a else "",
                        str(end_a) if end_a else "",
                        peptide_pair.b.sequence,
                        "",
                        str(start_b) if start_b else "",
                        str(end_b) if end_b else "",
                        mapped_prot_a,
                        "",
                        mapped_prot_b,
                        "",
                        ambiguity,
                        in_unique,
                        in_minimal,
                        in_expanded,
                    ]
                    data.append("\t".join(row))

    with open(output_path, "w") as w:
        w.write("\t".join(header) + "\n")
        w.write("\n".join(data) + "\n")


def mapping_peptide_pair_to_protein_pairs(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Write mapping table with peptide pairs and their corresponding protein pairs.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

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


def mapping_peptide_pair_to_protein_pairs_with_linkages(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Write mapping table with peptide pairs and their corresponding protein pairs.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_protein_pair_with_linkages.tsv")
    data = []
    for peptide_pair_id in sorted(data_set.peptide_pairs.keys()):
        peptide_pair = data_set.peptide_pairs[peptide_pair_id]
        row = [str(peptide_pair)]
        for prot_pair_id in sorted(peptide_pair.connections):
            linkage_pairs = peptide_pair.protein_pair_linkage_pairs.get(
                prot_pair_id, set()
            )
            if len(linkage_pairs) == 0:
                row.append(prot_pair_id)
                continue
            for linkage_pair in sorted(
                linkage_pairs, key=lambda linkage_pair: linkage_pair.sort_key()
            ):
                row.append(linkage_pair.protein_pair_with_linkages())
        data.append("\t".join(row))
    with open(output_path, "w") as w:
        w.write("\n".join(data))


def mapping_peptide_to_protein(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Write mapping table that shows a peptide and its corresponding mapped proteins.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_protein.tsv")
    data = []
    uniq_peptides = set()
    for peptide_pair in sorted(
        data_set.peptide_pairs.values(),
        key=lambda p: p.a.sequence + p.b.sequence,
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


def mapping_peptide_to_protein_with_linkages(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Write mapping table that shows a peptide and its corresponding mapped proteins.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_protein_with_linkages.tsv")
    data = []
    uniq_peptides = set()
    for peptide_pair in sorted(
        data_set.peptide_pairs.values(),
        key=lambda p: p.a.sequence + p.b.sequence,
    ):
        for peptide in (peptide_pair.a, peptide_pair.b):
            pep_str = str(peptide)
            if pep_str not in uniq_peptides:
                uniq_peptides.add(pep_str)
                row = [pep_str]
                for prot in peptide.mapped_proteins:
                    linkage = peptide.protein_linkages.get(prot, "")
                    prot_str = f"{prot}:{linkage}" if linkage else prot
                    row.append(prot_str)
                data.append("\t".join(row))
    with open(output_path, "w") as w:
        w.write("\n".join(data))


def mapping_peptide_pair_to_best_protein_pair(
    data_set: XLDataSet, output_path: str | None = None
) -> None:
    """Write mapping table showing the best selected protein pair for each peptide pair.

    Args:
        data_set (XLDataSet): XL data set to create mapping table from
        output_path (str | None, optional): Output path, if None will use default path.
            Defaults to None.

    """
    if output_path is None:
        output_folder = Path(config.output).joinpath("mapping")
        output_folder.mkdir(exist_ok=True)
        output_path = str(output_folder / "peptide_to_best_protein_pair.tsv")
    data = []
    for peptide_pair_id in sorted(data_set.peptide_pairs.keys()):
        peptide_pair = data_set.peptide_pairs[peptide_pair_id]

        best_prot_pair = None
        for prot_pair_id in peptide_pair.connections:
            prot_pair = data_set.protein_pairs[prot_pair_id]
            if best_prot_pair is None:
                best_prot_pair = prot_pair
            else:
                if prot_pair.report_status < best_prot_pair.report_status:
                    best_prot_pair = prot_pair
                elif prot_pair.report_status == best_prot_pair.report_status:
                    if prot_pair.score > best_prot_pair.score:
                        best_prot_pair = prot_pair

        if best_prot_pair is not None:
            data.append(f"{peptide_pair_id}\t{best_prot_pair.pair_id}")
        else:
            data.append(f"{peptide_pair_id}\t")

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
        if len(peptide_pair.linkage_pairs) > 0:
            for linkage_pair in peptide_pair.sorted_linkage_pairs():
                prot_pair = data_set.protein_pairs.get(
                    linkage_pair.protein_pair_id, None
                )
                if prot_pair is None:
                    continue
                rows.append(
                    DetailedRow(
                        peptide_pair=pair_id,
                        protein_pair=prot_pair.pair_id,
                        mapping_order=linkage_pair.mapping_order,
                        peptide_linkage_pair=linkage_pair.peptide_pair_with_linkages(),
                        protein_linkage_pair=linkage_pair.protein_pair_with_linkages(),
                        protein_pair_status=str(prot_pair.report_status),
                        selection_criteria=str(prot_pair.prioritization_status),
                        pair_score=float(prot_pair.score),
                    )
                )
            continue
        for mapped_prot_a in sorted(peptide_pair.a.mapped_proteins):
            for mapped_prot_b in sorted(peptide_pair.b.mapped_proteins):
                prot_pair = data_set.protein_pairs.get(
                    get_pair_id_from_str(mapped_prot_a, mapped_prot_b), None
                )
                if prot_pair is None:
                    continue
                rows.append(
                    DetailedRow(
                        peptide_pair=pair_id,
                        protein_pair=prot_pair.pair_id,
                        mapping_order=f"{mapped_prot_a}{config.advanced.pair_separator}{mapped_prot_b}",
                        protein_pair_status=str(prot_pair.report_status),
                        selection_criteria=str(prot_pair.prioritization_status),
                        pair_score=float(prot_pair.score),
                    )
                )
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
