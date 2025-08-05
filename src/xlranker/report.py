"""Report helper functions."""

from pathlib import Path

from xlranker.bio.pairs import ProteinPair
from xlranker.config import config
from xlranker.lib import write_pair_to_network
from xlranker.status import ReportStatus


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
