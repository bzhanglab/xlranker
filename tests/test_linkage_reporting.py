"""Tests for linkage provenance and reporting."""

import xlranker as xlr


def test_peptide_instances_do_not_share_linkage_state():
    """Protein linkage mappings should not leak between peptide instances."""
    peptide_a = xlr.bio.Peptide("PEPA")
    peptide_b = xlr.bio.Peptide("PEPB")

    peptide_a.protein_linkages["P1"] = "K10"
    peptide_a.mapped_proteins.append("P1")

    assert peptide_b.protein_linkages == {}
    assert peptide_b.mapped_proteins == []


def test_build_proteins_tracks_individual_linkage_pairs():
    """Each mapped protein linkage should be tracked as its own linkage pair."""
    peptide_pair = xlr.lib.PeptidePair(
        xlr.bio.Peptide(
            "PEPA",
            linkage=2,
            mapped_proteins=["ESR1"],
            protein_linkages={"ESR1": "K10"},
        ),
        xlr.bio.Peptide(
            "PEPB",
            linkage=3,
            mapped_proteins=["ESR2", "ESR2B"],
            protein_linkages={"ESR2": "K20", "ESR2B": "K25"},
        ),
    )
    data_set = xlr.XLDataSet({peptide_pair.pair_id: peptide_pair}, {})

    data_set.build_proteins()

    assert peptide_pair.connections == {"ESR1+ESR2", "ESR1+ESR2B"}
    assert len(peptide_pair.linkage_pairs) == 2

    linkage_strings = {
        linkage_pair.protein_pair_with_linkages()
        for linkage_pair in peptide_pair.linkage_pairs
    }
    assert linkage_strings == {"ESR1:K10+ESR2:K20", "ESR1:K10+ESR2B:K25"}

    for protein_pair_id in peptide_pair.connections:
        assert len(data_set.protein_pairs[protein_pair_id].linkage_pairs) == 1


def test_network_with_linkages_reports_one_row_per_linkage(tmp_path):
    """Protein linkage reports should expand one protein pair into one row per linkage."""  # noqa: E501
    protein_a = xlr.bio.Protein("ESR1", "ESR1")
    protein_b = xlr.bio.Protein("ESR2", "ESR2")
    protein_pair = xlr.lib.ProteinPair(protein_a, protein_b)
    protein_pair.add_linkage_pair(
        xlr.lib.LinkagePair(
            peptide_pair_id="PEPA:2+PEPB:3",
            protein_pair_id=protein_pair.pair_id,
            mapping_order="ESR1+ESR2",
            peptide_a="PEPA",
            peptide_b="PEPB",
            peptide_linkage_a=2,
            peptide_linkage_b=3,
            protein_a="ESR1",
            protein_b="ESR2",
            protein_linkage_a="K10",
            protein_linkage_b="K20",
        )
    )
    protein_pair.add_linkage_pair(
        xlr.lib.LinkagePair(
            peptide_pair_id="PEPC:4+PEPD:5",
            protein_pair_id=protein_pair.pair_id,
            mapping_order="ESR1+ESR2",
            peptide_a="PEPC",
            peptide_b="PEPD",
            peptide_linkage_a=4,
            peptide_linkage_b=5,
            protein_a="ESR1",
            protein_b="ESR2",
            protein_linkage_a="K15",
            protein_linkage_b="K25",
        )
    )

    output_path = tmp_path / "network.tsv"
    xlr.lib.write_pair_to_network_with_linkages([protein_pair], str(output_path))

    assert output_path.read_text().splitlines() == [
        "ESR1\tK10\tESR2\tK20",
        "ESR1\tK15\tESR2\tK25",
    ]


def test_mapping_table_with_linkages_expands_linkage_columns(tmp_path):
    """Peptide-pair mapping table should preserve each linkage pair separately."""
    peptide_pair = xlr.lib.PeptidePair(
        xlr.bio.Peptide(
            "PEPA",
            linkage=2,
            mapped_proteins=["ESR1"],
            protein_linkages={"ESR1": "K10"},
        ),
        xlr.bio.Peptide(
            "PEPB",
            linkage=3,
            mapped_proteins=["ESR2", "ESR2B"],
            protein_linkages={"ESR2": "K20", "ESR2B": "K25"},
        ),
    )
    data_set = xlr.XLDataSet({peptide_pair.pair_id: peptide_pair}, {})

    data_set.build_proteins()

    output_path = tmp_path / "mapping.tsv"
    xlr.report.mapping_peptide_pair_to_protein_pairs_with_linkages(
        data_set, str(output_path)
    )

    assert output_path.read_text().strip() == (
        "PEPA:2+PEPB:3\tESR1:K10+ESR2:K20\tESR1:K10+ESR2B:K25"
    )


def test_make_all_reports_skips_linkage_reports_when_not_provided(tmp_path, monkeypatch):
    """Standard reports should be the only report files when input has no linkages."""
    monkeypatch.setattr(xlr.config.config, "output", str(tmp_path))

    peptide_pair = xlr.lib.PeptidePair(
        xlr.bio.Peptide("PEPA", mapped_proteins=["ESR1"]),
        xlr.bio.Peptide("PEPB", mapped_proteins=["ESR2"]),
    )
    data_set = xlr.XLDataSet(
        {peptide_pair.pair_id: peptide_pair}, {}, has_linkages=False
    )
    data_set.build_proteins()

    xlr.report.make_all_reports(data_set)

    report_dir = tmp_path / "reports"
    assert (report_dir / "unique.tsv").exists()
    assert (report_dir / "minimal.tsv").exists()
    assert (report_dir / "expanded.tsv").exists()
    assert (report_dir / "all.tsv").exists()
    assert not (report_dir / "unique_with_linkages.tsv").exists()
    assert not (report_dir / "minimal_with_linkages.tsv").exists()
    assert not (report_dir / "expanded_with_linkages.tsv").exists()
    assert not (report_dir / "all_with_linkages.tsv").exists()
