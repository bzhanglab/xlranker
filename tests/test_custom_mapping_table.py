import logging

import xlranker

LOGGER = logging.getLogger(__name__)

DUPLICATE_PEPTIDE_SEQ = """SEQ1\tPROT1
SEQ2\tPROT2\tPROT3
SEQ2\tPROT4
"""


def test_custom_table(tmp_path):
    """tests mapping results with duplicated entries. first entry should be prioritized"""
    temp_file = tmp_path / "table.tsv"
    temp_file.write_text(DUPLICATE_PEPTIDE_SEQ)
    res = xlranker.util.readers.read_mapping_table_file(temp_file)
    assert res["SEQ1"] == ["PROT1"]
    assert res["SEQ2"] == ["PROT2", "PROT3"]


def test_duplicate_entry(tmp_path, caplog):
    """checks to see if log warns user of the duplicated sequence"""
    temp_file = tmp_path / "table.tsv"
    temp_file.write_text(DUPLICATE_PEPTIDE_SEQ)
    with caplog.at_level(logging.WARNING):
        xlranker.util.readers.read_mapping_table_file(temp_file)
    assert any(
        record.levelname == "WARNING"
        and "Peptide sequence SEQ2 duplicated!" in record.message
        for record in caplog.records
    )
