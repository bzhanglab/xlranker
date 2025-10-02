"""Tests for homolog loading functionality."""

import xlranker.data
from xlranker.config import config


def test_mouse_homologs() -> None:
    """Test loading mouse homologs."""
    config.species = "mmusculus"
    homologs = xlranker.data.load_homologs()
    assert isinstance(homologs, dict)
    assert "Trp53" in homologs
    assert homologs["Trp53"] == "TP53"
    assert "NoHomologGene" not in homologs
