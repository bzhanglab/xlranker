"""Test config related functions."""

import xlranker
from xlranker.config import config as xlr_config
import xlranker.config
import pytest
from pydantic import ValidationError


EXAMPLE_JSON = """
{
  "detailed": true,
  "advanced": {
    "save_model_files": true
  },
  "mapping": {
    "fasta_type": "GENCODE"
  }
}
"""


def test_loading_config_from_dict():
    """Test that config can be modified from a dict object."""
    # Confirm config is in right state.
    assert not xlr_config.advanced.intra_in_training
    assert xlr_config.output != "new_output"
    config_dict = {"advanced": {"intra_in_training": True}, "output": "new_output"}
    xlranker.config.set_config_from_dict(config_dict)
    assert xlr_config.advanced.intra_in_training
    assert xlr_config.output == "new_output"


def test_json_loading(tmp_path):
    """Test if json file loading sets correct config options."""
    temp_file = tmp_path / "xlranker_example.json"
    temp_file.write_text(EXAMPLE_JSON)
    global xlr_config
    xlr_config = xlranker.config.Config()
    assert not xlr_config.detailed
    assert not xlr_config.advanced.save_model_files
    assert xlr_config.mapping.fasta_type != "GENCODE"

    xlranker.config.load_from_json(temp_file)
    assert xlr_config.detailed
    assert xlr_config.advanced.save_model_files
    assert xlr_config.mapping.fasta_type == "GENCODE"


def test_bad_typing():
    """Test to verify that setting incorrect value type triggers ValidationError.

    Test:
        config.advanced.intra_in_training is a bool, but provide string
    """
    xlr_config = xlranker.config.Config()  # Reset back to default
    assert not xlr_config.advanced.intra_in_training
    assert xlr_config.output != "new_output"
    config_dict = {
        "advanced": {"intra_in_training": "BAD_TYPE"},
        "output": "new_output",
    }
    with pytest.raises(ValidationError):
        xlranker.config.set_config_from_dict(config_dict)
