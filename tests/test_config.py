"""Test config related functions."""

import xlranker
from xlranker.config import config as xlr_config


def test_loading_config_from_dict():
    """Test that config can be modified from a dict object."""
    # Confirm config is in right state.
    assert not xlr_config.advanced.intra_in_training
    assert xlr_config.output != "new_output"
    config_dict = {"advanced": {"intra_in_training": True}, "output": "new_output"}
    xlranker.config.set_config_from_dict(config_dict)
    assert xlr_config.advanced.intra_in_training
    assert xlr_config.output == "new_output"
