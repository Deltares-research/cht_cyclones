import tempfile
from pathlib import Path

import pytest


@pytest.fixture()
def tmp_dir():
    """Create and return a temporary directory.

    Upon exiting the context, the directory and everything contained in it are removed.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)
