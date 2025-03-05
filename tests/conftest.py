import tempfile
from pathlib import Path

import pytest

from cht_cyclones.cyclone_track_database import CycloneTrackDatabase
from cht_cyclones.tropical_cyclone import TropicalCyclone


@pytest.fixture()
def tmp_dir():
    """Create and return a temporary directory.

    Upon exiting the context, the directory and everything contained in it are removed.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


@pytest.fixture()
def track_idai() -> TropicalCyclone:
    database_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    db = CycloneTrackDatabase("ibtracs", file_name=database_file)
    ind = db.list_names().index("IDAI")
    tc = db.get_track(index=ind)

    return tc
