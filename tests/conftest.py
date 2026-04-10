import tempfile
from pathlib import Path

import pytest

from cht_cyclones.track_dataset import CycloneTrackDataset
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
    dataset_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    # Tests need a directory with metadata.tml; skip if not available
    if not (dataset_file.parent / "metadata.tml").exists():
        pytest.skip("IBTrACS test dataset not set up (missing metadata.tml)")
    db = CycloneTrackDataset("ibtracs", dataset_file.parent)
    ind = db.list_names().index("IDAI")
    tc = db.get_track(index=ind)

    return tc
