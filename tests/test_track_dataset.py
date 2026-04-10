from pathlib import Path

import pytest

from cht_cyclones.track_dataset import CycloneTrackDataset


def test_load_track_dataset():
    dataset_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    if not (dataset_file.parent / "metadata.tml").exists():
        pytest.skip("IBTrACS test dataset not set up (missing metadata.tml)")
    tc = CycloneTrackDataset("ibtracs", dataset_file.parent)

    assert tc.nstorms == 13330
