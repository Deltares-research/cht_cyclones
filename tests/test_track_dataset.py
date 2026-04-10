from pathlib import Path

from cht_cyclones.track_dataset import CycloneTrackDataset


def test_load_track_dataset():
    dataset_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    tc = CycloneTrackDataset("ibtracs", file_name=dataset_file)

    assert tc.nstorms == 13330
