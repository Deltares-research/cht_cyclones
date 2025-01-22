from pathlib import Path

from cht_cyclones.cyclone_track_database import CycloneTrackDatabase


def test_load_track_database():
    database_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    tc = CycloneTrackDatabase("ibtracs", file_name=database_file)

    assert tc.nstorms == 13330
