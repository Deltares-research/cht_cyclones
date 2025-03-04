# Load modules
from pathlib import Path

from cht_cyclones.cyclone_track_database import CycloneTrackDatabase
from cht_cyclones.tropical_cyclone import TropicalCyclone


def test_reading_writing(tmp_dir):
    # Create track file
    database_file = Path(__file__).parent / "IBTrACS.ALL.v04r00.nc"
    db = CycloneTrackDatabase("ibtracs", file_name=database_file)
    ind = db.list_names().index("IDAI")
    tc = db.get_track(index=ind)
    tc.write_track(tmp_dir / "best_track_idai.cyc", "ddb_cyc")

    # Define a track
    tc = TropicalCyclone(name="Idai")

    # Read a track
    tc.from_ddb_cyc(tmp_dir / "best_track_idai.cyc")

    # Define settings
    # Uses typical Wind Enhance Scheme (WES) formulations
    tc.include_rainfall = (
        True  # using empirical formulation to compute rainfall (IPET default)
    )
    tc.rainfall_factor = (
        2.0  # factor to calibrate rainfall 2.0 means 2x as much as relationship
    )

    # create (regular) ASCII spiderweb
    tc.to_spiderweb(tmp_dir / "best_track_idai.spw")

    # create (netcdf) spiderweb
    tc.to_spiderweb(tmp_dir / "best_track_idai.nc", format_type="netcdf")

    # Write out as shapefile
    tc.track.to_file(tmp_dir / "best_track_idai.shp")

    # Write out as geojson
    tc.track.to_file(tmp_dir / "best_track_idai.json", driver="GeoJSON")

    # Write as cyc again
    tc.convert_units_metric_imperial()
    tc.write_track(tmp_dir / "best_track_idai_cht.cyc", "ddb_cyc")

    suffixes = ["cpg", "cyc", "dbf", "json", "nc", "shp", "shx", "prj", "spw", "nc"]
    assert all(
        (Path(tmp_dir) / f"best_track_idai.{suffix}").exists() for suffix in suffixes
    )
