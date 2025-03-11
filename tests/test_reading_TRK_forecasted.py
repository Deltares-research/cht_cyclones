# Load modules
from pathlib import Path

from cht_cyclones.tropical_cyclone import (
    TropicalCyclone,
    TropicalCycloneEnsemble,
)


def test_read_TRK_forecasted(tmp_dir):
    # Define a track
    tc = TropicalCyclone(name="Maria_forecast")

    # Read a track
    file_path = Path(__file__).parent / "TRK_COAMPS_CTCX_3_2017092006_15L"
    tc.read_track(file_path, "trk")

    # Define settings
    # Uses typical Wind Enhance Scheme (WES) formulations
    tc.include_rainfall = (
        True  # using empirical formulation to compute rainfall (IPET default)
    )
    tc.rainfall_factor = (
        2.0  # factor to calibrate rainfall 2.0 means 2x as much as relationship
    )

    # Write as cyc
    tc.convert_units_metric_imperial()
    tc.write_track(tmp_dir / "forecasted_track_Maria.cyc", "ddb_cyc")

    # create (regular) ASCII spiderweb
    tc.to_spiderweb(tmp_dir / "forecasted_track_Maria.spw")

    # create (netcdf) spiderweb
    tc.to_spiderweb(tmp_dir / "forecasted_track_Maria.nc", format_type="netcdf")

    # Write out as shapefile
    tc.track.to_file(tmp_dir / "forecasted_track_Maria.shp")

    # Write out as geojson
    tc.track.to_file(tmp_dir / "forecasted_track_Maria.json", driver="GeoJSON")

    # TC-FF ensemble members
    tc2 = TropicalCycloneEnsemble(name="Maria_forecast_ensembles", TropicalCyclone=tc)
    tc2.dt = 3  # 3 hourly
    tc2.include_best_track = 0  # not included best-track
    tc2.compute_ensemble(100)
    tc2.to_shapefile(tmp_dir)

    suffixes = ["cpg", "cyc", "dbf", "json", "shp", "shx", "prj", "spw"]
    assert all(
        (Path(tmp_dir) / f"forecasted_track_Maria.{suffix}").exists()
        for suffix in suffixes
    )
