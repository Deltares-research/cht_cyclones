import geopandas as gpd
import pytest
from shapely.geometry import Point

from cht_cyclones.fileio import TropicalCycloneTrack
from cht_cyclones.tropical_cyclone import TropicalCyclone


def test_read_ddb_cyc(tmp_dir, track_idai: TropicalCyclone):
    file_path = tmp_dir / "ddb_cyc.cyc"
    track_idai.write_track(file_path, "ddb_cyc")

    track = TropicalCycloneTrack()
    track.read(file_path, "ddb_cyc")

    assert isinstance(track.track, gpd.GeoDataFrame), "track should be a GeoDataFrame"
    assert len(track.track) > 0, "track should contain data"

    # Check if track contains expected columns
    expected_columns = [
        "datetime",
        "geometry",
        "vmax",
        "pc",
        "RMW",
        "R35_NE",
        "R35_SE",
        "R35_SW",
        "R35_NW",
        "R50_NE",
        "R50_SE",
        "R50_SW",
        "R50_NW",
        "R65_NE",
        "R65_SE",
        "R65_SW",
        "R65_NW",
        "R100_NE",
        "R100_SE",
        "R100_SW",
        "R100_NW",
    ]

    for col in expected_columns:
        assert col in track.track.columns, f"{col} should be in the track DataFrame"

    # Validate first row values
    first_row = track.track.iloc[0]

    expected_first_line = {
        "datetime": "20190309 090000",
        "geometry": Point(-17.17, 40.78),
        "vmax": 35.0,
        "pc": 1003.0,
        "RMW": 40.0,
        "R35_NE": 75.0,
        "R35_SE": 75.0,
        "R35_SW": 57.0,
        "R35_NW": 57.0,
        "R50_NE": -999.0,
        "R50_SE": -999.0,
        "R50_SW": -999.0,
        "R50_NW": -999.0,
        "R65_NE": -999.0,
        "R65_SE": -999.0,
        "R65_SW": -999.0,
        "R65_NW": -999.0,
        "R100_NE": -999.0,
        "R100_SE": -999.0,
        "R100_SW": -999.0,
        "R100_NW": -999.0,
    }

    for key, value in expected_first_line.items():
        assert (
            first_row[key] == value
        ), f"{key} does not match expected: {first_row[key]} != {value}"


def test_read_invalid_format():
    track = TropicalCycloneTrack()
    with pytest.raises(
        Exception, match="This file format is not supported as read track!"
    ):
        track.read("non_existent_file.txt", "invalid_format")


def test_empty_file(tmp_dir):
    filename = tmp_dir / "empty.txt"
    with open(filename, "w") as f:
        pass

    track = TropicalCycloneTrack()
    with pytest.raises(Exception):
        track.read(filename, "ddb_cyc")
