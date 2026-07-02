"""
Reader for NHC ATCF forecast files (``.fst``).

The Automated Tropical Cyclone Forecasting (ATCF) system stores tropical-cyclone
forecasts as comma-delimited tables on the NHC server
(``https://ftp.nhc.noaa.gov/atcf/``).  Each ``fst/<basin><cy><year>.fst`` file
holds the National Hurricane Center's **official (OFCL) forecast**: one row per
forecast lead time (``TAU``) and per wind-radii threshold (34/50/64 kt), with the
quadrant wind radii given as NE/SE/SW/NW.

This module parses that table into a GeoDataFrame with the same columns used by
the rest of ``cht_cyclones`` (see :func:`cht_cyclones.jtwc.jmv30.to_gdf`), so an
NHC ATCF forecast track can be handled identically to a JTWC JMV 3.0 track.

ATCF column layout (0-based index) used here::

    0  BASIN        6  Lat (tenths, N/S)   12 WINDCODE      19 RMW
    1  CY           7  Lon (tenths, E/W)   13 RAD1 (NE)     20 GUSTS
    2  YYYYMMDDHH   8  VMAX (kt)           14 RAD2 (SE)     ...
    3  TECHNUM      9  MSLP (mb)           15 RAD3 (SW)     27 STORMNAME
    4  TECH         10 TY                  16 RAD4 (NW)
    5  TAU (hours)  11 RAD (34/50/64)
"""

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
from pyproj import CRS
from shapely.geometry import Point

# ATCF wind-radii thresholds (kt) mapped to the cht_cyclones column prefixes.
_RADII_MAP = {34: "r35", 50: "r50", 64: "r65", 100: "r100"}

# Quadrant order for the standard NEQ wind code.
_QUADRANTS = ("ne", "se", "sw", "nw")


def _parse_latlon(token: str) -> float:
    """Convert an ATCF lat/lon token (e.g. ``294N`` or ``1268W``) to degrees."""
    hemi = token[-1].upper()
    value = float(token[:-1]) * 0.1
    if hemi in ("S", "W"):
        value = -value
    return value


def _blank_record(time: datetime, lon: float, lat: float) -> dict:
    """Return a track record dict pre-filled with NaNs for a given point."""
    rec = {
        "time": time,
        "x": lon,
        "y": lat,
        "vmax": np.nan,
        "pc": np.nan,
        "rmw": np.nan,
    }
    for r in ("r35", "r50", "r65", "r100"):
        for quad in _QUADRANTS:
            rec[f"{r}_{quad}"] = np.nan
    return rec


def to_dict(fname: str, tech: str = "OFCL") -> dict:
    """
    Parse an ATCF ``.fst`` file into a dictionary.

    Rows are grouped by forecast lead time (``TAU``); the 34/50/64 kt wind-radii
    rows for a given lead time are merged into a single record.

    Parameters
    ----------
    fname : str
        Path to the ATCF forecast file.
    tech : str, optional
        Forecast technique/aid to extract (default ``"OFCL"``, the official
        forecast).  Rows for other techniques are ignored; this makes the reader
        usable on multi-aid a-deck files as well as single-aid ``.fst`` files.

    Returns
    -------
    dict
        Dictionary with keys ``"name"``, ``"basin"``, ``"cy"``, ``"init_time"``
        (a :class:`datetime`), and ``"records"`` (a list of per-lead-time record
        dictionaries in chronological order).
    """
    with open(fname, "r") as f:
        lines = f.readlines()

    name = None
    basin = None
    cy = None
    init_time = None
    records = {}  # keyed by tau to merge multiple radii rows

    for line in lines:
        if not line.strip():
            continue
        s = [field.strip() for field in line.split(",")]
        if len(s) < 20:
            continue
        if s[4].upper() != tech.upper():
            continue

        basin = s[0]
        cy = s[1]
        init_time = datetime.strptime(s[2], "%Y%m%d%H")
        tau = int(s[5])
        time = init_time + timedelta(hours=tau)

        # Create the record for this lead time on first encounter.
        if tau not in records:
            lat = _parse_latlon(s[6])
            lon = _parse_latlon(s[7])
            rec = _blank_record(time, lon, lat)
            vmax = float(s[8])
            rec["vmax"] = vmax if vmax > 0 else np.nan
            pc = float(s[9])
            rec["pc"] = pc if pc > 0 else np.nan
            rmw = float(s[19])
            rec["rmw"] = rmw if rmw > 0 else np.nan
            records[tau] = rec

        rec = records[tau]

        # Merge the wind-radii threshold carried by this row.
        threshold = int(s[11]) if s[11] else 0
        prefix = _RADII_MAP.get(threshold)
        if prefix is not None:
            windcode = s[12].upper()
            radii = [float(s[13]), float(s[14]), float(s[15]), float(s[16])]
            if windcode == "AAA":
                # Full-circle radius: apply to all quadrants.
                radii = [radii[0]] * 4
            for quad, value in zip(_QUADRANTS, radii):
                if value > 0:
                    rec[f"{prefix}_{quad}"] = value

        # Storm name, if present (a-deck field; usually absent in .fst files).
        if len(s) >= 28 and s[27] and s[27].isalpha():
            name = s[27]

    ordered = [records[tau] for tau in sorted(records)]

    return {
        "name": name,
        "basin": basin,
        "cy": cy,
        "init_time": init_time,
        "records": ordered,
    }


def to_gdf(fname: str, tech: str = "OFCL") -> tuple:
    """
    Read an ATCF ``.fst`` file and return a GeoDataFrame of the track.

    Parameters
    ----------
    fname : str
        Path to the ATCF forecast file.
    tech : str, optional
        Forecast technique/aid to extract (default ``"OFCL"``).

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        One row per lead time; columns include ``datetime``, ``geometry``,
        ``vmax``, ``pc``, ``rmw``, and quadrant wind radii, matching the layout
        produced by :func:`cht_cyclones.jtwc.jmv30.to_gdf`.
    name : str or None
        Storm name if present in the file, else ``None``.
    advisory : None
        Always ``None`` (ATCF forecast files carry no advisory number).
    """
    tc = to_dict(fname, tech=tech)

    gdf = GeoDataFrame()  # Initialize empty GeoDataFrame

    for rec in tc["records"]:
        tc_time_string = rec["time"].strftime("%Y%m%d %H%M%S")
        point = Point(rec["x"], rec["y"])
        gdf_point = GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [rec["vmax"]],
                "pc": [rec["pc"]],
                "rmw": [rec["rmw"]],
                "r35_ne": [rec["r35_ne"]],
                "r35_se": [rec["r35_se"]],
                "r35_sw": [rec["r35_sw"]],
                "r35_nw": [rec["r35_nw"]],
                "r50_ne": [rec["r50_ne"]],
                "r50_se": [rec["r50_se"]],
                "r50_sw": [rec["r50_sw"]],
                "r50_nw": [rec["r50_nw"]],
                "r65_ne": [rec["r65_ne"]],
                "r65_se": [rec["r65_se"]],
                "r65_sw": [rec["r65_sw"]],
                "r65_nw": [rec["r65_nw"]],
                "r100_ne": [rec["r100_ne"]],
                "r100_se": [rec["r100_se"]],
                "r100_sw": [rec["r100_sw"]],
                "r100_nw": [rec["r100_nw"]],
            }
        )
        gdf = pd.concat([gdf, gdf_point])

    # Replace -999.0 with NaN (consistency with the other readers)
    gdf = gdf.replace(-999.0, np.nan)
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(crs=CRS(4326), inplace=True)

    return gdf, tc["name"], None
