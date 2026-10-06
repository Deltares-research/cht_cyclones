"""
Reader for NHC tropical-cyclone Forecast/Advisory (TCM) products.

The Forecast/Advisory (product identifier ``TCMxxN``, also known as the marine
advisory) is the National Hurricane Center's machine-readable bulletin.  It
contains the current (analysis) storm position, intensity and quadrant wind
radii, followed by a series of forecast positions at 12/24/36/48/72/96/120 h.

This module parses that plain-text product into a GeoDataFrame with the same
columns used by the rest of ``cht_cyclones`` (see
:func:`cht_cyclones.jtwc.jmv30.to_gdf`), so an NHC forecast track can be handled
identically to a JTWC JMV 3.0 track.
"""

import re
from datetime import datetime

import numpy as np
import pandas as pd
from geopandas import GeoDataFrame
from pyproj import CRS
from shapely.geometry import Point

# Storm-type phrases that may precede the storm name in the header line.
_STORM_TYPES = (
    "POTENTIAL TROPICAL CYCLONE",
    "SUBTROPICAL DEPRESSION",
    "SUBTROPICAL STORM",
    "TROPICAL DEPRESSION",
    "TROPICAL STORM",
    "POST-TROPICAL CYCLONE",
    "REMNANTS OF",
    "HURRICANE",
    "TYPHOON",
)

# NHC wind-radii thresholds (kt) mapped to the cht_cyclones column prefixes.
_RADII_MAP = {34: "r35", 50: "r50", 64: "r65"}

# Regular expressions for the header / body of the product.
_RE_NAME = re.compile(
    r"(" + "|".join(_STORM_TYPES) + r")\s+(.*?)\s+"
    r"FORECAST/ADVISORY\s+NUMBER\s+(\S+)",
    re.IGNORECASE,
)
_RE_STORM_ID = re.compile(r"\b([A-Z]{2}\d{6})\b")
_RE_ISSUANCE = re.compile(
    r"^\s*(\d{4})\s+UTC\s+\w{3}\s+(\w{3})\s+(\d{1,2})\s+(\d{4})",
    re.IGNORECASE,
)
_RE_CENTER = re.compile(
    r"CENTER\s+LOCATED\s+NEAR\s+([\d.]+)([NS])\s+([\d.]+)([EW])\s+AT\s+(\d{2})/(\d{4})Z",
    re.IGNORECASE,
)
_RE_VALID = re.compile(
    r"(?:FORECAST|OUTLOOK)\s+VALID\s+(\d{2})/(\d{4})Z\s+([\d.]+)([NS])\s+([\d.]+)([EW])",
    re.IGNORECASE,
)
_RE_CP = re.compile(r"CENTRAL\s+PRESSURE\s+(\d+)\s*MB", re.IGNORECASE)
_RE_MAX_SUSTAINED = re.compile(r"MAX\s+SUSTAINED\s+WINDS\s+(\d+)\s*KT", re.IGNORECASE)
_RE_MAX_WIND = re.compile(r"MAX\s+WIND\s+(\d+)\s*KT", re.IGNORECASE)
_RE_RADII_LINE = re.compile(r"^\s*(\d+)\s*KT\b", re.IGNORECASE)
_RE_QUADRANT = re.compile(r"(\d+)\s*(NE|SE|SW|NW)", re.IGNORECASE)

_MONTHS = {
    "JAN": 1,
    "FEB": 2,
    "MAR": 3,
    "APR": 4,
    "MAY": 5,
    "JUN": 6,
    "JUL": 7,
    "AUG": 8,
    "SEP": 9,
    "OCT": 10,
    "NOV": 11,
    "DEC": 12,
}


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
        for quad in ("ne", "se", "sw", "nw"):
            rec[f"{r}_{quad}"] = np.nan
    return rec


def _reconstruct_time(day: int, hhmm: str, issue: datetime) -> datetime:
    """
    Build a full datetime from a ``DD/HHMM`` forecast stamp and the issuance time.

    The Forecast/Advisory only gives the day-of-month and hour for each forecast
    valid time; the year and month are taken from the advisory issuance time,
    accounting for a possible month/year rollover when the forecast day is
    smaller than the issuance day.
    """
    hour = int(hhmm[:2])
    minute = int(hhmm[2:])
    year = issue.year
    month = issue.month
    if day < issue.day:
        # Forecast valid time has rolled into the next month.
        month += 1
        if month > 12:
            month = 1
            year += 1
    return datetime(year, month, day, hour, minute)


def _parse_radii(line: str, rec: dict) -> None:
    """Parse a ``NN KT... <r>NE <r>SE <r>SW <r>NW`` line into a record."""
    m = _RE_RADII_LINE.match(line)
    if not m or "SEAS" in line.upper():
        # Not a wind-radii line (e.g. "12 FT SEAS..." or a MAX WIND line).
        return
    threshold = int(m.group(1))
    prefix = _RADII_MAP.get(threshold)
    if prefix is None:
        return
    for value, quad in _RE_QUADRANT.findall(line):
        v = float(value)
        # NHC reports 0 for a quadrant with no winds of that strength; keep as NaN.
        rec[f"{prefix}_{quad.lower()}"] = v if v > 0 else np.nan


def to_dict(fname: str) -> dict:
    """
    Parse an NHC Forecast/Advisory (TCM) text file into a dictionary.

    Parameters
    ----------
    fname : str
        Path to the plain-text Forecast/Advisory product.

    Returns
    -------
    dict
        Dictionary with keys ``"name"``, ``"advisorynumber"``, ``"storm_id"``,
        and ``"records"`` (a list of per-time-step record dictionaries in
        chronological order).
    """
    with open(fname, "r") as f:
        lines = f.readlines()

    name = ""
    advisory = None
    storm_id = None
    issue = None
    records = []
    current = None
    initial_done = False
    forecast_started = False

    for raw in lines:
        line = raw.rstrip("\n").rstrip()

        # --- Header: storm type / name / advisory number ---
        if name == "":
            m = _RE_NAME.search(line)
            if m:
                name = m.group(2).strip()
                advnum = re.sub(r"\D", "", m.group(3))
                advisory = int(advnum) if advnum else None

        # --- Header: storm id (e.g. AL092024, EP042026) ---
        if storm_id is None and "HURRICANE CENTER" in line.upper():
            m = _RE_STORM_ID.search(line)
            if m:
                storm_id = m.group(1).upper()

        # --- Header: issuance date/time ---
        if issue is None:
            m = _RE_ISSUANCE.match(line)
            if m:
                hhmm, mon, day, year = m.groups()
                month = _MONTHS.get(mon.upper())
                if month is not None:
                    issue = datetime(
                        int(year), month, int(day), int(hhmm[:2]), int(hhmm[2:])
                    )

        # --- Analysis (initial) position ---
        if not initial_done and not forecast_started:
            m = _RE_CENTER.search(line)
            if m and issue is not None:
                lat = float(m.group(1)) * (-1 if m.group(2).upper() == "S" else 1)
                lon = float(m.group(3)) * (-1 if m.group(4).upper() == "W" else 1)
                day = int(m.group(5))
                t = _reconstruct_time(day, m.group(6), issue)
                current = _blank_record(t, lon, lat)
                records.append(current)
                initial_done = True
                continue

        # --- Analysis intensity / pressure ---
        if current is not None and not forecast_started:
            mcp = _RE_CP.search(line)
            if mcp:
                current["pc"] = float(mcp.group(1))
            msw = _RE_MAX_SUSTAINED.search(line)
            if msw:
                current["vmax"] = float(msw.group(1))

        # --- Forecast / outlook position ---
        m = _RE_VALID.search(line)
        if m and issue is not None:
            forecast_started = True
            lat = float(m.group(3)) * (-1 if m.group(4).upper() == "S" else 1)
            lon = float(m.group(5)) * (-1 if m.group(6).upper() == "W" else 1)
            day = int(m.group(1))
            t = _reconstruct_time(day, m.group(2), issue)
            current = _blank_record(t, lon, lat)
            records.append(current)
            continue

        # --- Forecast intensity ---
        if current is not None and forecast_started:
            mw = _RE_MAX_WIND.search(line)
            if mw:
                current["vmax"] = float(mw.group(1))

        # --- Wind radii (both analysis and forecast blocks) ---
        if current is not None and _RE_RADII_LINE.match(line):
            _parse_radii(line, current)

    return {
        "name": name,
        "advisorynumber": advisory,
        "storm_id": storm_id,
        "records": records,
    }


def to_gdf(fname: str) -> tuple:
    """
    Read an NHC Forecast/Advisory file and return a GeoDataFrame of the track.

    Parameters
    ----------
    fname : str
        Path to the Forecast/Advisory (TCM) text file.

    Returns
    -------
    gdf : geopandas.GeoDataFrame
        One row per time step; columns include ``datetime``, ``geometry``,
        ``vmax``, ``pc``, ``rmw``, and quadrant wind radii, matching the layout
        produced by :func:`cht_cyclones.jtwc.jmv30.to_gdf`.
    name : str
        Storm name extracted from the advisory.
    advisory : int or None
        Advisory number extracted from the advisory.
    """
    tc = to_dict(fname)

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

    return gdf, tc["name"], tc["advisorynumber"]
