import numpy as np
from datetime import datetime, timedelta
from geopandas import GeoDataFrame
import pandas as pd
from shapely.geometry import Point
from pyproj import CRS

def datenum_to_datetime(dn):
    """Convert MATLAB datenum to Python datetime"""
    return datetime.fromordinal(int(dn)) + timedelta(days=dn%1) - timedelta(days=366)

def datetime_to_datenum(dt):
    """Convert Python datetime to MATLAB datenum"""
    return dt.toordinal() + 366 + (dt - datetime(dt.year, dt.month, dt.day)).total_seconds() / 86400

def jmv30_to_dict(fname):
    tc = {
        "name": "",
        "advisorynumber": None,
        "time": [],
        "x": [], "y": [], "vmax": [],
        "r35ne": [], "r35se": [], "r35sw": [], "r35nw": [],
        "r50ne": [], "r50se": [], "r50sw": [], "r50nw": [],
        "r65ne": [], "r65se": [], "r65sw": [], "r65nw": [],
        "r100ne": [], "r100se": [], "r100sw": [], "r100nw": [],
    }

    # --- Read file ---
    with open(fname, "r") as f:
        lines = f.readlines()

    # --- Read header ---
    str_line = lines[2].strip()
    if not str_line:
        str_line = lines[3].strip()
    parts = str_line.split()
    t0 = datetime.strptime(parts[0][:10], "%Y%m%d%H")
    tc["name"] = parts[2]
    tc["advisorynumber"] = int(parts[3])

    # --- Read forecast section ---
    i = 3
    while i < len(lines):
        line = lines[i].strip()
        parts = line.split()
        if not parts or not parts[0].startswith("T"):
            break
        tc["time"].append(t0 + timedelta(hours=float(parts[0][1:])))
        lat = float(parts[1][:-1]) * 0.1
        if parts[1].endswith("S"): lat = -lat
        lon = float(parts[2][:-1]) * 0.1
        if parts[2].endswith("W"): lon = -lon
        tc["y"].append(lat)
        tc["x"].append(lon)
        tc["vmax"].append(float(parts[3]))

        # Initialize radii
        radii = {"r35": [np.nan]*4, "r50": [np.nan]*4, "r65": [np.nan]*4, "r100": [np.nan]*4}

        # Parse radii
        for n in range(4, len(parts)):
            tag = parts[n].lower()
            if tag in ["r034", "r050", "r064", "r100"]:
                key = "r35" if tag == "r034" else "r50" if tag == "r050" else "r65" if tag == "r064" else "r100"
                radii[key] = [float(parts[n+1]), float(parts[n+4]), float(parts[n+7]), float(parts[n+10])]

        # Append radii
        for quad, idx in zip(["ne","se","sw","nw"], range(4)):
            tc[f"r35{quad}"].append(radii["r35"][idx])
            tc[f"r50{quad}"].append(radii["r50"][idx])
            tc[f"r65{quad}"].append(radii["r65"][idx])
            tc[f"r100{quad}"].append(radii["r100"][idx])

        i += 1

    tc["first_forecast_time"] = tc["time"][0] if tc["time"] else None

    # --- Read hindcast section ---
    istart = istop = None
    for n, line in enumerate(lines):
        parts = line.split()
        if parts and parts[-1].endswith("//"):
            istart = n + 1
        if parts and parts[-1] == "NNNN":
            istop = n - 2
            break

    if istart is None:
        # Did not find istart
        # istart should be the first line that looks something like: 2825100200 168N1508E  15
        for n, line in enumerate(lines):
            parts = line.split()
            if parts and parts[0].isdigit() and (parts[1].endswith("E") or parts[1].endswith("W")) and parts[2].isdigit():
                istart = n
                break

    tc0 = {key: [] for key in tc}
    tc0["name"] = tc["name"]
    tc0["advisorynumber"] = tc["advisorynumber"]

    for line in lines[istart:istop]:
        if not line.strip():
            continue
        f1, f2, f3 = line[:10].strip(), line[11:20].strip(), line[21:].strip()
        time = datetime.strptime(f1[2:], "%y%m%d%H")
        tc0["time"].append(time)

        # Latitude
        if "N" in f2:
            idx = f2.index("N") + 1
            lat = 0.1 * float(f2[:idx-1])
        else:
            idx = f2.index("S") + 1
            lat = -0.1 * float(f2[:idx-1])
        tc0["y"].append(lat)

        # Longitude
        lon_str = f2[idx:]
        lon = 0.1 * float(lon_str[:-1])
        if lon_str.endswith("W"):
            lon = -lon
        tc0["x"].append(lon)

        tc0["vmax"].append(float(f3))

        for quad in ["ne","se","sw","nw"]:
            for r in ["r35","r50","r65","r100"]:
                tc0[f"{r}{quad}"].append(np.nan)

    # --- Merge hindcast & forecast ---
    if abs((tc0["time"][-1] - tc["time"][0]).total_seconds()) < 1:
        for key in tc0:
            if isinstance(tc0[key], list):
                tc0[key] = tc0[key][:-1]

    for key in tc:
        if isinstance(tc[key], list):
            tc[key] = tc0[key] + tc[key]

    # --- Replace NaNs with -999 ---
    for key in tc:
        if isinstance(tc[key], list) and key.startswith("r"):
            tc[key] = [-999 if np.isnan(v) else v for v in tc[key]]

    return tc

def to_gdf(fname):

    tc_dict = jmv30_to_dict(fname)

    gdf = GeoDataFrame()  # Initialize empty GeoDataFrame

    for itime, time in enumerate(tc_dict["time"]):
        x = tc_dict["x"][itime]
        y = tc_dict["y"][itime]
        vmax = tc_dict["vmax"][itime]
        r35_ne = tc_dict["r35ne"][itime]
        r35_se = tc_dict["r35se"][itime]
        r35_sw = tc_dict["r35sw"][itime]
        r35_nw = tc_dict["r35nw"][itime]
        r50_ne = tc_dict["r50ne"][itime]
        r50_se = tc_dict["r50se"][itime]
        r50_sw = tc_dict["r50sw"][itime]
        r50_nw = tc_dict["r50nw"][itime]
        r65_ne = tc_dict["r65ne"][itime]
        r65_se = tc_dict["r65se"][itime]
        r65_sw = tc_dict["r65sw"][itime]
        r65_nw = tc_dict["r65nw"][itime]
        r100_ne = tc_dict["r100ne"][itime]
        r100_se = tc_dict["r100se"][itime]
        r100_sw = tc_dict["r100sw"][itime]
        r100_nw = tc_dict["r100nw"][itime]
        tc_time_string = time.strftime("%Y%m%d %H%M%S")
        # Create point geometry
        point = Point(x, y)
        # Create GeoDataFrame for this point
        gdf_point = GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [vmax],
                "pc": [np.nan],  # Placeholder for pressure
                "rmw": [np.nan],  # Placeholder for RMW
                "r35_ne": [r35_ne],
                "r35_se": [r35_se],
                "r35_sw": [r35_sw],
                "r35_nw": [r35_nw],
                "r50_ne": [r50_ne],
                "r50_se": [r50_se],
                "r50_sw": [r50_sw],
                "r50_nw": [r50_nw],
                "r65_ne": [r65_ne],
                "r65_se": [r65_se],
                "r65_sw": [r65_sw],
                "r65_nw": [r65_nw],
                "r100_ne": [r100_ne],
                "r100_se": [r100_se],
                "r100_sw": [r100_sw],
                "r100_nw": [r100_nw],
            }
        )

        # Append self
        gdf = pd.concat([gdf, gdf_point])

    # Replace -999.0 with NaN
    gdf = gdf.replace(-999.0, np.nan)
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(crs=CRS(4326), inplace=True)

    name = tc_dict["name"]
    advisory = tc_dict["advisorynumber"]

    return gdf, name, advisory
