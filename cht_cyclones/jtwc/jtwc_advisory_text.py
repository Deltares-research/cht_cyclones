import os
import math
from datetime import datetime, timedelta
import re
from geopandas import GeoDataFrame
from shapely.geometry import Point
import pandas as pd
import numpy as np
from pyproj import CRS

def to_gdf(filename):
    """
    Parse a JTWC advisory file and extract cyclone parameters.
    Returns a dictionary with extracted values.
    """
    result = {
        "tc_name": "",
        "date": [],
        "time": [],
        "lat": [],
        "lon": [],
        "mslp": [],
        "rmw": [],
        "msw_kt": [],
        "tc_dir": [],
        "tc_speed": [],
        "radii": [],
        "034NE": [],
        "034SE": [],
        "034SW": [],
        "034NW": [],
        "050NE": [],
        "050SE": [],
        "050SW": [],
        "050NW": [],
        "064NE": [],
        "064SE": [],
        "064SW": [],
        "064NW": [],
        "100NE": [],
        "100SE": [],
        "100SW": [],
        "100NW": [],
    }

    def wind_pressure_relationship(msw):
        # Simplified version of the wind-pressure relationship
        # Table values from Perl code, truncated for brevity
        WPR = [
            (15, 1010.5), (20, 1006.8), (25, 1003.1), (30, 999.5),
            (35, 995.8), (40, 992.1), (45, 988.4), (50, 984.7),
            (55, 981.0), (60, 977.3), (65, 973.6), (70, 969.9),
            (75, 966.2), (80, 962.5), (85, 958.8), (90, 955.0),
            (95, 951.3), (100, 947.6), (105, 943.8), (110, 940.1),
            (115, 936.3), (120, 932.6), (125, 928.8), (130, 925.1),
            (135, 921.3), (140, 917.5), (145, 913.8), (150, 910.0),
            (155, 906.2), (160, 902.4), (165, 898.6), (170, 894.8),
        ]
        for i in range(len(WPR)):
            if WPR[i][0] == msw:
                return WPR[i][1]
            elif i > 0 and WPR[i][0] > msw and WPR[i-1][0] < msw:
                # Linear interpolation
                x0, y0 = WPR[i-1]
                x1, y1 = WPR[i]
                return y0 + (y1-y0)*(msw-x0)/(x1-x0)
        return 1013.25  # Default

    def storm_dir_speed(lon1, lat1, lon2, lat2, tau1, tau2):
        # Calculate storm direction and speed
        if lat1 == lat2 and lon1 == lon2:
            return 0.0, 0.0
        # Haversine formula for distance
        R = 6371.0
        phi1, phi2 = math.radians(lat1), math.radians(lat2)
        dphi = math.radians(lat2 - lat1)
        dlambda = math.radians(lon2 - lon1)
        a = math.sin(dphi/2)**2 + math.cos(phi1)*math.cos(phi2)*math.sin(dlambda/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        km = R * c
        # Bearing
        y = math.sin(dlambda) * math.cos(phi2)
        x = math.cos(phi1)*math.sin(phi2) - math.sin(phi1)*math.cos(phi2)*math.cos(dlambda)
        bearing = (math.degrees(math.atan2(y, x)) + 360) % 360
        # Speed in knots
        nm2km = 1.852
        speed = (km/nm2km) / (tau2 - tau1) if (tau2 - tau1) != 0 else 0.0
        return bearing, speed

    with open(filename, "r") as f:
        lines = f.readlines()

    nrec = -1
    is_analysis = False
    tau = []
    for idx, line in enumerate(lines):
        line = line.strip()
        fields = line.split()
        # Storm name
        if line.startswith("1.") and "WARNING NR" in line:
            result["tc_name"] = " ".join([w for w in fields if not any(c.isdigit() for c in w)])
        # Date/time/position
        if fields and fields[0].endswith("Z") and len(fields) >= 3 and "POSITION" not in line:
            date_str = fields[0]
            lat_str = fields[-2]
            lon_str = fields[-1]
            # Parse date/time
            dt = datetime.strptime(date_str[:6], "%d%H%M")
            date_val = dt.strftime("%Y%m%d")
            time_val = dt.strftime("%H%M")
            result["date"].append(date_val)
            result["time"].append(time_val)
            # Parse lat/lon
            lat = float(lat_str[:-1])
            if lat_str[-1] == "S":
                lat *= -1
            lon = float(lon_str[:-1])
            if lon_str[-1] == "W":
                lon *= -1
            result["lat"].append(lat)
            result["lon"].append(lon)
            nrec += 1
            # Append NaN to radii
            result["034NE"].append(np.nan)
            result["034SE"].append(np.nan)
            result["034SW"].append(np.nan)
            result["034NW"].append(np.nan)
            result["050NE"].append(np.nan)
            result["050SE"].append(np.nan)
            result["050SW"].append(np.nan)
            result["050NW"].append(np.nan)
            result["064NE"].append(np.nan)
            result["064SE"].append(np.nan)
            result["064SW"].append(np.nan)
            result["064NW"].append(np.nan)
            result["100NE"].append(np.nan)
            result["100SE"].append(np.nan)
            result["100SW"].append(np.nan)
            result["100NW"].append(np.nan)
        # Max sustained winds
        if "MAX SUSTAINED WINDS -" in line:
            msw = float(fields[4])
            result["msw_kt"].append(msw)
            result["mslp"].append(wind_pressure_relationship(msw))
        # Movement
        if "MOVEMENT PAST SIX HOURS" in line:
            result["tc_dir"].append(float(fields[5]))
            result["tc_speed"].append(float(fields[8]))
            tau.append(0)
        # Wind radii (simplified, only one set)
        if "RADIUS OF" in line and "QUADRANT" in line:
            # active radius wind speed
            wspradius = fields[2]
        if "NORTHEAST QUADRANT" in line:
            tmpstr = wspradius + "NE"
            result[tmpstr][nrec] = float(fields[-4])
        if "SOUTHEAST QUADRANT" in line:
            tmpstr = wspradius + "SE"
            result[tmpstr][nrec] = float(fields[-4])
        if "SOUTHWEST QUADRANT" in line:
            tmpstr = wspradius + "SW"
            result[tmpstr][nrec] = float(fields[-4])
        if "NORTHWEST QUADRANT" in line:
            tmpstr = wspradius + "NW"
            result[tmpstr][nrec] = float(fields[-4])
        # if "RADIUS OF" in line and "QUADRANT" in line:
        #     # Example: RADIUS OF 034 KT WINDS - 120 NM NORTHEAST QUADRANT
        #     try:
        #         speed = int(fields[2])
        #         radius = int(fields[6])
        #         result["radii"].append((speed, radius))
        #     except Exception:
        #         pass

    # Calculate direction/speed for forecasts
    for i in range(1, len(result["lat"])):
        dir, spd = storm_dir_speed(
            result["lon"][i-1], result["lat"][i-1],
            result["lon"][i], result["lat"][i],
            tau[i-1] if i-1 < len(tau) else 0,
            tau[i] if i < len(tau) else 6
        )
        result["tc_dir"].append(dir)
        result["tc_speed"].append(spd)
 
    gdf = GeoDataFrame()  # Initialize empty GeoDataFrame

    for itime, timestr in enumerate(result["time"]):
        x = result["lon"][itime]
        y = result["lat"][itime]
        vmax = result["msw_kt"][itime] if itime < len(result["msw_kt"]) else np.nan
        pc = result["mslp"][itime] if itime < len(result["mslp"]) else np.nan
        # Wind radii
        r35_ne = result["034NE"][itime]
        r35_se = result["034SE"][itime]
        r35_sw = result["034SW"][itime]
        r35_nw = result["034NW"][itime]
        r50_ne = result["050NE"][itime]
        r50_se = result["050SE"][itime]
        r50_sw = result["050SW"][itime]
        r50_nw = result["050NW"][itime]
        r65_ne = result["064NE"][itime]
        r65_se = result["064SE"][itime]
        r65_sw = result["064SW"][itime]
        r65_nw = result["064NW"][itime]
        r100_ne = result["100NE"][itime]
        r100_se = result["100SE"][itime]
        r100_sw = result["100SW"][itime]
        r100_nw = result["100NW"][itime]
        # r35_ne = result["034NE"][itime] if itime < len(result["034NE"]) else np.nan
        # r35_ne = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r35_se = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r35_sw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r35_nw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r50_ne = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r50_se = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r50_sw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r50_nw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r65_ne = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r65_se = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r65_sw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r65_nw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r100_ne = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r100_se = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r100_sw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        # r100_nw = result["radii"][itime][1] if itime < len(result["radii"]) else np.nan
        datestr = result["date"][itime] if itime < len(result["date"]) else ""
        time = datetime.strptime(datestr + " " + timestr, "%Y%m%d %H%M")
        tc_time_string = time.strftime("%Y%m%d %H%M%S")
        # Create point geometry
        point = Point(x, y)
        # Create GeoDataFrame for this point
        gdf_point = GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [vmax],
                "pc": [pc],
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

    return gdf, result["tc_name"]
