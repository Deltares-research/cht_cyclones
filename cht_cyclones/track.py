import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from shapely.geometry import LineString, MultiLineString, Point, mapping
from geojson import Feature, FeatureCollection
from geopandas import GeoDataFrame
from scipy.interpolate import CubicSpline, interp1d
import os

from .wind_profiles import wpr_holland2008, wind_radii_nederhoff
from .utils import gdf_to_geojson_js

knots_to_ms = float(0.51444)
nm_to_km = float(1.852)
nm_to_m = float(1.852) * 1000
dateformat_module = "%Y%m%d %H%M%S"

class TropicalCycloneTrack:
    def __init__(self):
        self.gdf = GeoDataFrame()
        self.epsg = 4326
        self.unit_intensity = "knots"
        self.unit_radii = "nm"

    def read(self, filename, fmt):
        if fmt == "ddb_cyc":
            config = self.read_ddb_cyc(filename)
            return config
        elif fmt == "trk":
            self.read_trk(filename)

    # 3A. convert_units_imperial_metric
    def convert_units_imperial_to_metric(self, config):
        # Convert wind speeds
        # from  knots   - typically 1-minute averaged
        # to    m/s     - we account for conversion here        
        if (self.unit_intensity == "knots") and (self.unit_radii == "nm"):

            # Convert wind radii
            for it in range(len(self.gdf)):
                if self.gdf.vmax[it] > 0.0:
                    self.gdf.vmax[it] = self.gdf.vmax[it] * knots_to_ms * config["wind_conversion_factor"]
                else:
                    self.gdf.vmax[it] = np.NaN    
                if self.gdf.rmw[it] > 0.0:
                    self.gdf.rmw[it] = self.gdf.rmw[it] * nm_to_km
                else:
                    self.gdf.rmw[it] = np.NaN    
                # R35
                if self.gdf.r35_ne[it] > 0.0:
                    self.gdf.r35_ne[it] = self.gdf.r35_ne[it] * nm_to_km
                else:
                    self.gdf.r35_ne[it] = np.NaN

                if self.gdf.r35_se[it] > 0.0:
                    self.gdf.r35_se[it] = self.gdf.r35_se[it] * nm_to_km
                else:
                    self.gdf.r35_se[it] = np.NaN

                if self.gdf.r35_sw[it] > 0.0:
                    self.gdf.r35_sw[it] = self.gdf.r35_sw[it] * nm_to_km
                else:
                    self.gdf.r35_sw[it] = np.NaN

                if self.gdf.r35_nw[it] > 0.0:
                    self.gdf.r35_nw[it] = self.gdf.r35_nw[it] * nm_to_km
                else:
                    self.gdf.r35_nw[it] = np.NaN

                # r50
                if self.gdf.r50_ne[it] > 0.0:
                    self.gdf.r50_ne[it] = self.gdf.r50_ne[it] * nm_to_km
                else:
                    self.gdf.r50_ne[it] = np.NaN

                if self.gdf.r50_se[it] > 0.0:
                    self.gdf.r50_se[it] = self.gdf.r50_se[it] * nm_to_km
                else:
                    self.gdf.r50_se[it] = np.NaN

                if self.gdf.r50_sw[it] > 0.0:
                    self.gdf.r50_sw[it] = self.gdf.r50_sw[it] * nm_to_km
                else:
                    self.gdf.r50_sw[it] = np.NaN

                if self.gdf.r50_nw[it] > 0.0:
                    self.gdf.r50_nw[it] = self.gdf.r50_nw[it] * nm_to_km
                else:
                    self.gdf.r50_nw[it] = np.NaN

                # r65
                if self.gdf.r65_ne[it] > 0.0:
                    self.gdf.r65_ne[it] = self.gdf.r65_ne[it] * nm_to_km
                else:
                    self.gdf.r65_ne[it] = np.NaN

                if self.gdf.r65_se[it] > 0.0:
                    self.gdf.r65_se[it] = self.gdf.r65_se[it] * nm_to_km
                else:
                    self.gdf.r65_se[it] = np.NaN

                if self.gdf.r65_sw[it] > 0.0:
                    self.gdf.r65_sw[it] = self.gdf.r65_sw[it] * nm_to_km
                else:
                    self.gdf.r65_sw[it] = np.NaN

                if self.gdf.r65_nw[it] > 0.0:
                    self.gdf.r65_nw[it] = self.gdf.r65_nw[it] * nm_to_km
                else:
                    self.gdf.r65_nw[it] = np.NaN

                # r100
                if self.gdf.r100_ne[it] > 0.0:
                    self.gdf.r100_ne[it] = self.gdf.r100_ne[it] * nm_to_km
                else:
                    self.gdf.r100_ne[it] = np.NaN

                if self.gdf.r100_se[it] > 0.0:
                    self.gdf.r100_se[it] = self.gdf.r100_se[it] * nm_to_km
                else:
                    self.gdf.r100_se[it] = np.NaN

                if self.gdf.r100_sw[it] > 0.0:
                    self.gdf.r100_sw[it] = self.gdf.r100_sw[it] * nm_to_km
                else:
                    self.gdf.r100_sw[it] = np.NaN

                if self.gdf.r100_nw[it] > 0.0:
                    self.gdf.r100_nw[it] = self.gdf.r100_nw[it] * nm_to_km
                else:
                    self.gdf.r100_nw[it] = np.NaN

            # Done, so set variable
            self.unit_intensity = "ms"
            self.unit_radii = "km"


    def compute_forward_speed(self, config):
        """Computes forward speed and dpcdt and store these in track geodataframe"""

        # Assign variables to geopandas dataframe
        zero_list   = [0.0 for _ in range(len(self.gdf))]
        self.gdf = self.gdf.assign(vtx=zero_list)
        self.gdf = self.gdf.assign(vty=zero_list)
        self.gdf = self.gdf.assign(dpcdt=zero_list)

        # Go over time steps
        for it in range(len(self.gdf)):

            # Get basics
            datetime_it = datetime.strptime(self.gdf.datetime[it], dateformat_module)
            coords_it = self.gdf.geometry[it]
            geofacx = 1
            geofacy = 1

            # Determine geo factors (currently only geo support)
            if self.gdf.crs.name == "WGS 84":
                geofacy = 110540
                geofacx = 111320 * np.cos(coords_it.y * np.pi / 180)

            if it == 0:

                # Forward
                datetime_forward = datetime.strptime(
                    self.gdf.datetime[it + 1], dateformat_module
                )
                coords_forward = self.gdf.geometry[it + 1]
                dt = datetime_forward - datetime_it
                dt = dt.total_seconds()
                dx = (coords_forward.x - coords_it.x) * geofacx
                dy = (coords_forward.y - coords_it.y) * geofacy
                dpc = self.gdf.pc[it + 1] - self.gdf.pc[it]

            elif it == len(self.gdf) - 1:

                # Backward
                datetime_backward = datetime.strptime(
                    self.gdf.datetime[it - 1], dateformat_module
                )
                coords_backward = self.gdf.geometry[it - 1]
                dt = datetime_it - datetime_backward
                dt = dt.total_seconds()
                dx = (coords_it.x - coords_backward.x) * geofacx
                dy = (coords_it.y - coords_backward.y) * geofacy
                dpc = self.gdf.pc[it] - self.gdf.pc[it - 1]

            else:

                # Forward
                datetime_forward = datetime.strptime(
                    self.gdf.datetime[it + 1], dateformat_module
                )
                coords_forward = self.gdf.geometry[it + 1]
                dt1 = datetime_forward - datetime_it
                dt1 = dt1.total_seconds()
                dx1 = (coords_forward.x - coords_it.x) * geofacx
                dy1 = (coords_forward.y - coords_it.y) * geofacy
                dpc1 = self.gdf.pc[it + 1] - self.gdf.pc[it]

                # Backward
                datetime_backward = datetime.strptime(
                    self.gdf.datetime[it - 1], dateformat_module
                )
                coords_backward = self.gdf.geometry[it - 1]
                dt2 = datetime_it - datetime_backward
                dt2 = dt2.total_seconds()
                dx2 = (coords_it.x - coords_backward.x) * geofacx
                dy2 = (coords_it.y - coords_backward.y) * geofacy
                dpc2 = self.gdf.pc[it] - self.gdf.pc[it - 1]

                # Combined yields central differences
                dx = np.mean([dx1, dx2])
                dy = np.mean([dy1, dy2])
                dt = np.mean([dt1, dt2])
                dpc = np.mean([dpc1, dpc2])

            # Check dt so we do not divide by zero (always assume a minimum of 10 min)
            if dt < 600.0:
                dt = 600.0

            # Compute variables
            ux = dx / dt  # speed in meter per seconds
            uy = dy / dt
            dpcdt = dpc / (dt / 3600)  # pressure change per hour

            # Check to limit dpc which happens when pc is not know
            if dpcdt > 100.0 or dpcdt < -100.0 or np.isnan(dpcdt):
                dpcdt = 0.0

            # Save as part of the track
            self.gdf.vtx[it] = ux  # forward speed in x
            self.gdf.vty[it] = uy  # forward speed in y
            self.gdf.dpcdt[it] = dpcdt  # pressure difference in time


    # Support functions for creating spiderweb
    # 1. estimate_missing_values => still assuming imperial system
    def estimate_missing_values(self, config):
        # Go over the track and determine missing values

        for it in range(len(self.gdf)):
            # Get coordinates
            coords_it = self.gdf.geometry[it]

            # determine wind speed
            if np.isnan(self.gdf.vmax[it]):
                if config["wind_pressure_relation"] == "holland2008":
                    # estimate this: vmax is in m/s
                    vmax = wpr_holland2008(
                        pc=self.gdf.pc[it],
                        pn=config["background_pressure"],
                        phi=coords_it.y,
                        dpcdt=self.gdf.dpcdt[it],
                        vt=np.sqrt(self.gdf.vtx[it] ** 2 + self.gdf.vty[it] ** 2),
                        rhoa=config["rho_air"],
                    )

                    # # place this
                    # if self.unit_intensity == "knots":
                    #     self.gdf.vmax[it] = vmax / knots_to_ms
                    # else:
                    self.gdf.vmax[it] = vmax

            # determine pressure
            if np.isnan(self.gdf.pc[it]):
                if config["wind_pressure_relation"] == "holland2008":
                    # estimate this
                    if self.unit_intensity == "knots":
                        pc = wpr_holland2008(
                            vmax=self.gdf.vmax[it] * knots_to_ms,
                            pn=self.background_pressure,
                            phi=coords_it.y,
                            dpcdt=self.gdf.dpcdt[it],
                            vt=np.sqrt(
                                self.gdf.vtx[it] ** 2 + self.gdf.vty[it] ** 2
                            ),
                            rhoa=self.rho_air,
                        )
                    else:
                        pc = wpr_holland2008(
                            vmax=self.gdf.vmax[it],
                            pn=self.background_pressure,
                            phi=coords_it.y,
                            dpcdt=self.gdf.dpcdt[it],
                            vt=np.sqrt(
                                self.gdf.vtx[it] ** 2 + self.gdf.vty[it] ** 2
                            ),
                            rhoa=self.rho_air,
                        )

                    # place this
                    self.gdf.pc[it] = pc

            # radius of maximum winds (RMW)
            if np.isnan(self.gdf.rmw[it]):
                # Nederhoff et al. 2019
                if config["rmw_relation"] == "nederhoff2019":
                    # Estimate: relationship assumes m/s
                    if self.unit_intensity == "knots":
                        [rmax, dr35] = wind_radii_nederhoff(
                            self.gdf.vmax[it] * knots_to_ms, coords_it.y, 7, 0
                        )
                    else:
                        [rmax, dr35] = wind_radii_nederhoff(
                            self.gdf.vmax[it], coords_it.y, 7, 0
                        )

                    # Place value: output is in km
                    if self.unit_radii == "nm":
                        self.gdf.rmw[it] = rmax["mode"] / nm_to_km
                    else:
                        self.gdf.rmw[it] = rmax["mode"]

                # Gross et al. 2004
                elif config["rmw_relation"] == "gross2004":
                    # Estimate: relationship assume knots
                    if self.unit_radii == "knots":
                        rmax = (
                            35.37
                            - 0.11100 * self.gdf.vmax[it]
                            + 0.5700 * (abs(coords_it.y) - 25)
                        )
                    else:
                        rmax = (
                            35.37
                            - 0.11100 * self.gdf.vmax[it] / knots_to_ms
                            + 0.5700 * (abs(coords_it.y) - 25)
                        )

                    # Place value: output is in nm
                    if self.unit_radii == "nm":
                        self.gdf.rmw[it] = rmax
                    else:
                        self.gdf.rmw[it] = rmax * nm_to_km

                # Simple constant value of 25 nm
                elif config["rmw_relation"] == "constant_25nm":
                    # Place 25 nm
                    if self.unit_radii == "nm":
                        self.gdf.rmw[it] = float(25)
                    else:
                        self.gdf.rmw[it] = float(25) * nm_to_km

            # radius of gale force winds (R35)
            if config["wind_profile"] == "holland2010":
                if self.gdf.vmax[it] >= 35:
                    if (
                        (self.gdf.r35_ne[it] == -999)
                        and (self.gdf.r35_se[it] == -999)
                        and (self.gdf.r35_sw[it] == -999)
                        and (self.gdf.r35_nw[it] == -999)
                    ):
                        # Estimate values
                        if self.unit_radii == "knots":
                            [rmax, dr35] = wind_radii_nederhoff(
                                self.gdf.vmax[it] * knots_to_ms, coords_it.y, 7, 0
                            )
                        else:
                            [rmax, dr35] = wind_radii_nederhoff(
                                self.gdf.vmax[it], coords_it.y, 7, 0
                            )

                        if self.unit_radii == "nm":
                            self.gdf.r35_NE[it] = (
                                dr35["mode"] / nm_to_km + self.gdf.rmw[it]
                            )
                            self.gdf.r35_se[it] = (
                                dr35["mode"] / nm_to_km + self.gdf.rmw[it]
                            )
                            self.gdf.r35_sw[it] = (
                                dr35["mode"] / nm_to_km + self.gdf.rmw[it]
                            )
                            self.gdf.r35_nw[it] = (
                                dr35["mode"] / nm_to_km + self.gdf.rmw[it]
                            )
                        else:
                            self.gdf.r35_NE[it] = dr35["mode"] + self.gdf.rmw[it]
                            self.gdf.r35_se[it] = dr35["mode"] + self.gdf.rmw[it]
                            self.gdf.r35_sw[it] = dr35["mode"] + self.gdf.rmw[it]
                            self.gdf.r35_nw[it] = dr35["mode"] + self.gdf.rmw[it]

    # 2A. cut_off_low_wind_speeds
    def cut_off_low_wind_speeds(self):
        # Only apply this when the cut_off wind is defined
        if self.low_wind_speeds_cut_off > 0.0:
            # Find first
            ifirst = []
            for it in range(len(self.track)):
                if self.track.vmax[it] >= self.low_wind_speeds_cut_off and not ifirst:
                    ifirst = it
                    break
            if ifirst > 0:
                self.track = self.track.drop(list(range(0, ifirst)))
                self.track = self.track.reset_index(drop=True)

            # Find last
            ilast = []
            for it in range(len(self.track) - 1, 0 - 1, -1):
                if self.track.vmax[it] >= self.low_wind_speeds_cut_off and not ilast:
                    ilast = it
                    break
            if ilast:
                self.track = self.track.drop(list(range(ilast + 1, len(self.track))))
                self.track = self.track.reset_index(drop=True)

        else:
            if self.debug == 1:
                print("No cut_off_low_wind_speeds since wind speed is zero or lower")

    # 2B. Extent track with certain number of days
    def adjust_track_extent(self):
        # Only apply this when the extend days are defined
        if self.extend_track > 0.0:
            # Compute last gradient
            it_last = len(self.track) - 1
            coords2 = self.track.geometry[it_last]  # last
            datetime2 = datetime.strptime(
                self.track.datetime[it_last], dateformat_module
            )
            it = len(self.track) - 2
            coords1 = self.track.geometry[it]  # first
            datetime1 = datetime.strptime(self.track.datetime[it], dateformat_module)
            dx = coords2.x - coords1.x
            dy = coords2.y - coords1.y
            dt = datetime2 - datetime1

            # Extending the track
            for i in range(1, int(self.extend_track)):
                # Get location
                dt_factor = 86400 / dt.seconds
                x = coords2.x + dx * dt_factor * i
                y = coords2.y + dy * dt_factor * i
                point = Point(x, y)

                # Make time
                date_format = "%Y%m%d %H%M%S"
                tc_time = datetime.strptime(self.track.datetime[it_last], date_format)
                tc_time = tc_time + timedelta(days=i)
                tc_time_string = tc_time.strftime(date_format)

                # Make GeoDataFrame
                gdf = gpd.GeoDataFrame(
                    {
                        "datetime": [tc_time_string],
                        "geometry": [point],
                        "vmax": [self.track.vmax[it_last]],
                        "pc": [self.track.pc[it_last]],
                        "RMW": [self.track.RMW[it_last]],
                        "R35_NE": [0],
                        "R35_SE": [0],
                        "R35_SW": [0],
                        "R35_NW": [0],
                        "R50_NE": [0],
                        "R50_SE": [0],
                        "R50_SW": [0],
                        "R50_NW": [0],
                        "R65_NE": [0],
                        "R65_SE": [0],
                        "R65_SW": [0],
                        "R65_NW": [0],
                        "R100_NE": [0],
                        "R100_SE": [0],
                        "R100_SW": [0],
                        "R100_NW": [0],
                    }
                )
                gdf.set_crs(epsg=self.EPSG, inplace=True)

                # Append self
                self.track = pd.concat([self.track, gdf])

            # Done
            self.track = self.track.reset_index(drop=True)
            if self.debug == 1:
                print("Successfully extended track")

        else:
            if self.debug == 1:
                print("No extending since number of days is zero or lower")

    def resample(self, dt, method="spline"):
        """Resample the track to a new time step. By default uses a spline interpolation for the track points."""
        
        # Get the first and last time
        t0 = datetime.strptime(self.gdf.datetime[0], dateformat_module)
        t1 = datetime.strptime(self.gdf.datetime[len(self.gdf) - 1], dateformat_module)
        dt = timedelta(hours=dt)
        track1 = {}
        # Determine number of time steps
        nt = int((t1 - t0) / dt) + 1
        track1["datetime"] = [t0 + x * dt for x in range(nt)]

        track0 = {}
        track0["datetime"] = []
        track0["x"]        = []
        track0["y"]        = []
        # Loop through columns in gdf
        for column in self.gdf.items():
            if column[0] == "datetime":
                continue
            elif column[0] == "geometry":
                continue
            else:
                track0[column[0]] = []

        # Now loop through the times
        for index, row in self.gdf.iterrows():
            # Convert to datetime
            track0["datetime"].append(datetime.strptime(row.datetime, dateformat_module))
            # Get x and y
            track0["x"].append(row.geometry.x)
            track0["y"].append(row.geometry.y)
            # Loop through columns in gdf
            for column in self.gdf.items():
                if column[0] == "datetime":
                    continue
                elif column[0] == "geometry":
                    continue
                else:
                    track0[column[0]].append(row[column[0]])

        # Convert to timestamp for interpolation
        track0["timestamp"] = [date.timestamp() for date in track0["datetime"]]
        track1["timestamp"] = [date.timestamp() for date in track1["datetime"]]

        # Iterate through items in track0
        for key, value in track0.items():
            if key == "datetime" or key == "timestamp":
                continue
            elif key == "x" or key == "y":
                if method == "linear":
                    f = interp1d(track0["timestamp"], value)
                else:
                    f = CubicSpline(track0["timestamp"], value)
                track1[key] = f(track1["timestamp"])
            else:
                # Interpolate linearly
                f = interp1d(track0["timestamp"], value)
                track1[key] = f(track1["timestamp"])

        # And now create a new gdf
        gdf = GeoDataFrame()
        for index, time in enumerate(track1["datetime"]):
            point = Point(track1["x"][index], track1["y"][index])
            gdf_point = GeoDataFrame()
            gdf_point["geometry"] = [Point(track1["x"][index], track1["y"][index])]  # Assign the new geometry
            gdf_point["datetime"] = track1["datetime"][index].strftime(dateformat_module)
            # Loop through items in track1
            for key, value in track1.items():
                if key == "datetime" or key == "timestamp" or key == "x" or key == "y":
                    continue
                else:
                    gdf_point[key] = track1[key][index]
            gdf = pd.concat([gdf, gdf_point], ignore_index=True)        

        crs = self.gdf.crs
        self.gdf = gdf
        self.gdf.set_crs(crs, inplace=True)



    # # 3B. convert_units_metric_imperial
    # def convert_units_metric_imperial(self):
    #     # Convert wind speeds
    #     # from  m/s     - we account for conversion here
    #     # from  knots   - typically 1-minute averaged

    #     if (self.unit_intensity == "ms") and (self.unit_radii == "km"):
    #         # Intensity first
    #         self.track.vmax = (
    #             self.track.vmax / knots_to_ms / self.wind_conversion_factor
    #         )

    #         # Convert radius of maximum winds
    #         self.track.RMW = self.track.RMW / nm_to_km

    #         # Convert wind radii
    #         for it in range(len(self.track)):
    #             # R35
    #             if self.track.R35_NE[it] > 0.0:
    #                 self.track.R35_NE[it] = self.track.R35_NE[it] / nm_to_km
    #             else:
    #                 self.track.R35_NE[it] = -999.0

    #             if self.track.R35_SE[it] > 0.0:
    #                 self.track.R35_SE[it] = self.track.R35_SE[it] / nm_to_km
    #             else:
    #                 self.track.R35_SE[it] = -999.0

    #             if self.track.R35_SW[it] > 0:
    #                 self.track.R35_SW[it] = self.track.R35_SW[it] / nm_to_km
    #             else:
    #                 self.track.R35_SW[it] = -999.0

    #             if self.track.R35_NW[it] > 0.0:
    #                 self.track.R35_NW[it] = self.track.R35_NW[it] / nm_to_km
    #             else:
    #                 self.track.R35_NW[it] = -999.0

    #             # R50
    #             if self.track.R50_NE[it] > 0.0:
    #                 self.track.R50_NE[it] = self.track.R50_NE[it] / nm_to_km
    #             else:
    #                 self.track.R50_NE[it] = -999.0

    #             if self.track.R50_SE[it] > 0.0:
    #                 self.track.R50_SE[it] = self.track.R50_SE[it] / nm_to_km
    #             else:
    #                 self.track.R50_SE[it] = -999.0

    #             if self.track.R50_SW[it] > 0.0:
    #                 self.track.R50_SW[it] = self.track.R50_SW[it] / nm_to_km
    #             else:
    #                 self.track.R50_SW[it] = -999.0

    #             if self.track.R50_NW[it] > 0.0:
    #                 self.track.R50_NW[it] = self.track.R50_NW[it] / nm_to_km
    #             else:
    #                 self.track.R50_NW[it] = -999.0

    #             # R65
    #             if self.track.R65_NE[it] > 0.0:
    #                 self.track.R65_NE[it] = self.track.R65_NE[it] / nm_to_km
    #             else:
    #                 self.track.R65_NE[it] = -999.0

    #             if self.track.R65_SE[it] > 0.0:
    #                 self.track.R65_SE[it] = self.track.R65_SE[it] / nm_to_km
    #             else:
    #                 self.track.R65_SE[it] = -999.0

    #             if self.track.R65_SW[it] > 0.0:
    #                 self.track.R65_SW[it] = self.track.R65_SW[it] / nm_to_km
    #             else:
    #                 self.track.R65_SW[it] = -999.0

    #             if self.track.R65_NW[it] > 0.0:
    #                 self.track.R65_NW[it] = self.track.R65_NW[it] / nm_to_km
    #             else:
    #                 self.track.R65_NW[it] = -999.0

    #             # R100
    #             if self.track.R100_NE[it] > 0.0:
    #                 self.track.R100_NE[it] = self.track.R100_NE[it] / nm_to_km
    #             else:
    #                 self.track.R100_NE[it] = -999.0

    #             if self.track.R100_SE[it] > 0.0:
    #                 self.track.R100_SE[it] = self.track.R100_SE[it] / nm_to_km
    #             else:
    #                 self.track.R100_SE[it] = -999.0

    #             if self.track.R100_SW[it] > 0.0:
    #                 self.track.R100_SW[it] = self.track.R100_SW[it] / nm_to_km
    #             else:
    #                 self.track.R100_SW[it] = -999.0

    #             if self.track.R100_NW[it] > 0.0:
    #                 self.track.R100_NW[it] = self.track.R100_NW[it] / nm_to_km
    #             else:
    #                 self.track.R100_NW[it] = -999.0

    #         # Done, so set variable
    #         self.unit_intensity = "knots"
    #         self.unit_radii = "nm"
    #         if self.debug == 1:
    #             print("convert units to imperial system")

    #     else:
    #         if self.debug == 1:
    #             print("units are already in the imperial system: no action")

    # Providing the track (from other source)
    def provide_track(
        self, datetimes=[], lons=[], lats=[], winds=[], pressures=[], rmw=[], r35=[]
    ):
        # All variables are Python lists and we simply place them in the right structure
        # note that value -999 is a NaN for this module
        # we assume datetimes but convert them internally
        # r35 is a matrix with time and NE, SE, sw, nw Symmetric_Circle

        # Loop over the list and place them
        for index, vmax in enumerate(winds):
            # Get time ready
            date_format = "%Y%m%d %H%M%S"
            tc_time = datetimes[index]
            tc_time_string = tc_time.strftime(date_format)

            # Make GeoDataFrame
            point = Point(lons[index], lats[index])
            gdf = gpd.GeoDataFrame(
                {
                    "datetime": [tc_time_string],
                    "geometry": [point],
                    "vmax": [vmax],
                    "pc": [pressures[index]],
                    "RMW": [rmw[index]],
                    "R35_NE": [r35[index][0]],
                    "R35_SE": [r35[index][1]],
                    "R35_sw": [r35[index][2]],
                    "R35_nw": [r35[index][3]],
                    "R50_NE": [-999],
                    "R50_SE": [-999],
                    "R50_sw": [-999],
                    "R50_nw": [-999],
                    "R65_NE": [-999],
                    "R65_SE": [-999],
                    "R65_sw": [-999],
                    "R65_nw": [-999],
                    "R100_NE": [-999],
                    "R100_SE": [-999],
                    "R100_sw": [-999],
                    "R100_nw": [-999],
                }
            )
            gdf.set_crs(epsg=self.EPSG, inplace=True)

            # Append self
            self.track = pd.concat([self.track, gdf])

        # Done with this
        self.track = self.track.reset_index(drop=True)
        self.track = self.track.drop([0])  # remove the dummy
        self.track = self.track.reset_index(drop=True)


    def delete_point(self, index=None, time=None):
        # Remove point from the track
        pass    

    def interpolate(time=None, dt=None, t0=None, t1=None):
        pass

    def add_point(self, time, lon, lat, vmax=-999.0, pc=-999.0, rmw=-999.0):
        # Insert point at the right place, depending on the time        
        # point = Point(lon, lat)
        # gdf = gpd.GeoDataFrame(
        #     {
        #         "datetime": [time],
        #         "geometry": [point],
        #         "vmax": [vmax],
        #         "pc": [pc],
        #         "RMW": [rmw],
        #         "R35_NE": [-999],
        #         "R35_SE": [-999],
        #         "R35_SW": [-999],
        #         "R35_NW": [-999],
        #         "R50_NE": [-999],
        #         "R50_SE": [-999],
        #         "R50_SW": [-999],
        #         "R50_NW": [-999],
        #         "R65_NE": [-999],
        #         "R65_SE": [-999],
        #         "R65_SW": [-999],
        #         "R65_NW": [-999],
        #         "R100_NE": [-999],
        #         "R100_SE": [-999],
        #         "R100_SW": [-999],
        #         "R100_NW": [-999],
        #     }
        # )
        self.gdf.set_crs(epsg=self.EPSG, inplace=True)
        # pass  

    def read_ddb_config(self, filename):

        config = {}

        # Read all the lines first
        with open(filename, "r") as f:
            lines = f.readlines()

        # Define the name first
        for line in lines:
            if line[0:4] == "Name":
                string_value = line[5:]
                string_value = "".join(ch for ch in string_value if ch.isalnum())
                config[name] = string_value

        # Define other variables names (if they exist)
        for i in range(len(lines)):
            line = lines[i]
            if line[0:11] == "WindProfile":
                string_value = line[23:]
                string_value = "".join(ch for ch in string_value if ch.isalnum())
                config[wind_profile] = string_value
            if line[0:20] == "WindPressureRelation":
                string_value = line[23:]
                string_value = "".join(ch for ch in string_value if ch.isalnum())
                config[wind_pressure_relation] = string_value
            if line[0:12] == "RMaxRelation":
                string_value = line[23:]
                string_value = "".join(ch for ch in string_value if ch.isalnum())
                config[rmw_relation] = string_value
            if line[0:18] == "Backgroundpressure":
                string_value = line[23:]
                config[background_pressure] = float(string_value)
            if line[0:9] == "PhiSpiral":
                string_value = line[23:]
                config[phi_spiral] = float(string_value)
            if line[0:20] == "WindConversionFactor":
                string_value = line[23:]
                config[wind_conversion_factor] = float(string_value)
            if line[0:15] == "SpiderwebRadius":
                string_value = line[23:]
                config[spiderweb_radius] = float(string_value)
            if line[0:12] == "NrRadialBins":
                string_value = line[23:]
                config[nr_radial_bins] = int(string_value)
            if line[0:17] == "NrDirectionalBins":
                string_value = line[23:]
                config[nr_directional_bins] = int(string_value)

        # Read the track
        for i in range(len(lines)):
            line = lines[i]
            if line[0:15] == "##    Datetime " or line[0:15] == "#   Date   Time":
                break

        # Place coordinates in Tropical Cyclone Track
        for j in range(i + 2, len(lines)):
            # Get values
            line = lines[j]
            line = line.split()
            date_format = "%Y%m%d %H%M%S"
            date_string = line[0] + " " + line[1]
            tc_time = datetime.strptime(date_string, date_format)
            tc_time_string = tc_time.strftime(date_format)
            y = float(line[2])
            x = float(line[3])
            vmax = float(line[4])
            pc = float(line[5])
            RMW = float(line[6])
            R35_NE = float(line[7])
            R35_SE = float(line[8])
            R35_SW = float(line[9])
            R35_NW = float(line[10])
            R50_NE = float(line[11])
            R50_SE = float(line[12])
            R50_SW = float(line[13])
            R50_NW = float(line[14])
            R65_NE = float(line[15])
            R65_SE = float(line[16])
            R65_SW = float(line[17])
            R65_NW = float(line[18])
            R100_NE = float(line[19])
            R100_SE = float(line[20])
            R100_SW = float(line[21])
            R100_NW = float(line[22])

            # Make GeoDataFrame
            point = Point(x, y)
            gdf = GeoDataFrame(
                {
                    "datetime": [tc_time_string],
                    "geometry": [point],
                    "vmax": [vmax],
                    "pc": [pc],
                    "RMW": [RMW],
                    "R35_NE": [R35_NE],
                    "R35_SE": [R35_SE],
                    "R35_SW": [R35_SW],
                    "R35_NW": [R35_NW],
                    "R50_NE": [R50_NE],
                    "R50_SE": [R50_SE],
                    "R50_SW": [R50_SW],
                    "R50_NW": [R50_NW],
                    "R65_NE": [R65_NE],
                    "R65_SE": [R65_SE],
                    "R65_SW": [R65_SW],
                    "R65_NW": [R65_NW],
                    "R100_NE": [R100_NE],
                    "R100_SE": [R100_SE],
                    "R100_SW": [R100_SW],
                    "R100_NW": [R100_NW],
                }
            )
            gdf.set_crs(epsg=self.EPSG, inplace=True)

            # Append self
            self.track = pd.concat([self.track, gdf])

        # Done with this
        self.gdf = self.gdf.reset_index(drop=True)
        self.gdf = self.gdf.drop([0])  # remove the dummy
        self.gdf = self.gdf.reset_index(drop=True)

    def read_trk(self, filename):

        # Initialize variables
        lasttime = np.datetime64(datetime.strptime("2000010118", "%Y%m%d%H"))

        # Create a new empty GDF
        self.gdf = GeoDataFrame()

        # Open the file
        with open(filename, "r") as fid:

            # Read line by line
            for s0 in fid:
                s0 = s0.strip()
                if not s0:
                    continue

                vmax = np.nan
                pc   = np.nan
                rmax = np.nan
                r35_ne = np.nan
                r35_se = np.nan
                r35_sw = np.nan
                r35_nw = np.nan
                r50_ne = np.nan
                r50_se = np.nan
                r50_sw = np.nan
                r50_nw = np.nan
                r65_ne = np.nan
                r65_se = np.nan
                r65_sw = np.nan
                r65_nw = np.nan
                r100_ne = np.nan
                r100_se = np.nan
                r100_sw = np.nan
                r100_nw = np.nan

                # Start
                s = s0.split(",")
                tc = {"name": "not_known"}
                tc["basin"] = s[0]
                tc["storm_number"] = int(s[1])

                # Read time
                tstr = str(s[2])
                tstr = tstr.replace(" ", "")
                time = datetime.strptime(tstr, "%Y%m%d%H")
                newtime = np.datetime64(time)

                # Add forecasted time
                hrs = float(s[5])
                newtime = newtime + np.timedelta64(int(hrs * 3600), "s")

                # Latitude
                if s[6][-1] == "N":
                    y = 0.1 * float(s[6][:-1])
                else:
                    y = -0.1 * float(s[6][:-1])
                # Longitude
                if s[7][-1] == "E":
                    x = 0.1 * float(s[7][:-1])
                else:
                    x = -0.1 * float(s[7][:-1])

                # vmax
                if float(s[8]) > 0:
                    vmax = float(s[8])

                # Pressure and raddi
                if len(s) > 9:
                    # Pressure
                    if float(s[9]) > 0:
                        pc = float(s[9])

                    # Radii
                    r = float(s[11])
                    if r in [34, 35]:
                        r35_ne = float(s[13])
                        r35_se = float(s[14])
                        r35_sw = float(s[15])
                        r35_nw = float(s[16])
                    elif r == 50:
                        r50_ne = float(s[13])
                        r50_se = float(s[14])
                        r50_sw = float(s[15])
                        r50_nw = float(s[16])
                    elif r in [64, 65]:
                        r65_ne = float(s[13])
                        r65_se = float(s[14])
                        r65_sw = float(s[15])
                        r65_nw = float(s[16])
                    elif r == 100:
                        r100_ne = float(s[13])
                        r100_se = float(s[14])
                        r100_sw = float(s[15])
                        r100_nw = float(s[16])

                    # Other things
                    if len(s) > 17:
                        if s[17]:
                            pressure_last_closed_isobar = float(s[17])
                        if s[18]:
                            radius_last_closed_isobar = float(s[18])
                        if s[19]:
                            try:
                                if float(s[19]) > 0:
                                    rmax = float(s[19])
                            except ValueError:
                                # print("Error: Unable to convert to float for RMW")
                                rmax = np.nan
                        if len(s) >= 28:
                            if s[27]:
                                name = s[27]

                if newtime > lasttime + np.timedelta64(1, "s"):
                    # New time point found
                    # create new gdf
                    new_point = True
                    gdf_point = GeoDataFrame()
                    gdf_point["geometry"] = [Point(x, y)]  # Assign the new geometry
                    gdf_point["datetime"] = newtime.astype("O").strftime("%Y%m%d %H%M%S")
                    gdf_point["vmax"]   = np.NaN
                    gdf_point["pc"]     = np.NaN
                    gdf_point["rmw"]    = np.NaN
                    gdf_point["r35_ne"] = np.NaN
                    gdf_point["r35_se"] = np.NaN
                    gdf_point["r35_sw"] = np.NaN
                    gdf_point["r35_nw"] = np.NaN
                    gdf_point["r50_ne"] = np.NaN
                    gdf_point["r50_se"] = np.NaN
                    gdf_point["r50_sw"] = np.NaN
                    gdf_point["r50_nw"] = np.NaN
                    gdf_point["r65_ne"] = np.NaN
                    gdf_point["r65_se"] = np.NaN
                    gdf_point["r65_sw"] = np.NaN
                    gdf_point["r65_nw"] = np.NaN
                    gdf_point["r100_ne"] = np.NaN
                    gdf_point["r100_se"] = np.NaN
                    gdf_point["r100_sw"] = np.NaN
                    gdf_point["r100_nw"] = np.NaN
                else:
                    new_point = False
                        
                # Update all variables in gdf (when it is not the first one)                
                gdf_point["vmax"] = vmax
                gdf_point["pc"]   = pc
                gdf_point["rmw"]  = rmax
                if r35_ne > 0.0:
                    gdf_point["r35_ne"] = r35_ne
                if r35_se > 0.0:
                    gdf_point["r35_se"] = r35_se
                if r35_sw > 0.0:
                    gdf_point["r35_sw"] = r35_sw
                if r35_nw > 0.0:
                    gdf_point["r35_nw"] = r35_nw
                if r50_ne > 0.0:
                    gdf_point["r50_ne"] = r50_ne
                if r50_se > 0.0:
                    gdf_point["r50_se"] = r50_se
                if r50_sw > 0.0:
                    gdf_point["r50_sw"] = r50_sw
                if r50_nw > 0.0:
                    gdf_point["r50_nw"] = r50_nw
                if r65_ne > 0.0:
                    gdf_point["r65_ne"] = r65_ne
                if r65_se > 0.0:
                    gdf_point["r65_se"] = r65_se
                if r65_sw > 0.0:
                    gdf_point["r65_sw"] = r65_sw
                if r65_nw > 0.0:
                    gdf_point["r65_nw"] = r65_nw
                if r100_ne > 0.0:
                    gdf_point["r100_ne"] = r100_ne
                if r100_se > 0.0:
                    gdf_point["r100_se"] = r100_se
                if r100_sw > 0.0:
                    gdf_point["r100_sw"] = r100_sw
                if r100_nw > 0.0:
                    gdf_point["r100_nw"] = r100_nw

                if new_point:
                    # Add new point
                    self.gdf = pd.concat([self.gdf, gdf_point])
                else:
                    # Replace row in dataframe with row from other dataframe
                    self.gdf.iloc[-1] = gdf_point.iloc[0]

                lasttime = newtime

        # Done with this => rest so track looks good
        self.gdf = self.gdf.reset_index(drop=True)
        self.gdf.set_crs(epsg=self.epsg, inplace=True)

    def write(self, filename, fmt="ddb_cyc"):
        # If ddb_cyc
        if fmt == "ddb_cyc":
            # Open file
            with open(filename, "wt") as f:
                # # Print header
                # f.writelines(
                #     "# Tropical Cyclone Toolbox - Coastal Hazards Toolkit - "
                #     + self.creation_time.strftime(dateformat_module)
                #     + "\n"
                # )

                # # Print rest
                # f.writelines('Name                   "' + self.name + '"\n')
                # f.writelines("WindProfile            " + self.wind_profile + "\n")
                # f.writelines(
                #     "WindPressureRelation   " + self.wind_pressure_relation + "\n"
                # )
                # f.writelines("RMaxRelation           " + self.rmw_relation + "\n")
                # f.writelines(
                #     "Backgroundpressure     " + str(self.background_pressure) + "\n"
                # )
                # f.writelines("PhiSpiral              " + str(self.phi_spiral) + "\n")
                # f.writelines(
                #     "WindConversionFactor   " + str(self.wind_conversion_factor) + "\n"
                # )
                # f.writelines(
                #     "SpiderwebRadius        " + str(self.spiderweb_radius) + "\n"
                # )
                # f.writelines(
                #     "NrRadialBins           " + str(self.nr_radial_bins) + "\n"
                # )
                # f.writelines(
                #     "NrDirectionalBins      " + str(self.nr_directional_bins) + "\n"
                # )
                # epsg = self.track.crs.name
                # f.writelines("EPSG                   " + epsg + "\n")
                # f.writelines(
                #     "UnitIntensity          " + str(self.unit_intensity) + "\n"
                # )
                # f.writelines("UnitWindRadii          " + str(self.unit_radii) + "\n")

                # # Print header for the track
                # f.writelines("#  \n")
                f.writelines(
                    "#   Date   Time               Lat        Lon         Vmax       Pc          Rmax         R35(NE)      R35(SE)     R35(SW)     R35(NW)     R50(NE)     R50(SE)    R50(SW)    R50(NW)     R65(NE)     R65(SE)     R65(SW)     R65(NW)    R100(NE)    R100(SE)    R100(SW)    R100(NE)  \n"
                )
                f.writelines("#  \n")

                # Print the actual track
                for i in range(len(self.gdf)):

                    f.writelines(self.gdf.datetime[i].rjust(20))
                    coords = self.gdf.geometry[i]
                    f.writelines(str(round(coords.y, 2)).rjust(12))
                    f.writelines(str(round(coords.x, 2)).rjust(12))

                    f.writelines(str(round(self.gdf.vmax[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.pc[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.rmw[i], 1)).rjust(12))

                    f.writelines(str(round(self.gdf.r35_ne[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r35_se[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r35_sw[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r35_nw[i], 1)).rjust(12))

                    f.writelines(str(round(self.gdf.r50_ne[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r50_se[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r50_sw[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r50_nw[i], 1)).rjust(12))

                    f.writelines(str(round(self.gdf.r65_ne[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r65_se[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r65_sw[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r65_nw[i], 1)).rjust(12))

                    f.writelines(str(round(self.gdf.r100_ne[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r100_se[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r100_sw[i], 1)).rjust(12))
                    f.writelines(str(round(self.gdf.r100_nw[i], 1)).rjust(12))

                    f.writelines("\n")

            # if self.debug == 1:
            #     print("Successfully written track - ddb_cyc")
        else:
            print(
                'For other methods of writing the track; please used the "tc.track.to_file" option'
            )

    def to_gdf(self, filename=None):
        """Make track GeoDataFrame and optionally write to file"""

        categories = {64.0: "TS", 83.0: "1", 96.0: "2", 113.0: "3", 137.0: "4"}

        features = []
        points = []

        # First the track
        for ip in range(np.size(self.gdf.geometry.x)):
            points.append([self.gdf.geometry.x[ip], self.gdf.geometry.y[ip]])
        features.append(Feature(geometry=LineString(coordinates=points), properties={"name": "No name"}))

        # Then the points
        for ip in range(np.size(self.gdf.geometry.x)):
            point = Point((self.gdf.geometry.x[ip], self.gdf.geometry.y[ip]))
            tmptime = datetime.strptime(self.gdf.datetime[ip], "%Y%m%d %H%M%S")
            vmax = self.gdf.vmax[ip]
            for threshold, cat in categories.items():
                if vmax < threshold:
                    break
            else:
                cat = "5"
            features.append(
                Feature(
                    geometry=point,
                    properties={
                        "time": tmptime.strftime("%Y/%m/%d %H:%M") + " UTC",
                        "lon": self.gdf.geometry.x[ip],
                        "lat": self.gdf.geometry.y[ip],
                        "vmax": self.gdf.vmax[ip],
                        "pc": self.gdf.pc[ip],
                        "category": cat,
                    },
                )
            )

        gdf = GeoDataFrame().from_features(FeatureCollection(features))

        # Write to file 
        if filename is not None:
            # Get extension of filename
            ext = os.path.splitext(filename)[-1]
            # Get path of filename
            folder_path = os.path.dirname(filename)
            # Make path (if needed)
            if folder_path != "":
                os.makedirs(folder_path, exist_ok=True)
            if ext == ".shp":
                gdf.to_file(filename)
            elif ext == ".geojson":
                gdf.to_file(filename, driver="GeoJSON")
            elif ext == ".js":
                gdf_to_geojson_js(gdf, filename, varname="track_ensemble")   

        return gdf    
