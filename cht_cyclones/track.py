import os
from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from geojson import Feature, FeatureCollection
from geopandas import GeoDataFrame
from pyproj import CRS, Geod
from scipy.interpolate import CubicSpline, interp1d
from shapely.geometry import LineString, MultiLineString, Point, mapping

from .utils import gdf_to_geojson_js, gdf_to_pli
from .wind_profiles import wind_radii_nederhoff, wpr_holland2008

knots_to_ms = float(0.51444)
nm_to_km = float(1.852)
nm_to_m = float(1.852) * 1000
dateformat_module = "%Y%m%d %H%M%S"

geodesic = Geod(ellps="WGS84")


class TropicalCycloneTrack:
    def __init__(self):
        self.gdf = GeoDataFrame()
        self.crs = CRS.from_epsg(4326)
        self.unit_intensity = "knots"
        self.unit_radii = "nm"

    def read(self, filename, tau=0):
        """Read a tropical cyclone track from a file."""

        # Try to read cyc file new format
        try:
            gdf = read_cyc(filename)
            self.gdf = gdf
            self.apply_tau(tau)
            return None
        except:
            pass

        # Try to read trk file (COAMPS-TC)
        try:
            gdf = read_trk(filename)
            self.gdf = gdf
            # Get rid of NaN values in vmax
            self.fix_nan_vmax()
            self.apply_tau(tau)
            return None
        except:
            pass

        # Try to read ddb cyc file
        try:
            gdf, config, name = read_ddb_cyc(filename)
            self.gdf = gdf
            self.apply_tau(tau)
            # Do not do anything with name
            return config
        except:
            pass

        raise Exception("Error: could not read file " + filename)

    def apply_tau(self, tau):
        """Apply tau to the track, i.e. start later if timec in track is smaller than tau"""
        if tau > 0:
            tstart = datetime.strptime(self.gdf.datetime[0], dateformat_module) + timedelta(hours=tau)
            self.shorten(tstart=tstart)

    def fix_nan_vmax(self):
        # Sometimes COAMPS-TC data has 0.0 for vmax, which is replaced by NaN
        # This is causes problems later on. Try to fix it by doing an interpolation where necessary. If it happens at the beginning or end, use nearest neighbor.
        okay = True
        for it in range(len(self.gdf)):
            if np.isnan(self.gdf.vmax[it]):
                okay = False
                break

        if not okay:
            # Make pandas series
            s = pd.Series(self.gdf.vmax)
            # Interpolate
            s = s.interpolate(method="linear")
            # Put values back in gdf
            for it in range(len(self.gdf)):
                if np.isnan(self.gdf.vmax[it]):
                    print("Warning! Replacing vmax = NaN in track with interpolated value at it = " + str(it))
                    self.gdf.loc[it, "vmax"] = s.values[it]

    def write(self,
              filename,
              format="cyc",
              config=None,
              name=None,
              include_header=True):

        # Make a copy of the track and replace NaN with -999.0
        gdf = self.gdf.fillna(-999.0)

        if format == "ddb_cyc":
            if config is None:
                print("Error: no configuration file provided")
                return
            write_ddb_cyc(filename, gdf, config, name, self.unit_intensity, self.unit_radii)
        elif format == "cyc":
            write_cyc(filename, gdf, include_header=include_header)
        else:
            raise Exception("Error: track format " + format + " not supported")

    # 3A. convert_units_imperial_metric
    def convert_units_imperial_to_metric(self, config):
        # Convert wind speeds
        # from  knots   - typically 1-minute averaged
        # to    m/s     - we account for conversion here        
        if (self.unit_intensity == "knots") and (self.unit_radii == "nm"):

            # Convert wind radii
            for it in range(len(self.gdf)):
                if self.gdf.vmax[it] > 0.0:
                    self.gdf.loc[it, "vmax"] = self.gdf.vmax[it] * knots_to_ms * config["wind_conversion_factor"]
                else:
                    self.gdf.loc[it, "vmax"] = np.nan    
                if self.gdf.rmw[it] > 0.0:
                    self.gdf.loc[it, "rmw"] = self.gdf.rmw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "rmw"] = np.nan    
                # R35
                if self.gdf.r35_ne[it] > 0.0:
                    self.gdf.loc[it, "r35_ne"] = self.gdf.r35_ne[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r35_ne"] = np.nan

                if self.gdf.r35_se[it] > 0.0:
                    self.gdf.loc[it, "r35_se"] = self.gdf.r35_se[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r35_se"] = np.nan

                if self.gdf.r35_sw[it] > 0.0:
                    self.gdf.loc[it, "r35_sw"] = self.gdf.r35_sw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r35_sw"] = np.nan

                if self.gdf.r35_nw[it] > 0.0:
                    self.gdf.loc[it, "r35_nw"] = self.gdf.r35_nw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r35_nw"] = np.nan

                # r50
                if self.gdf.r50_ne[it] > 0.0:
                    self.gdf.loc[it, "r50_ne"] = self.gdf.r50_ne[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r50_ne"] = np.nan

                if self.gdf.r50_se[it] > 0.0:
                    self.gdf.loc[it, "r50_se"] = self.gdf.r50_se[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r50_se"] = np.nan

                if self.gdf.r50_sw[it] > 0.0:
                    self.gdf.loc[it, "r50_sw"] = self.gdf.r50_sw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r50_sw"] = np.nan

                if self.gdf.r50_nw[it] > 0.0:
                    self.gdf.loc[it, "r50_nw"] = self.gdf.r50_nw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r50_nw"] = np.nan

                # r65
                if self.gdf.r65_ne[it] > 0.0:
                    self.gdf.loc[it, "r65_ne"] = self.gdf.r65_ne[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r65_ne"]= np.nan

                if self.gdf.r65_se[it] > 0.0:
                    self.gdf.loc[it, "r65_se"] = self.gdf.r65_se[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r65_se"] = np.nan

                if self.gdf.r65_sw[it] > 0.0:
                    self.gdf.loc[it, "r65_sw"] = self.gdf.r65_sw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r65_sw"] = np.nan

                if self.gdf.r65_nw[it] > 0.0:
                    self.gdf.loc[it, "r65_nw"] = self.gdf.r65_nw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r65_nw"] = np.nan

                # r100
                if self.gdf.r100_ne[it] > 0.0:
                    self.gdf.loc[it, "r100_ne"] = self.gdf.r100_ne[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r100_ne"] = np.nan

                if self.gdf.r100_se[it] > 0.0:
                    self.gdf.loc[it, "r100_se"] = self.gdf.r100_se[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r100_se"] = np.nan

                if self.gdf.r100_sw[it] > 0.0:
                    self.gdf.loc[it, "r100_sw"] = self.gdf.r100_sw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r100_sw"] = np.nan

                if self.gdf.r100_nw[it] > 0.0:
                    self.gdf.loc[it, "r100_nw"] = self.gdf.r100_nw[it] * nm_to_km
                else:
                    self.gdf.loc[it, "r100_nw"] = np.nan

            # Done, so set variable
            self.unit_intensity = "ms"
            self.unit_radii = "km"

    def compute_forward_speed(self):
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
                            pn=config["background_pressure"],
                            phi=coords_it.y,
                            dpcdt=self.gdf.dpcdt[it],
                            vt=np.sqrt(
                                self.gdf.vtx[it] ** 2 + self.gdf.vty[it] ** 2
                            ),
                            rhoa=config["rho_air"],
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

    # # 2A. cut_off_low_wind_speeds (should not do this anymore)
    # def cut_off_low_wind_speeds(self):
    #     # Only apply this when the cut_off wind is defined
    #     if self.low_wind_speeds_cut_off > 0.0:
    #         # Find first
    #         ifirst = []
    #         for it in range(len(self.track)):
    #             if self.track.vmax[it] >= self.low_wind_speeds_cut_off and not ifirst:
    #                 ifirst = it
    #                 break
    #         if ifirst > 0:
    #             self.track = self.track.drop(list(range(0, ifirst)))
    #             self.track = self.track.reset_index(drop=True)

    #         # Find last
    #         ilast = []
    #         for it in range(len(self.track) - 1, 0 - 1, -1):
    #             if self.track.vmax[it] >= self.low_wind_speeds_cut_off and not ilast:
    #                 ilast = it
    #                 break
    #         if ilast:
    #             self.track = self.track.drop(list(range(ilast + 1, len(self.track))))
    #             self.track = self.track.reset_index(drop=True)

    #     else:
    #         if self.debug == 1:
    #             print("No cut_off_low_wind_speeds since wind speed is zero or lower")

    def extend(self, hours_before=None, hours_after=None):
        """Extend the track with a certain number of hours before and after the track."""
        if hours_before is not None:
            t0 = datetime.strptime(self.gdf.datetime[0], dateformat_module)
            t = t0 - timedelta(hours=hours_before)
            self.add_point(t)
        if hours_after is not None:
            t1 = datetime.strptime(self.gdf.datetime[len(self.gdf) - 1], dateformat_module)
            t = t1 + timedelta(hours=hours_after)
            self.add_point(t)    

    def shorten(self, tstart=None, tend=None):
        """Shorten the track to a certain time period."""
        times_to_remove = []
        if tstart is not None:
            for it, row in self.gdf.iterrows():
                t = datetime.strptime(row.datetime, dateformat_module)
                if t < tstart:
                    times_to_remove.append(t)
                    # self.remove_point(time=t)
        if tend is not None:
            for it, row in self.gdf.iterrows():
                t = datetime.strptime(row.datetime, dateformat_module)
                if t > tend:
                    times_to_remove.append(t)
                    # self.remove_point(time=t)            
        if len(times_to_remove) > 0:
            for t in times_to_remove:
                self.remove_point(time=t)

    def resample(self, dt, method="spline"):

        """Resample the track to a new time step. Returns new gdf. By default uses a spline interpolation for the track points."""
        
        # Get the first and last time
        t0 = datetime.strptime(self.gdf.datetime[0], dateformat_module)
        t1 = datetime.strptime(self.gdf.datetime[len(self.gdf) - 1], dateformat_module)
        dt = timedelta(hours=dt)

        # Determine number of time steps
        nt = int((t1 - t0) / dt) + 1

        # Create a new track gdf
        gdf = GeoDataFrame()

        # Copy first point of existing track
        for it in range(nt):
            gdf = pd.concat([gdf, self.gdf.iloc[[0]]], ignore_index=True)
            gdf.loc[it, "datetime"] = (t0 + it * dt).strftime(dateformat_module)

        self.gdf = interpolate_track(self.gdf, gdf, method=method)

        return gdf

    def add_point(self, time: datetime,
                  lon: float=-9999.0,
                  lat: float=-9999.0,
                  method: str="spline"):
        """Insert a point at a certain time and location (optional)"""
            
        # Find index where to insert
        # Check first if time is before the first time
        t0 = datetime.strptime(self.gdf.datetime[0], dateformat_module)
        if time < t0:
            index = 0
        else:
            for it, row in self.gdf.iterrows():
                if time == datetime.strptime(row.datetime, dateformat_module):
                    raise Exception("Error: time already exists in track")
                if time > datetime.strptime(row.datetime, dateformat_module):
                    index = it + 1

        if index == 0 or index == len(self.gdf):
            # Need to compute the forward speed to estimate track position
            self.compute_forward_speed()

        gdf_point, index = interpolate_to_point(self.gdf, time, method=method)

        # Check if lon and lat are provided
        if lon == -9999.0 or lat == -9999.0:
            # Not provided so use interpolated values
            pass
        else:
            # Use provided coordinates
            point = Point(lon, lat)
            gdf_point["geometry"] = [point]

        gdf_point.replace(-999.0, np.nan, inplace=True)
        crs = self.gdf.crs

        if index == 0:
            self.gdf = pd.concat([gdf_point, self.gdf], ignore_index=True)
        elif index == len(self.gdf):
            self.gdf = pd.concat([self.gdf, gdf_point], ignore_index=True)
        else:
            self.gdf = pd.concat([self.gdf.iloc[:index], gdf_point, self.gdf.iloc[index:]], ignore_index=True)            

        self.gdf = self.gdf.set_crs(crs)    
        self.gdf = self.gdf.reset_index(drop=True)

    def remove_point(self, index=None, time=None):
        """Remove point from the track by index or time"""

        if index is not None:
            self.gdf = self.gdf.drop(index)
            self.gdf = self.gdf.reset_index(drop=True)

        if time is not None:
            # Loop through the track and find the index
            for index, row in self.gdf.iterrows():
                t = datetime.strptime(row.datetime, dateformat_module)
                if t == time:
                    self.gdf = self.gdf.drop(index)
                    self.gdf = self.gdf.reset_index(drop=True)
                    break

    def add(self, track1):
        """Add another track to the current track"""

        # If self is empty, just copy track1
        if len(self.gdf) == 0:
            self.gdf = track1.gdf
            return

        # Get first time of track1
        t0 = datetime.strptime(track1.gdf.datetime[0], dateformat_module)

        # Remove all points in self that are at or after t0
        for it, row in self.gdf.iterrows():
            t = datetime.strptime(row.datetime, dateformat_module)
            if t >= t0:
                self.remove_point(time=t)

        # Concatenate the two tracks
        self.gdf = pd.concat([self.gdf, track1.gdf], ignore_index=True).reset_index(drop=True)


    def to_gdf(self, filename=None):
        """Make track GeoDataFrame and optionally write to file"""

        categories = {33.0: "TD", 64.0: "TS", 83.0: "1", 96.0: "2", 113.0: "3", 137.0: "4"}

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
            # if vmax is NaN, set vmax to 1.0
            if np.isnan(vmax):
                vmax = 1.0
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
                gdf_to_geojson_js(gdf, filename, varname="track_data")   
            elif ext == ".pli":
                gdf_to_pli(gdf, filename)   

        return gdf    


def read_cyc(filename):
    # cyc format

    gdf = GeoDataFrame()

    # Read all the lines first
    with open(filename, "r") as f:
        lines = f.readlines()

    for line in lines:
        if line[0] == "#":
            continue
        # Get values
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
        gdf_point = GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [vmax],
                "pc": [pc],
                "rmw": [RMW],
                "r35_ne": [R35_NE],
                "r35_se": [R35_SE],
                "r35_sw": [R35_SW],
                "r35_nw": [R35_NW],
                "r50_ne": [R50_NE],
                "r50_se": [R50_SE],
                "r50_sw": [R50_SW],
                "r50_nw": [R50_NW],
                "r65_ne": [R65_NE],
                "r65_se": [R65_SE],
                "r65_sw": [R65_SW],
                "r65_nw": [R65_NW],
                "r100_ne": [R100_NE],
                "r100_se": [R100_SE],
                "r100_sw": [R100_SW],
                "r100_nw": [R100_NW],
            }
        )

        # Append self
        gdf = pd.concat([gdf, gdf_point])

    # Replace -999.0 with NaN
    gdf = gdf.replace(-999.0, np.nan)
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(crs=CRS(4326), inplace=True)

    return gdf

def read_ddb_cyc(filename):
    # Old ddb cyc format

    gdf = GeoDataFrame()
    config = {}

    # Read all the lines first
    with open(filename, "r") as f:
        lines = f.readlines()

    # Define the name first
    for line in lines:
        if line[0:4] == "Name":
            string_value = line[5:]
            string_value = "".join(ch for ch in string_value if ch.isalnum())
            name = string_value

    # Define other variables names (if they exist)
    for i in range(len(lines)):
        line = lines[i]
        if line[0:11] == "WindProfile":
            string_value = line[23:]
            string_value = "".join(ch for ch in string_value if ch.isalnum())
            config["wind_profile"] = string_value
        if line[0:20] == "WindPressureRelation":
            string_value = line[23:]
            string_value = "".join(ch for ch in string_value if ch.isalnum())
            config["wind_pressure_relation"] = string_value
        if line[0:12] == "RMaxRelation":
            string_value = line[23:]
            string_value = "".join(ch for ch in string_value if ch.isalnum())
            config["rmw_relation"] = string_value
        if line[0:18] == "Backgroundpressure":
            string_value = line[23:]
            config["background_pressure"] = float(string_value)
        if line[0:9] == "PhiSpiral":
            string_value = line[23:]
            config["phi_spiral"] = float(string_value)
        if line[0:20] == "WindConversionFactor":
            string_value = line[23:]
            config["wind_conversion_factor"] = float(string_value)
        if line[0:15] == "SpiderwebRadius":
            string_value = line[23:]
            config["spiderweb_radius"] = float(string_value)
        if line[0:12] == "NrRadialBins":
            string_value = line[23:]
            config["nr_radial_bins"] = int(string_value)
        if line[0:17] == "NrDirectionalBins":
            string_value = line[23:]
            config["nr_directional_bins"] = int(string_value)

    # Read the track and find last comment line
    for i in range(len(lines)):
        line = lines[i]
        if line[0] == "#":
            icomment = i

    # Place coordinates in Tropical Cyclone Track
    for j in range(icomment + 1, len(lines)):
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
        gdf_point = GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [vmax],
                "pc": [pc],
                "rmw": [RMW],
                "r35_ne": [R35_NE],
                "r35_se": [R35_SE],
                "r35_sw": [R35_SW],
                "r35_nw": [R35_NW],
                "r50_ne": [R50_NE],
                "r50_se": [R50_SE],
                "r50_sw": [R50_SW],
                "r50_nw": [R50_NW],
                "r65_ne": [R65_NE],
                "r65_se": [R65_SE],
                "r65_sw": [R65_SW],
                "r65_nw": [R65_NW],
                "r100_ne": [R100_NE],
                "r100_se": [R100_SE],
                "r100_sw": [R100_SW],
                "r100_nw": [R100_NW],
            }
        )

        # Append self
        gdf = pd.concat([gdf, gdf_point])

    # Replace -999.0 with NaN
    gdf = gdf.replace(-999.0, np.nan)
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(crs=CRS(4326), inplace=True)

    return gdf, config, name

def read_trk(filename):
    # Read the track from a TRK file

    # Create a new empty GDF
    gdf = GeoDataFrame()

    # Initialize variables
    lasttime = np.datetime64(datetime.strptime("2000010118", "%Y%m%d%H"))

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
                gdf_point = GeoDataFrame(geometry=[Point(x, y)])
                # gdf_point["geometry"] = [Point(x, y)]  # Assign the new geometry
                gdf_point["datetime"] = newtime.astype("O").strftime("%Y%m%d %H%M%S")
                gdf_point["vmax"]   = np.nan
                gdf_point["pc"]     = np.nan
                gdf_point["rmw"]    = np.nan
                gdf_point["r35_ne"] = np.nan
                gdf_point["r35_se"] = np.nan
                gdf_point["r35_sw"] = np.nan
                gdf_point["r35_nw"] = np.nan
                gdf_point["r50_ne"] = np.nan
                gdf_point["r50_se"] = np.nan
                gdf_point["r50_sw"] = np.nan
                gdf_point["r50_nw"] = np.nan
                gdf_point["r65_ne"] = np.nan
                gdf_point["r65_se"] = np.nan
                gdf_point["r65_sw"] = np.nan
                gdf_point["r65_nw"] = np.nan
                gdf_point["r100_ne"] = np.nan
                gdf_point["r100_se"] = np.nan
                gdf_point["r100_sw"] = np.nan
                gdf_point["r100_nw"] = np.nan
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
                gdf = pd.concat([gdf, gdf_point])
            else:
                # Replace row in dataframe with row from other dataframe
                gdf.iloc[-1] = gdf_point.iloc[0]

            lasttime = newtime

    # Done with this => rest so track looks good
    gdf = gdf.reset_index(drop=True)
    gdf = gdf.set_crs(crs=CRS(4326), inplace=True)

    return gdf

def write_cyc(filename, gdf, include_header=True):   
    """Write to cyc format"""
    with open(filename, "wt") as f:
        # Print header
        if include_header:
            f.writelines(
                "#   Date   Time      Lat      Lon     Vmax       Pc     Rmax  R35(NE)  R35(SE)  R35(SW)  R35(NW)  R50(NE)  R50(SE)  R50(SW)  R50(NW)  R65(NE)  R65(SE)  R65(SW)  R65(NW) R100(NE) R100(SE) R100(SW) R100(NE)\n"
            )

        # Print the actual track
        for i in range(len(gdf)):

            f.writelines(gdf.datetime[i].rjust(15))
            coords = gdf.geometry[i]
            f.writelines(str(round(coords.y, 2)).rjust(9))
            f.writelines(str(round(coords.x, 2)).rjust(9))

            f.writelines(str(round(gdf.vmax[i], 1)).rjust(9))
            f.writelines(str(round(gdf.pc[i], 1)).rjust(9))
            f.writelines(str(round(gdf.rmw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r35_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r50_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r65_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r100_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_nw[i], 1)).rjust(9))

            f.writelines("\n")

def write_ddb_cyc(filename, gdf, config, name, unit_intensity, unit_radii):
    """Write to 'old' Delft Dashboard cyc format"""
    with open(filename, "wt") as f:
        # Print header
        f.writelines(
            "# Tropical Cyclone Toolbox - Coastal Hazards Toolkit - "
            + datetime.now().strftime(dateformat_module)
            + "\n"
        )
        # Print rest
        f.writelines('Name                   "' + name + '"\n')
        f.writelines("WindProfile            " + config["wind_profile"] + "\n")
        f.writelines(
            "WindPressureRelation   " + config["wind_pressure_relation"] + "\n"
        )
        f.writelines("RMaxRelation           " + config["rmw_relation"] + "\n")
        f.writelines(
            "Backgroundpressure     " + str(config["background_pressure"]) + "\n"
        )
        f.writelines("PhiSpiral              " + str(config["phi_spiral"]) + "\n")
        f.writelines(
            "WindConversionFactor   " + str(config["wind_conversion_factor"]) + "\n"
        )
        f.writelines(
            "SpiderwebRadius        " + str(config["spiderweb_radius"]) + "\n"
        )
        f.writelines(
            "NrRadialBins           " + str(config["nr_radial_bins"]) + "\n"
        )
        f.writelines(
            "NrDirectionalBins      " + str(config["nr_directional_bins"]) + "\n"
        )
        f.writelines("EPSG                   WGS84\n")
        f.writelines(
            "UnitIntensity          " + str(unit_intensity) + "\n"
        )
        f.writelines("UnitWindRadii          " + str(unit_radii) + "\n")

        # Print header for the track
        f.writelines("#  \n")

        f.writelines(
            "#   Date   Time      Lat      Lon     Vmax       Pc     Rmax  R35(NE)  R35(SE)  R35(SW)  R35(NW)  R50(NE)  R50(SE)  R50(SW)  R50(NW)  R65(NE)  R65(SE)  R65(SW)  R65(NW) R100(NE) R100(SE) R100(SW) R100(NE)\n"
        )
        f.writelines("#  \n")

        # Print the actual track
        for i in range(len(gdf)):

            f.writelines(gdf.datetime[i].rjust(15))
            coords = gdf.geometry[i]
            f.writelines(str(round(coords.y, 2)).rjust(9))
            f.writelines(str(round(coords.x, 2)).rjust(9))

            f.writelines(str(round(gdf.vmax[i], 1)).rjust(9))
            f.writelines(str(round(gdf.pc[i], 1)).rjust(9))
            f.writelines(str(round(gdf.rmw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r35_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r35_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r50_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r50_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r65_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r65_nw[i], 1)).rjust(9))

            f.writelines(str(round(gdf.r100_ne[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_se[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_sw[i], 1)).rjust(9))
            f.writelines(str(round(gdf.r100_nw[i], 1)).rjust(9))

            f.writelines("\n")

def interpolate_track(gdf0, gdf1, method="linear"):
    """Interpolate track0 to track1. Only the times need to be provided in track1. The rest will be interpolated. Method is either 'linear' or 'spline'"""

    track0 = gdf_to_dict(gdf0)
    track1 = gdf_to_dict(gdf1)

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
        gdf_point = GeoDataFrame(geometry=[point])
        # gdf_point["geometry"] = [Point(track1["x"][index], track1["y"][index])]  # Assign the new geometry
        gdf_point["datetime"] = track1["datetime"][index].strftime(dateformat_module)
        # Loop through items in track1
        for key, value in track1.items():
            if key == "datetime" or key == "timestamp" or key == "x" or key == "y":
                continue
            else:
                gdf_point[key] = track1[key][index]
        gdf = pd.concat([gdf, gdf_point], ignore_index=True)        
    
    gdf = gdf.set_crs(gdf0.crs, inplace=True)

    return gdf

def interpolate_to_point(gdf, time: datetime, method="spline"):
    """Return interpolated track point at a certain time"""            
    # Find index where to insert
    # Check first if time is before the first time
    t0 = datetime.strptime(gdf.datetime[0], dateformat_module)
    if time < t0:
        index = 0
    else:
        for it, row in gdf.iterrows():
            if time == datetime.strptime(row.datetime, dateformat_module):
                raise Exception("Error: time already exists in track")
            if time > datetime.strptime(row.datetime, dateformat_module):
                index = it + 1

    if index == 0:
        ipos = 0
        vfac = -1.0
    elif index == len(gdf):       
        ipos = len(gdf) - 1
        vfac = 1.0

    if index == 0 or index == len(gdf):
        # Need to estimate position from forward speed and direction

        # Take starting point of original track
        gdf_point = gdf.iloc[[ipos]].copy()

        # Get time difference
        dt = np.abs((time - datetime.strptime(gdf.datetime[ipos], dateformat_module)).total_seconds())

        # Get forward speed at the start
        vtx = vfac * gdf.vtx[ipos]
        vty = vfac * gdf.vty[ipos]

        # Starting point
        lon0 = gdf.geometry[ipos].x
        lat0 = gdf.geometry[ipos].y

        # Determine azimuth and direction
        az = 90.0 - np.arctan2(vty, vtx) * 180 / np.pi
        dst = np.sqrt(vtx ** 2 + vty ** 2) * dt

        # Along track error shift
        lon1, lat1, backaz = geodesic.fwd(
            lon0,
            lat0,
            az,
            dst,
            radians=False,
        )
        point = Point(lon1, lat1)
        gdf_point["geometry"] = [point]
        gdf_point["datetime"] = [time.strftime(dateformat_module)]

    else:
        # Need to do full interpolation
        gdf_point = gdf.iloc[[0]].copy()
        gdf_point["datetime"] = [time.strftime(dateformat_module)]
        gdf_point = interpolate_track(gdf, gdf_point, method=method)

    return gdf_point, index    

def gdf_to_dict(gdf):

    track0 = {}

    track0["datetime"] = []
    track0["x"]        = []
    track0["y"]        = []
    # Loop through columns in gdf
    for column in gdf.items():
        if column[0] == "datetime":
            continue
        elif column[0] == "geometry":
            continue
        else:
            track0[column[0]] = []

    # Now loop through the times
    for index, row in gdf.iterrows():
        # Convert to datetime
        track0["datetime"].append(datetime.strptime(row.datetime, dateformat_module))
        # Get x and y
        track0["x"].append(row.geometry.x)
        track0["y"].append(row.geometry.y)
        # Loop through columns in gdf
        for column in gdf.items():
            if column[0] == "datetime":
                continue
            elif column[0] == "geometry":
                continue
            else:
                track0[column[0]].append(row[column[0]])

    # Convert to timestamp for interpolation
    track0["timestamp"] = [date.timestamp() for date in track0["datetime"]]

    return track0
