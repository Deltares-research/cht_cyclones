import datetime

import numpy as np
import pandas as pd
import xarray as xr
from geopandas import GeoDataFrame
from pyproj import CRS
from shapely.geometry import LineString, MultiLineString, Point, mapping

knots_to_ms = float(0.51444)
nm_to_km = float(1.852)
nm_to_m = float(1.852) * 1000
dateformat_module = "%Y%m%d %H%M%S"

from .track import TropicalCycloneTrack


class TropicalCycloneSpiderweb:
    def __init__(self):
        self.ds = xr.Dataset()
        self.ds.attrs["description"] = "Tropical cyclone spiderweb by cht_cyclones"

    def read(self, filename):
        # Read spw file
        # Get file extension
        fmt = filename.split(".")[-1]        
        if fmt == "nc":
            self.ds = xr.open_dataset(filename)
        elif fmt == "spw":
            self.read_spiderweb_ascii(filename)

    def initialize_grid(self, track, nrad, ndir, spiderweb_radius):
        # Make Xarray dataset

        nt = len(track.gdf)

        # Initialize arrays
        t = np.zeros(nt, dtype="datetime64[s]")
        lon = np.zeros(nt)
        lat = np.zeros(nt)

        for it in range(nt):
            # convert date string with YYYYMMDD HHMMSS to np.datetime64 object
            t[it] = datetime.datetime.strptime(track.gdf.datetime[it], dateformat_module)
            lon[it] = track.gdf["geometry"][it].coords[0][0]
            lat[it] = track.gdf["geometry"][it].coords[0][1]

        dx = spiderweb_radius / nrad
        r = np.arange(dx, spiderweb_radius + dx, dx)
        dphi = 360 / ndir
        phi = np.arange(90, -270, -dphi) # Phi is ordered in nautical convention

        # Compute the grid
        lon_grid = np.zeros((len(t), len(r), len(phi)))
        lat_grid = np.zeros((len(t), len(r), len(phi)))
        for it in range(len(t)):
            # Create circular grid
            x0 = lon[it]
            y0 = lat[it]
            pp, rr = np.meshgrid(phi, r)
            xx = rr * np.cos(pp * np.pi / 180)
            yy = rr * np.sin(pp * np.pi / 180)
            # Convert to degrees
            lon_grid[it, :, :] = x0 + xx / (111320 * np.cos(y0 * np.pi / 180)) 
            lat_grid[it, :, :] = y0 + yy / 111320
        
        # Time
        da = xr.DataArray(t, coords=[("time", t)])
        # self.ds["time"] = xr.DataArray(t, coords=[("time", t)])
        da.assign_attrs(standard_name="time")
        self.ds["time"] = da

        # Longitude of eye
        da = xr.DataArray(lon, coords=[("time", t)])
        da.attrs["standard_name"] = "longitude_of_eye"
        da.attrs["long_name"] = "longitude of eye"
        da.attrs["units"] = "degrees_east"
        self.ds["longitude_eye"] = da

        # Latitude of eye
        da = xr.DataArray(lat, coords=[("time", t)])
        da.attrs["standard_name"] = "latitude_of_eye"
        da.attrs["long_name"] = "latitude of eye"
        da.attrs["units"] = "degrees_north"
        self.ds["latitude_eye"] = da

        # Longitude of grid
        da = xr.DataArray(lon_grid, coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "longitude"
        da.attrs["long_name"] = "longitude"
        da.attrs["units"] = "degrees_east"
        self.ds["lon"] = da

        # Latitude of grid
        da = xr.DataArray(lat_grid, coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "latitude"
        da.attrs["long_name"] = "latitude"
        da.attrs["units"] = "degrees_north"
        self.ds["lat"] = da

        # And now for the data (set to zeros)

        # Wind U
        da = xr.DataArray(np.zeros((nt, nrad, ndir)), coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "eastward_wind"
        da.attrs["units"] = "m/s"
        self.ds["wind_x"] = da

        # Wind V
        da = xr.DataArray(np.zeros((nt, nrad, ndir)), coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "northward_wind"
        da.attrs["units"] = "m/s"
        self.ds["wind_y"] = da

        # Pressure
        da = xr.DataArray(np.zeros((nt, nrad, ndir)), coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "surface_air_pressure"
        da.attrs["long_name"] = "pressure at the bottom of the atmosphere"
        da.attrs["units"] = "Pa"
        self.ds["pressure"] = da

        da = xr.DataArray(np.zeros((nt, nrad, ndir)), coords=[("time", t), ("range", r), ("azimuth", phi)])
        da.attrs["standard_name"] = "precipitation"
        da.attrs["units"] = "mm/h"
        self.ds["precipitation"] = da


    def write(self, filename, format="netcdf", background_pressure=1013.0, tref=None, include_rainfall=False, merge_frac=0.5):
        if format == "netcdf":
            if not include_rainfall:
                # Drop precipitation if not needed
                self.ds = self.ds.drop_vars("precipitation").to_netcdf(filename)
            else:
                self.ds.to_netcdf(filename)    
        elif format == "ascii":
            self.write_spiderweb_ascii(filename,
                                       background_pressure=background_pressure,
                                       tref=tref,
                                       merge_frac=merge_frac,
                                       include_rainfall=include_rainfall)

    def read_spiderweb_ascii(self, filename):
        # Read in ASCII

        # Open file
        fid = open(filename, "r")

        # First loop through all lines and count number of times that a line starts with "TIME"
        n_times = 0
        while True:
            line = fid.readline()
            if not line:
                # End of file
                break
            if line.split()[0] == "n_rows":
                n_rows = int(line.split()[2])
            if line.split()[0] == "n_cols":
                n_cols = int(line.split()[2])
            if line.split()[0] == "n_quantity":
                n_quantity = int(line.split()[2])    
            if line.split()[0] == "TIME":
                n_times += 1
        # Rewind file
        fid.seek(0)

        # Initialize arrays
        t = np.zeros(n_times, dtype="datetime64[s]")
        lon = np.zeros(n_times)
        lat = np.zeros(n_times)

        # Now find the number of rows and columns
        n_rows = 0
        n_cols = 0
        for i in range(10):
            line = fid.readline()
            if line.split()[0] == "n_rows":
                n_rows = int(line.split()[2])
            if line.split()[0] == "n_cols":
                n_cols = int(line.split()[2])
        # Rewind file
        fid.seek(0)


        # Read header information
        vsn = fid.readline().split()[2]
        filetype = fid.readline().split()[2]
        nodata_value = float(fid.readline().split()[2])
        n_cols = int(fid.readline().split()[2])
        n_rows = int(fid.readline().split()[2])
        grid_unit = fid.readline().split()[2]
        spiderweb_radius = float(fid.readline().split()[2])
        spw_rad_unit = fid.readline().split()[2]
        spw_merge_frac = float(fid.readline().split()[2])
        n_quantity = int(fid.readline().split()[2])
        quantity1 = fid.readline().split()[2]
        quantity2 = fid.readline().split()[2]
        quantity3 = fid.readline().split()[2]
        unit1 = fid.readline().split()[2]
        unit2 = fid.readline().split()[2]
        unit3 = fid.readline().split()[2]

        if n_quantity == 4:
            quantity4 = fid.readline().split()[2]
            unit4 = fid.readline().split()[2]

        # Initialize arrays
        wind_speed = np.zeros((n_times, n_rows, n_cols))
        wind_from_direction = np.zeros((n_times, n_rows, n_cols))
        pressure_drop = np.zeros((n_times, n_rows, n_cols))
        if n_quantity == 4:
            precipitation = np.zeros((n_times, n_rows, n_cols))

        for it in range(n_times):
            line = fid.readline()
            parts = line.split("=")
            # Get the second part
            line = parts[1]
            t0 = float(line.split()[0])
            tunits = line.split()[1]
            trefstrs = line.split()[3:-1]
            trefstr = " ".join(trefstrs)
            tref = np.datetime64(trefstr)
            if tunits[0:2] == "mi":
                t = tref + np.timedelta64(int(t0), "m")
            elif tunits[0:2] == "ho":
                t = tref + np.timedelta64(int(t0), "h")
            elif tunits[0:2] == "se":
                t = tref + np.timedelta64(int(t0), "s")
            x_spw_eye = float(fid.readline().split()[2])
            y_spw_eye = float(fid.readline().split()[2])
            p_drop_spw_eye = float(fid.readline().split()[2])


            # Read the data
            for i in range(n_rows):
                wind_speed[it, i, :] = np.array(fid.readline().split(), dtype=float)
            for i in range(n_rows):
                wind_from_direction[it, i, :] = np.array(fid.readline().split(), dtype=float)
            for i in range(n_rows):
                pressure_drop[it, i, :] = np.array(fid.readline().split(), dtype=float)
            if n_quantity == 4:
                for i in range(n_rows):
                    precipitation[it, i, :] = np.array(fid.readline().split(), dtype=float)

        # Close file
        fid.close()

        # # Create Xarray dataset
        # self.ds = xr.Dataset(
        #     {
        #         "wind_speed": (["time", "range", "azimuth"], wind_speed),
        #         "wind_from_direction": (["time", "range", "azimuth"], wind_from_direction),
        #         "pressure": (["time", "range", "azimuth"], 101200.0 - pressure_drop),
        #     },
        #     coords={
        #         "time": np.arange(0, n_times),
        #         "range": np.arange(0, n_rows),
        #         "azimuth": np.arange(0, n_cols),
        #     },
        # )
        # xxx=1

    def write_spiderweb_ascii(self, filename, background_pressure=1013.0, tref=None, merge_frac=0.5, include_rainfall=False):
        # Write in ASCII

        # Make sure tref is in np.datetime64 format
        if tref is not None:
            # tref comes in as a string, so convert to datetime object
            tref = datetime.datetime.strptime(tref, dateformat_module) 
            tref = np.datetime64(tref)
        else:
            # Round tref to January 1st of the year of tref
            tref = np.datetime64(np.datetime64(self.ds["time"].values[0], "Y"), "D")

        # Header information
        vsn = "1.03"
        gridunit = "degree"

        # The rows and columns need to be switched for python
        nt    = np.shape(self.ds["wind_x"])[0]
        nrows = np.shape(self.ds["wind_x"])[1]
        ncols = np.shape(self.ds["wind_x"])[2]

        spiderweb_radius = self.ds["range"].values[-1]

        # Create output
        fid = open(filename, "w")
        format1 = "{:<14} {:1} {:<} {:>} \n"
        format2 = "{:<14} {:1} {:<} \n"

        # Write header information
        fid.write(
            format1.format(
                "FileVersion",
                "=",
                vsn,
                "                            # Version of meteo input file, to check if the newest file format is used",
            )
        )
        fid.write(
            format2.format(
                "filetype",
                "=",
                "meteo_on_spiderweb_grid          # from TRACK file: trackfile.trk",
            )
        )
        fid.write(format2.format("NODATA_value", "=", "-999.000"))
        fid.write(
            format1.format(
                "n_cols",
                "=",
                str(ncols),
                "                   # Number of columns used for wind datafield",
            )
        )
        fid.write(
            format1.format(
                "n_rows",
                "=",
                str(nrows),
                "                           # Number of rows used for wind datafield",
            )
        )
        fid.write(format2.format("grid_unit", "=", gridunit))
        fid.write(format2.format("spw_radius", "=", str(spiderweb_radius)))
        fid.write(format2.format("spw_rad_unit", "=", "m"))
        if merge_frac:
            fid.write(format2.format("spw_merge_frac", "=", str(merge_frac)))
        # Check is    
        if include_rainfall:
            fid.write(format2.format("n_quantity", "=", "4"))
        else:
            fid.write(format2.format("n_quantity", "=", "3"))
        fid.write(format2.format("quantity1", "=", "wind_speed"))
        fid.write(format2.format("quantity2", "=", "wind_from_direction"))
        fid.write(format2.format("quantity3", "=", "p_drop"))
        if include_rainfall:
            fid.write(format2.format("quantity4", "=", "precipitation"))
        fid.write(format2.format("unit1", "=", "m s-1"))
        fid.write(format2.format("unit2", "=", "degree"))
        fid.write(format2.format("unit3", "=", "Pa"))
        if include_rainfall:
            fid.write(format2.format("unit4", "=", "mm/h"))

        # Go over the time steps
        for it in range(nt):

            # Convert to minutes
            dt = self.ds["time"].values[it] - tref
            dt = dt.astype('timedelta64[m]').astype(int)

            # Get main variables
            wu = self.ds["wind_x"].values[it, :, :]
            wv = self.ds["wind_y"].values[it, :, :]
            wind_speed = np.sqrt(wu**2 + wv**2)
            wind_from_direction = 270.0 - np.arctan2(wv, wu) * 180 / np.pi
            pressure_drop = background_pressure * 100 - self.ds["pressure"].values[it, :, :]
            if include_rainfall:
                rainfall_rate = self.ds["precipitation"].values[it, :, :]

            # Replace NaN with -999
            wind_speed = np.nan_to_num(wind_speed, nan=-999)
            wind_from_direction = np.nan_to_num(wind_from_direction, nan=-999)
            pressure_drop = np.nan_to_num(pressure_drop, nan=-999)
            if include_rainfall:
                rainfall_rate = np.nan_to_num(rainfall_rate, nan=-999)

            # Get coordinates
            lon = self.ds["longitude_eye"].values[it]
            lat = self.ds["latitude_eye"].values[it]

            # Print
            timestamp = ((tref - np.datetime64('1970-01-01T00:00:00'))
                 / np.timedelta64(1, 's'))
            tref_datetime = datetime.datetime.utcfromtimestamp(timestamp)
            reference_time_str = tref_datetime.strftime("%Y-%m-%d %H:%M:%S")
            fid.write(
                "{:<14} {:1} {:^10.2f} {:} {:} {:<6}".format(
                    "TIME", "=", int(dt), "minutes since", reference_time_str, "+00:00"
                )
            )
            fid.write("\n")
            fid.write(format2.format("x_spw_eye", "=", (round(lon, 2))))
            fid.write(format2.format("y_spw_eye", "=", (round(lat, 2))))
            fid.write(
                format2.format(
                    "p_drop_spw_eye", "=", (round(np.amax(pressure_drop), 2))
                )
            )

            # Save main variables
            np.savetxt(fid, wind_speed, fmt="%9.2f")
            np.savetxt(fid, wind_from_direction, fmt="%9.2f")
            np.savetxt(fid, pressure_drop, fmt="%9.2f")
            if include_rainfall == True:
                np.savetxt(fid, rainfall_rate, fmt="%9.2f")

        # We are done here
        fid.close()

    def get_track(self, config):

        track = TropicalCycloneTrack()
        r = self.ds["range"].values

        # Loop through the time steps
        for it in range(len(self.ds["time"])):

            # Get the time
            t = self.ds["time"].values[it]
            # Convert numpy time to datetime
            tc_time_string = np.datetime_as_string(t, unit="s")
            tc_time_string = tc_time_string.replace("T", " ")
            tc_time_string = tc_time_string.replace(":", "")
            tc_time_string = tc_time_string.replace("-", "")

            # Get the eye coordinates
            lon = self.ds["longitude_eye"].values[it]
            lat = self.ds["latitude_eye"].values[it]

            # Obtain values
            # Wind speed
            vx = self.ds["wind_x"].values[it, :, :]
            vy = self.ds["wind_y"].values[it, :, :]
            vmag = np.sqrt(vx**2 + vy**2)
            p = self.ds["pressure"].values[it, :, :]

            # Interpolate to higher resolution radial bins
            rmax = np.amax(r)
            dr = 1000.0
            rnew = np.arange(r[0], rmax + dr, dr)
            # Loop over directional bins
            vmagr = np.zeros((len(rnew), np.shape(vmag)[1]))
            for i in range(np.shape(vmagr)[1]):
                # Interpolate to new radii
                vv = np.interp(rnew, r, vmag[:, i])
                vmagr[:, i] = np.interp(rnew, r, vmag[:, i])

            # Get the maximum wind speed
            vmax = np.amax(vmagr)
            # Get the minimum pressure
            pc = np.amin(p) / 100
            # Get radius of maximum wind
            # r,theta indices where vmax is found            
            idx = np.where(vmag == vmax)
            rmw = 0.001 * r[idx[0][0]] / nm_to_km
            vmax = vmax / knots_to_ms / config["wind_conversion_factor"]

            # And now for the radii
            # Take first column of vmag and append to the end
            vmagkts = np.zeros((np.shape(vmagr)[0], np.shape(vmagr)[1] + 1))
            vmagkts[:, :-1] = vmagr
            vmagkts[:, -1] = vmagr[:, 0]
            # Turn into knots and divide by 0.9
            vmagkts = vmagkts / knots_to_ms / config["wind_conversion_factor"]

            itheta0 = [0, 9, 18, 27]
            itheta1 = [10, 19, 28, 37]
            vv = [35.0, 50.0, 65.0, 100.0]
            rr = np.zeros((4, 4)) - 999.0

            for i in range(4):
                vq = vmagkts[:, itheta0[i]:itheta1[i]]
                for iv in range(len(vv)):
                    # Find the first radius where the wind speed is greater than the threshold
                    idx = np.where(vq > vv[iv])
                    # If there are no values greater than the threshold, continue
                    if len(idx[0]) == 0:
                        continue
                    irmax = np.max(idx[0])
                    rr[i, iv] = 0.001 * rnew[irmax] / nm_to_km

            # Make GeoDataFrame
            point = Point(lon, lat)
            gdf_point = GeoDataFrame(
                {
                    "datetime": [tc_time_string],
                    "geometry": [point],
                    "vmax": [vmax],
                    "pc": [pc],
                    "rmw": [rmw],
                    "r35_ne": [rr[0, 0]],
                    "r35_se": [rr[1, 0]],
                    "r35_sw": [rr[2, 0]],
                    "r35_nw": [rr[3, 0]],
                    "r50_ne": [rr[0, 1]],
                    "r50_se": [rr[1, 1]],
                    "r50_sw": [rr[2, 1]],
                    "r50_nw": [rr[3, 1]],
                    "r65_ne": [rr[0, 2]],
                    "r65_se": [rr[1, 2]],
                    "r65_sw": [rr[2, 2]],
                    "r65_nw": [rr[3, 2]],
                    "r100_ne": [rr[0, 3]],
                    "r100_se": [rr[1, 3]],
                    "r100_sw": [rr[2, 3]],
                    "r100_nw": [rr[3, 3]],
                }
            )

            # Append self
            track.gdf = pd.concat([track.gdf, gdf_point])

        # Replace -999 with NaN
        track.gdf = track.gdf.replace(-999, np.nan)
        track.gdf = track.gdf.reset_index(drop=True)
        track.gdf = track.gdf.set_crs(crs=CRS(4326), inplace=True)

        return track

