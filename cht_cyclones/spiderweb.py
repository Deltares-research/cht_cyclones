import xarray as xr
import numpy as np
import datetime

knots_to_ms = float(0.51444)
nm_to_km = float(1.852)
nm_to_m = float(1.852) * 1000
dateformat_module = "%Y%m%d %H%M%S"

class TropicalCycloneSpiderweb:
    def __init__(self):
        self.ds = xr.Dataset()
        self.ds.attrs["description"] = "Tropical cyclone spiderweb by cht_cyclones"

    def read(self, filename, fmt):
        pass

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
