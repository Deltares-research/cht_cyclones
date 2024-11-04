# -*- coding: utf-8 -*-
"""
Tropical Cyclone Module in the Coastal Hazards Toolkit

Module supports two classes
    TropicalCyclone:             deterministic simulations
    TropicalCycloneEnsemble      probabilistic simulations using the simplified DeMaria et al. (2009) approach

    To do list - 'nice to haves'
        make reading of ddb_cyc not file size related but using actual keywords (since format are changing)
        => also remove spiderweb keywords and make it only related to the actual track (2 files?)
        add more reading formats (e.g. NHC, JTWC, etc.)
        enable coordinate conversions; now it is all WGS 84
        add more rainfall methods + allow for scaling of rainfall in ensembles
"""

# Standard Library Imports
import copy
import os
from datetime import datetime

# Third-Party Library Imports
import geopandas as gpd

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyproj
import toml

from .ensemble import TropicalCycloneEnsemble
from .fit_holland_2010 import fit_wind_field_holland2010
from .spiderweb import TropicalCycloneSpiderweb
from .track import TropicalCycloneTrack
from .wind_profiles import holland2010

geodesic = pyproj.Geod(ellps="WGS84")

# Settings
dateformat_module = "%Y%m%d %H%M%S"
knots_to_ms = float(0.51444)
nm_to_km = float(1.852)
nm_to_m = float(1.852) * 1000
pd.options.mode.chained_assignment = None

# Classes of the Tropical Cyclone module
class TropicalCyclone:
    # Init
    def __init__(self,
                 name="no_name",
                 track_file=None,
                 config_file=None,
                 tau=0):
        """Initialize the Tropical Cyclone object"""

        # Name         
        self.name = name

        # Configuration  
        self.config_file = config_file
        self.set_default_config()
        if config_file is not None:
            self.read_config(config_file)

        # Track    
        self.track = TropicalCycloneTrack()
        self.track_file = track_file
        if track_file is not None:
            # Track file can be a string or a list of strings
            if isinstance(track_file, list):
                self.read_track_files_and_merge(track_file, tau=tau)
            else:
                self.read_track(track_file)

        # Spiderweb wind field
        self.spiderweb = TropicalCycloneSpiderweb()

    def set_default_config(self):
        self.config = {}
        self.config["wind_profile"] = "holland2010"
        self.config["wind_pressure_relation"] = "holland2008"
        self.config["rmw_relation"] = "nederhoff2019"
        self.config["rainfall_relationship"] = "ipet"
        self.config["background_pressure"] = 1012.0
        self.config["wind_conversion_factor"] = 0.93
        self.config["phi_spiral"] = 20.0
        self.config["asymmetry_option"] = "schwerdt1979"
        self.config["phi_trans"] = 45.0
        self.config["include_rainfall"] = True
        self.config["rho_air"] = 1.15  # used in the determination of parametric wind field
        self.config["spiderweb_radius"] = 400.0  # radius in km
        self.config["nr_radial_bins"] = 100
        self.config["nr_directional_bins"] = 36
        self.config["tref"] = "20000101 000000"

    def read_config(self, config_file):
        """Read the configuration file (toml format)""" 
        with open(config_file, "r") as f:
            cfg = toml.load(f)
        self.set_default_config()
        # Now loop over keys in cfg and update self.config
        for key in cfg:
            self.config[key] = cfg[key]

    def write_config(self, config_file):
        """Write the configuration file (toml format)"""        
        # Get path of the file
        path = os.path.dirname(config_file)
        if path and not os.path.exists(path):
            os.makedirs(path)
        with open(config_file, "w") as file:            
            toml.dump(self.config, file)

    def read_track(self, filename):
        # Read track file (automatically determine format)
        config = self.track.read(filename)
        if config is not None:
            # Old ddb_cyc format with config in it
            self.set_default_config()
            # Now loop over keys in cfg and update self.config
            for key in config:
                self.config[key] = config[key]

    def read_track_files_and_merge(self, track_file_list, tau=0):
        # Read track files and merge them
        self.track = TropicalCycloneTrack()
        for itrack, track_file in enumerate(track_file_list):
            track1 = TropicalCycloneTrack()
            if itrack == 0:
                # Do not use tau for the first track
                tau1 = 0
            else:
                tau1 = tau
            track1.read(track_file, tau=tau1)
            self.track.add(track1)

    # def resample_track(self, dt, method="spline"):
    #     """Resample the track to a new time step"""
    #     self.track.resample(dt, method=method)        

    def write_track(self, filename, **kwargs):
        # Write track file
        kwargs["name"] = self.name
        self.track.write(filename, **kwargs)

    def write_spiderweb(self, filename, **kwargs):
        # Write the spiderweb wind field to file
        self.spiderweb.write(filename, **kwargs) 

    def extend_track(self, **kwargs):
        """Extend the track"""
        self.track.extend(**kwargs)

    def shorten_track(self, **kwargs):
        """Shorten the track"""
        self.track.shorten(**kwargs)

    def add_track_point(self, datetime, **kwargs):
        """Add a track point"""
        self.track.add_point(datetime, **kwargs)    

    def compute_metric_track(self):
        """Make a copy of the track and convert units to metric"""
        self.track_metric = copy.deepcopy(self.track)
        self.track_metric.convert_units_imperial_to_metric(self.config)

    def compute_wind_field(self, filename=None, progress_bar=None, format="ascii", error_stats=False):
        """Compute the wind field using the specified wind profile"""

        # 1. First make a copy of self.track. Units in this copy will be converted to metric. Also missing values will be estimated.
        # 2. Convert units to km and m/s
        self.compute_metric_track()

        # 3. Compute forward speed
        self.track_metric.compute_forward_speed()

        # 4. Estimate missing values
        self.track_metric.estimate_missing_values(self.config) 

        # # 3. cut off track points with low wind speeds at beginning and end of track + extent track
        # self.cut_off_low_wind_speeds()
        # self.adjust_track_extent()

        # # 4. account for forward speed (computes several derivative values)
        # self.account_for_forward_speed()

        # 5. Initialize spiderweb grid
        self.spiderweb.initialize_grid(self.track,
                                       self.config["nr_radial_bins"],
                                       self.config["nr_directional_bins"],
                                       self.config["spiderweb_radius"] * 1000.0) # Creates the spiderweb arrays

        r = self.spiderweb.ds["range"].values / 1000 # in kilometres
        phi = self.spiderweb.ds["azimuth"].values

        # 6. Go over time steps in the track
        for it in range(len(self.track_metric.gdf)):
            print(it)

            # Progress bar
            if progress_bar:
                progress_bar.set_value(it)
                if progress_bar.was_canceled():
                    return

            # Copy radii data from vectors to matrix
            quadrants_speed = np.array([35.0, 50.0, 65.0, 100.0]) * knots_to_ms * self.config["wind_conversion_factor"]
            quadrants_radii = np.zeros((4, 4))
            quadrants_radii[0, 0] = self.track_metric.gdf.r35_ne[it]
            quadrants_radii[0, 1] = self.track_metric.gdf.r35_se[it]
            quadrants_radii[0, 2] = self.track_metric.gdf.r35_sw[it]
            quadrants_radii[0, 3] = self.track_metric.gdf.r35_nw[it]
            quadrants_radii[1, 0] = self.track_metric.gdf.r50_ne[it]
            quadrants_radii[1, 1] = self.track_metric.gdf.r50_se[it]
            quadrants_radii[1, 2] = self.track_metric.gdf.r50_sw[it]
            quadrants_radii[1, 3] = self.track_metric.gdf.r50_nw[it]
            quadrants_radii[2, 0] = self.track_metric.gdf.r65_ne[it]
            quadrants_radii[2, 1] = self.track_metric.gdf.r65_se[it]
            quadrants_radii[2, 2] = self.track_metric.gdf.r65_sw[it]
            quadrants_radii[2, 3] = self.track_metric.gdf.r65_nw[it]
            quadrants_radii[3, 0] = self.track_metric.gdf.r100_ne[it]
            quadrants_radii[3, 1] = self.track_metric.gdf.r100_se[it]
            quadrants_radii[3, 2] = self.track_metric.gdf.r100_sw[it]
            quadrants_radii[3, 3] = self.track_metric.gdf.r100_nw[it]

            # Get values ready
            dp = self.config["background_pressure"] - self.track_metric.gdf.pc[it]
            vmax = self.track_metric.gdf.vmax[it]
            pc = self.track_metric.gdf.pc[it]
            rmax = self.track_metric.gdf.rmw[it]
            pn = self.config["background_pressure"]
            coords = self.track_metric.gdf.geometry[it]
            latitude = coords.y

            # Get derivative values
            vt = np.sqrt(
                self.track_metric.gdf.vtx[it] ** 2 + self.track_metric.gdf.vty[it] ** 2
            )  # forward speed - magnitude
            phit = (
                np.arctan2(self.track_metric.gdf.vty[it], self.track_metric.gdf.vtx[it]) * 180 / np.pi
            )  # angle
            dpcdt = self.track_metric.gdf.dpcdt[it]  # pressure gradient in time

            # initialize temp arrays
            wind_speed = np.zeros((len(r), len(phi)))
            wind_to_direction_cart = np.zeros((len(r), len(phi)))
            wind_x = np.zeros((len(r), len(phi)))
            wind_y = np.zeros((len(r), len(phi)))
            air_pressure = np.zeros((len(r), len(phi)))
            precipitation = np.zeros((len(r), len(phi)))

            # We first determine vtcor and phicor. xn is set to 0.5. These may be overwritten if we have observations. By default we set phicor to 45 degrees.
            phia = self.config["phi_trans"] * np.pi / 180
            xn   = 0.5
            # phi = self.config["phicor"] * np.pi / 180
            if self.config["asymmetry_option"] == "schwerdt1979":
                # Schwerdt (1979)
                a = vt / knots_to_ms # convert back to kts
                vtcor = 1.5 * a**0.63 * knots_to_ms
            elif self.config["asymmetry_option"] == "mvo":
                vtcor = vt * 0.6
            elif self.config["asymmetry_option"] == "none":
                vtcor = 0.0
            else:
                raise Exception("This asymmetry_option is not supported")

            # Now check if fitting of the wind field is necessary
            # Holland et al. 2010
            if self.config["wind_profile"] == "holland2010":

                # Count the number of observations, so that we know if we need to do a fitting of with Holland (2010)
                fit_wind = False
                n = 0
                if vmax > 20:
                    for iq in range(np.size(quadrants_radii, 0)):
                        for irad in range(np.size(quadrants_radii, 1)):
                            if not np.isnan(quadrants_radii[iq, irad]):
                                n = n + 1
                if n > 0:
                    fit_wind = True

                # Do fitting of Holland 2010
                if fit_wind:
                    obs = {}
                    obs["quadrants_radii"] = quadrants_radii
                    obs["quadrants_speed"] = quadrants_speed
                    [xn, vtcor, phia] = fit_wind_field_holland2010(
                        vmax,
                        rmax,
                        pc,
                        vt,
                        phit,
                        pn,
                        self.config["phi_spiral"],
                        latitude,
                        dpcdt,
                        obs,
                    )

            if self.config["wind_profile"] == "holland2008":
                # Assume constant xn that follows a relationship described in 2008 paper
                xn = 0.6 * (1 - dp / 215)

            # Now compute the wind field (in goes vmax and vtcor). vr does not include asymmetry!
            if self.config["wind_profile"] == "holland1980" or self.config["wind_profile"] == "holland2010":
                vrel = vmax - vtcor
                [vr, pr] = holland2010(r, vrel, pc, pn, rmax, dpcdt, latitude, vt, xn)
            else:
                raise Exception("This wind_profile is not supported")

            # Loop over the directions
            for iphi in range(len(phi)):

                # Place wind speed
                wind_speed[:, iphi] = vr

                # Check which hemisphere we are
                if latitude >= 0:
                    # northern hemisphere
                    dr = 90.0 + phi[iphi] + self.config["phi_spiral"]
                else:
                    # southern hemisphere
                    dr = -90.0 + phi[iphi] - self.config["phi_spiral"]

                # Wind direction and pressure drop
                wind_to_direction_cart[:, iphi] = dr
                air_pressure[:, iphi] = 100 * pr

            # Wind speed with asymmetry
            # vnorm = wind_speed / np.max(wind_speed)
            vnorm = np.zeros(np.shape(wind_speed)) + 1.0
            u_prop = vtcor * np.cos((phit + phia) * np.pi / 180)
            v_prop = vtcor * np.sin((phit + phia) * np.pi / 180)
            wind_x = wind_speed * np.cos(wind_to_direction_cart * np.pi / 180) + u_prop * vnorm
            wind_y = wind_speed * np.sin(wind_to_direction_cart * np.pi / 180) + v_prop * vnorm

            # Rainfall (if required)
            if self.config["include_rainfall"] == True:
                # Add rainfall
                if self.config["rainfall_relationship"] == "ipet":
                    # IPET is a simple rainfall model relating pressure to rainfall rate
                    pdef = (self.config["background_pressure"] * 100 - air_pressure[0, 0]) / 100  # % hPa to Pa
                    pdef = max(pdef, 0.0)
                    for ip in range(len(r)):
                        if r[ip] < rmax:
                            precipitation[ip, :] = 1.14 + (0.12 * pdef)
                        else:
                            precipitation[ip, :] = (1.14 + (0.12 * pdef)) * np.exp(
                                -0.3 * ((r[ip] - rmax) / rmax)
                            )
                # More options to be added later
                # Bader, Bacla, etc.

            # Update data in spiderweb array
            self.spiderweb.ds["wind_x"].values[it, :, :] = wind_x
            self.spiderweb.ds["wind_y"].values[it, :, :] = wind_y
            self.spiderweb.ds["pressure"].values[it, :, :] = air_pressure
            if self.config["include_rainfall"]:
                self.spiderweb.ds["precipitation"].values[it, :, :] = precipitation


            # # Scale wind speeds back to ensure we always reach vmax
            # missing_factor = vmax / np.max(wind_speed)
            # wind_speed = wind_speed * missing_factor

        # Default is spiderweb in ascii
        if filename is not None:
            self.spiderweb.write(filename,
                                 format=format,
                                 background_pressure=self.config["background_pressure"],
                                 tref=self.config["tref"],
                                 include_rainfall=self.config["include_rainfall"])

        # If error stats are required
        if error_stats:
            # Compute the error stats
            error_stats = self.compute_track_error_stats()


    def get_track_from_spiderweb(self, replace=False):
        """Get the track from the spiderweb."""
        track = self.spiderweb.get_track(self.config)
        if replace:
            self.track = track
        return track    

    def compute_track_error_stats(self):
        """Compute error statistics"""
        resulting_track = self.spiderweb.get_track(self.config)
        r_best   = np.empty((len(resulting_track.gdf), 4, 4)) * np.nan
        r_result = np.empty((len(resulting_track.gdf), 4, 4)) * np.nan
        r_diff   = np.empty((len(resulting_track.gdf), 4, 4)) * np.nan
        r_diff_rel  = np.empty((len(resulting_track.gdf), 4, 4)) * np.nan
        r_diff_radius_mean = np.zeros((len(resulting_track.gdf), 4))
        r_diff_radius_rel_mean = np.zeros((len(resulting_track.gdf), 4))
        r_diff_mean = np.zeros((len(resulting_track.gdf)))
        r_diff_rel_mean = np.zeros((len(resulting_track.gdf)))
        rstr = ["r35", "r50", "r65", "r100"]
        vstr = ["ne", "se", "sw", "nw"]
        # Loop through time
        for it in range(len(resulting_track.gdf)):
            print("Time step : " + str(it))
            for ir in range(4):
                iok = 0
                ibst = 0
                ires = 0
                for iv in range(4):
                    r_best[it, ir, iv]   = self.track.gdf[rstr[ir] + "_" + vstr[iv]][it]
                    r_result[it, ir, iv] = resulting_track.gdf[rstr[ir] + "_" + vstr[iv]][it]
                    # Compute the error statistics for this wind speed
                    r_diff[it, ir, iv] = r_result[it, ir, iv] - r_best[it, ir, iv]
                    r_diff_rel[it, ir, iv] = 100 * (r_result[it, ir, iv] - r_best[it, ir, iv]) / r_best[it, ir, iv]
                    if not np.isnan(r_best[it, ir, iv]) and not np.isnan(r_result[it, ir, iv]):
                        iok += 1
                    if not np.isnan(r_best[it, ir, iv]) and np.isnan(r_result[it, ir, iv]):
                        ibst += 1
                    if np.isnan(r_best[it, ir, iv]) and not np.isnan(r_result[it, ir, iv]):
                        ires += 1    
                r_diff_radius_mean[it, ir] = np.nanmean(r_diff[it, ir, :])
                r_diff_radius_rel_mean[it, ir] = np.nanmean(r_diff_rel[it, ir, :])

#                print("Mean error for " + rstr[ir] + " : " + '{0:.1f}'.format(r_diff_radius_mean[it, ir]) + " NM (" + '{0:.0f}'.format(r_diff_radius_rel_mean[it, ir]) + " %)" + "    - " + str(iok) + " OK, " + str(ibst) + " Best, " + str(ires) + " Result")
                print("Mean error for " + rstr[ir] + " : " + '{0:.1f}'.format(r_diff_radius_mean[it, ir]) + " NM (" + '{0:.0f}'.format(r_diff_radius_rel_mean[it, ir]) + " %)")

            r_diff_mean[it] = np.nanmean(r_diff_radius_mean[it, :])
            r_diff_rel_mean[it] = np.nanmean(r_diff_radius_rel_mean[it, :])

            print("Mean error for all : " + '{0:.1f}'.format(r_diff_mean[it]) + " NM (" + '{0:.0f}'.format(r_diff_rel_mean[it]) + " %)")

    def make_ensemble(self, **kwargs):
        """Make ensemble"""
        ensemble = TropicalCycloneEnsemble(self, **kwargs)        
        ensemble.generate()
        return ensemble
    
    def merge_with_meteo_dataset(self):
        """Merge with meteo dataset. Not implemented yet."""
        pass

    def get_track_from_meteo_dataset(self, meteo_dataset):
        """Get the track from a meteo dataset. Not implemented yet."""
        # Find the track(s)
        pass

    def get_wind_field_from_meteo_dataset(self, meteo_dataset, filename=None, format="ascii"):
        """Get the wind field from a meteo dataset"""

        # Track has already been defined
        self.spiderweb.initialize_grid(self.track,
                                       self.config["nr_radial_bins"],
                                       self.config["nr_directional_bins"],
                                       self.config["spiderweb_radius"] * 1000) # Creates the spiderweb arrays

        # Loop over time steps in the track
        for it in range(len(self.track.gdf)):
            # Get spiderweb grid coordinates
            tt = datetime.strptime(self.track.gdf.datetime[it], dateformat_module)
            xx = self.spiderweb.ds["lon"].values[it, :, :]
            yy = self.spiderweb.ds["lat"].values[it, :, :]
            # Get the data at the specified time (really should start doing this with xarray!)
            self.spiderweb.ds["wind_x"].values[it, :, :] = meteo_dataset.interpolate_variable("wind_u", tt, xx, yy)
            self.spiderweb.ds["wind_y"].values[it, :, :] = meteo_dataset.interpolate_variable("wind_v", tt, xx, yy)
            self.spiderweb.ds["pressure"].values[it, :, :] = meteo_dataset.interpolate_variable("barometric_pressure", tt, xx, yy)
            self.spiderweb.ds["precipitation"].values[it, :, :] = meteo_dataset.interpolate_variable("precipitation", tt, xx, yy)

        # Default is spiderweb in ascii
        if filename is not None:
            self.spiderweb.write(filename,
                                 format=format,
                                 background_pressure=self.config["background_pressure"],
                                 tref=self.config["tref"],
                                 include_rainfall=self.config["include_rainfall"])


    def get_wind_field_from_coamps(self, coamps_path, coamps_cycle, filename=None, format="ascii"):
        """Get the wind field from coamps-tc"""

        import xarray as xr
        from cht_utils.misc_tools import interp2

        # # Track has already been defined
        self.spiderweb.initialize_grid(self.track,
                                       self.config["nr_radial_bins"],
                                       self.config["nr_directional_bins"],
                                       self.config["spiderweb_radius"] * 1000) # Creates the spiderweb arrays

        # 6. Loop over time steps in the track
        for it in range(len(self.track.gdf)):

            # Get time
            tt = datetime.strptime(self.track.gdf.datetime[it], dateformat_module)

            # Get spiderweb grid coordinates
            xx = self.spiderweb.ds["lon"].values[it, :, :]
            yy = self.spiderweb.ds["lat"].values[it, :, :]

            # Convert coamps_cycle to time
            coamps_cycle_time = datetime.strptime(coamps_cycle, "%Y%m%d%H")

            # Get the time difference between tt and coamps_cycle_time in hours
            tau = int(((tt - coamps_cycle_time).total_seconds() / 3600))
            taustr = "tau" + str(tau).zfill(3)

            print("Time step : " + str(it) + " - " + str(tt) + " - " + taustr)

            # First do high-res d03

            coamps_file = os.path.join(coamps_path, "coamps-tc_d03_" + coamps_cycle + "_" + taustr + ".nc")

            if not os.path.exists(coamps_file):
                print("File does not exist : " + coamps_file)
                continue
                                       
            ds = xr.open_dataset(coamps_file)
            # ds = ds.sel(lon=slice(xx.min(), xx.max()), lat=slice(yy.min(), yy.max()))

            lon = ds["lon"].values
            lat = ds["lat"].values
            uu = ds["uuwind"].values
            vv = ds["vvwind"].values
            pp = ds["slpres"].values * 100
            rr = ds["precip"].values

            lon = lon[0, :] - 360.0
            lat = lat[:, 0]

            # Get the data at the specified time (really should start doing this with xarray!)
            self.spiderweb.ds["wind_x"].values[it, :, :] = interp2(lon, lat, uu, xx, yy)
            self.spiderweb.ds["wind_y"].values[it, :, :] = interp2(lon, lat, vv, xx, yy)
            self.spiderweb.ds["pressure"].values[it, :, :] = interp2(lon, lat, pp, xx, yy)
            self.spiderweb.ds["precipitation"].values[it, :, :] = interp2(lon, lat, rr, xx, yy)

            # Now do low-res d02
            # But only where we have missing values

            coamps_file = os.path.join(coamps_path, "coamps-tc_d02_" + coamps_cycle + "_" + taustr + ".nc")

            if not os.path.exists(coamps_file):
                print("File does not exist : " + coamps_file)
                continue

            ds = xr.open_dataset(coamps_file)

            lon = ds["lon"].values
            lat = ds["lat"].values
            uu = ds["uuwind"].values
            vv = ds["vvwind"].values
            pp = ds["slpres"].values * 100
            rr = ds["precip"].values

            lon = lon[0, :] - 360.0
            lat = lat[:, 0]

            # Get the data at the specified time (really should start doing this with xarray!)
            self.spiderweb.ds["wind_x"].values[it, :, :][np.isnan(self.spiderweb.ds["wind_x"].values[it, :, :])] = interp2(lon, lat, uu, xx, yy)[np.isnan(self.spiderweb.ds["wind_x"].values[it, :, :])]
            self.spiderweb.ds["wind_y"].values[it, :, :][np.isnan(self.spiderweb.ds["wind_y"].values[it, :, :])] = interp2(lon, lat, vv, xx, yy)[np.isnan(self.spiderweb.ds["wind_y"].values[it, :, :])]
            self.spiderweb.ds["pressure"].values[it, :, :][np.isnan(self.spiderweb.ds["pressure"].values[it, :, :])] = interp2(lon, lat, pp, xx, yy)[np.isnan(self.spiderweb.ds["pressure"].values[it, :, :])]
            self.spiderweb.ds["precipitation"].values[it, :, :][np.isnan(self.spiderweb.ds["precipitation"].values[it, :, :])] = interp2(lon, lat, rr, xx, yy)[np.isnan(self.spiderweb.ds["precipitation"].values[it, :, :])]

        # Replace nan with 0.0 in wind field
        self.spiderweb.ds["wind_x"].values[np.isnan(self.spiderweb.ds["wind_x"].values)] = 0.0
        self.spiderweb.ds["wind_y"].values[np.isnan(self.spiderweb.ds["wind_y"].values)] = 0.0  
        # Replace nan with background pressure        
        self.spiderweb.ds["pressure"].values[np.isnan(self.spiderweb.ds["pressure"].values)] = self.config["background_pressure"]
        # Replace nan with 0.0 in precipitation
        self.spiderweb.ds["precipitation"].values[np.isnan(self.spiderweb.ds["precipitation"].values)] = 0.0

        # Default is spiderweb in ascii
        if filename is not None:
            self.spiderweb.write(filename,
                                 format=format,
                                 background_pressure=self.config["background_pressure"],
                                 tref=self.config["tref"],
                                 include_rainfall=self.config["include_rainfall"])


    def to_gdf(self, filename=None):
        """Convert track to GeoDataFrame, and optionally write to file."""
        return self.track.to_gdf(filename=filename)


# Definition to find if this is a land point
def analyze_points_with_shapefile(
    shapefile_polygon, shapefile_polyline, latitudes, longitudes
):
    from shapely.geometry import Point
    from shapely.ops import nearest_points

    # Read the shapefile
    shapefile_polygon = gpd.read_file(shapefile_polygon)
    shapefile_polyline = gpd.read_file(shapefile_polyline)

    # Empty results
    results = []

    for latitude, longitude in zip(latitudes, longitudes):
        # Create a Point object from the latitude and longitude
        point = Point(longitude, latitude)

        # Check if the point falls within any polygon of the shapefile
        is_within = shapefile_polygon.contains(point).any()

        if is_within:
            # Find the nearest vertex from the point
            nearest_vertex = nearest_points(point, shapefile_polyline.unary_union)[1]

            # Calculate the distance between the point and the nearest vertex
            distance_km = point.distance(nearest_vertex) * 111.32

            result = {
                "latitude": latitude,
                "longitude": longitude,
                "is_within": True,
                "distance_km": distance_km,
            }
        else:
            result = {
                "latitude": latitude,
                "longitude": longitude,
                "is_within": False,
                "distance_km": None,
            }

        results.append(result)

    return results
