import datetime
from functools import reduce

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr

from .tropical_cyclone import TropicalCyclone


class CycloneTrackDatabase:
    """
    CycloneTrackDatabase is a class that reads in a database of tropical cyclone
    tracks and provides methods to filter the database and return a geopandas
    dataframe of the tracks.

    Parameters
    ----------
    name : str
        Name of the database. Currently only "ibtracs" is supported.
    file_name : str
        File name of the database. Currently only "ibtracs" is supported.

    Attributes
    ----------
    name : str
        Name of the database. Currently only "ibtracs" is supported.
    ds : xarray.Dataset
        Dataset of the database.
    lon : numpy.ndarray
        Longitude of the tracks.
    lat : numpy.ndarray
        Latitude of the tracks.
    basin : list
        Basin of the tracks.
    name : list
        Name of the tracks.
    year : numpy.ndarray
        Year of the tracks.
    nstorms : int
        Number of storms in the database.
    ntimes : int
        Number of time steps in the database.
    """

    def __init__(self, name, file_name=None):
        self.name = name
        self.file_name = file_name
        self.ds = None
        self.lon = None
        self.lat = None
        self.basin = None
        self.name = None
        self.year = None
        self.nstorms = None
        self.ntimes = None

        if name == "ibtracs":
            self._read_ibtracs(file_name)
        else:
            pass

    def _read_ibtracs(self, file_name):
        """
        Read in IBTrACS database

        Parameters
        ----------
        file_name : str
            File name of the IBTrACS database.

        Returns
        -------
        None
        """

        # Read in database
        self.ds = xr.open_dataset(file_name)

        # Convert to numpy arrays
        self.lon = self.ds["lon"].values[:]
        self.lat = self.ds["lat"].values[:]
        self.basin = self.ds["basin"].values[:, 0].astype(str).tolist()
        self.name = self.ds["name"].values[:].astype(str).tolist()
        self.year = self.ds["season"].values[:].astype(int)
        self.nstorms = np.shape(self.lon)[0]
        self.ntimes = np.shape(self.lon)[1]

    def get_track(self, index):
        """
        Get a single track from the database.

        Parameters
        ----------
        index : int
            Index of the track in the database.

        Returns
        -------
        tc : TropicalCyclone
            TropicalCyclone object of the track.
        """

        # Create a TropicalCyclone object
        tc = TropicalCyclone(name=self.name[index])

        # Add track
        for it in range(self.ntimes):
            # Check if track is finite
            if not np.isfinite(self.lon[index, it]):
                break

            # Create a shapely point
            point = shapely.Point(self.lon[index, it], self.lat[index, it])

            # Initialize data for geopandas dataframe
            t = self.ds["time"].values[index, it]
            tc_time_string = datetime.datetime.utcfromtimestamp(
                t.item() / 10**9
            ).strftime("%Y%m%d %H%M%S")
            vmax = self.ds["usa_wind"].values[index, it]
            pc = self.ds["usa_pres"].values[index, it]
            RMW = self.ds["usa_rmw"].values[index, it]
            R35_NE = self.ds["usa_r34"].values[index, it, 0]
            R35_SE = self.ds["usa_r34"].values[index, it, 1]
            R35_SW = self.ds["usa_r34"].values[index, it, 2]
            R35_NW = self.ds["usa_r34"].values[index, it, 3]
            R50_NE = self.ds["usa_r50"].values[index, it, 0]
            R50_SE = self.ds["usa_r50"].values[index, it, 1]
            R50_SW = self.ds["usa_r50"].values[index, it, 2]
            R50_NW = self.ds["usa_r50"].values[index, it, 3]
            R65_NE = self.ds["usa_r64"].values[index, it, 0]
            R65_SE = self.ds["usa_r64"].values[index, it, 1]
            R65_SW = self.ds["usa_r64"].values[index, it, 2]
            R65_NW = self.ds["usa_r64"].values[index, it, 3]
            # R100_NE = self.ds["usa_r100"].values[index, it, 0]
            # R100_SE = self.ds["usa_r100"].values[index, it, 1]
            # R100_SW = self.ds["usa_r100"].values[index, it, 2]
            # R100_NW = self.ds["usa_r100"].values[index, it, 3]
            vmax = np.nan_to_num(vmax, copy=True, nan=-999.0)
            pc = np.nan_to_num(pc, copy=True, nan=-999.0)
            RMW = np.nan_to_num(RMW, copy=True, nan=-999.0)
            R35_NE = np.nan_to_num(R35_NE, copy=True, nan=-999.0)
            R35_SE = np.nan_to_num(R35_SE, copy=True, nan=-999.0)
            R35_SW = np.nan_to_num(R35_SW, copy=True, nan=-999.0)
            R35_NW = np.nan_to_num(R35_NW, copy=True, nan=-999.0)
            R50_NE = np.nan_to_num(R50_NE, copy=True, nan=-999.0)
            R50_SE = np.nan_to_num(R50_SE, copy=True, nan=-999.0)
            R50_SW = np.nan_to_num(R50_SW, copy=True, nan=-999.0)
            R50_NW = np.nan_to_num(R50_NW, copy=True, nan=-999.0)
            R65_NE = np.nan_to_num(R65_NE, copy=True, nan=-999.0)
            R65_SE = np.nan_to_num(R65_SE, copy=True, nan=-999.0)
            R65_SW = np.nan_to_num(R65_SW, copy=True, nan=-999.0)
            R65_NW = np.nan_to_num(R65_NW, copy=True, nan=-999.0)
            R100_NE = -999.0
            R100_SE = -999.0
            R100_SW = -999.0
            R100_NW = -999.0

            # Create a geopandas dataframe
            gdf = gpd.GeoDataFrame(
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

            # Set CRS coordinate system
            gdf.set_crs(epsg=4326, inplace=True)

            # Append self
            tc.track = pd.concat([tc.track, gdf])

        # Done with this
        tc.track = tc.track.reset_index(drop=True)
        tc.track = tc.track.drop([0])  # remove the dummy
        tc.track = tc.track.reset_index(drop=True)

        return tc

    def list_names(self, index=None):
        if index is not None:
            return [self.name[i] for i in index]
        else:
            return self.name

    def to_gdf(self, index=None):
        """
        Returns a geopandas dataframe of the tracks.

        Parameters
        ----------
        index : list
            List of indices of the tracks to return. If None, all tracks are returned.

        Returns
        -------
        gdf : geopandas.GeoDataFrame
            Geopandas dataframe of the tracks.
        """

        # Initialize
        geom = []
        iok = []
        if index is None:
            index = range(self.nstorms)
        description = []

        # Loop over tracks
        for ind in index:
            lon = self.lon[ind, :]
            lat = self.lat[ind, :]
            lon = lon[np.isfinite(lon)]
            lat = lat[np.isfinite(lat)]
            if len(lon) > 1 and len(lat) > 1:
                geom.append(
                    shapely.geometry.LineString(np.transpose(np.stack((lon, lat))))
                )
                iok.append(ind)
                description.append(self.name[ind] + " (" + str(self.year[ind]) + ")")

        # Create geopandas dataframe
        gdf = gpd.GeoDataFrame(crs=4326, geometry=geom)
        names = [self.name[i] for i in iok]
        years = self.year[iok]
        gdf["name"] = names
        gdf["year"] = years
        gdf["description"] = description
        gdf["database_index"] = iok
        return gdf

    def filter(
        self,
        name=None,
        distance=None,
        lon=None,
        lat=None,
        basin=None,
        year=None,
        year_min=None,
        year_max=None,
        vmax_min=None,
        vmax_max=None,
    ):
        """
        Filter the database and return the indices of the filtered tracks.

        Parameters
        ----------
        name : str
            Name of the track to filter by.
        distance : float
            Distance (km) of the track to the point (lon, lat).
        lon : list
            Longitude range of the tracks.
        lat : list
            Latitude range of the tracks.
        basin : str
            Basin of the tracks.
        year : int
            Year of the tracks.
        year_min : int
            Minimum year of the tracks.
        year_max : int
            Maximum year of the tracks.

        Returns
        -------
        index : list
            List of indices of the filtered tracks.
        """

        # Initialize year_min and year_max
        if year:
            year_min = year
            year_max = year
        else:
            if not year_min:
                year_min = 0
            if not year_max:
                year_max = 9999

        if not vmax_min:
            vmax_min = 0.0
        if not vmax_max:
            vmax_max = 9999.0

        # Filter by basin
        if basin:
            ibasin = np.array(
                [
                    i
                    for i in range(len(self.basins))
                    if self.basin[i].lower() == basin.lower()
                ]
            )
        else:
            ibasin = np.arange(0, self.nstorms)

        # Filter by year
        if year_min and year_max:
            iyear = np.where((self.year >= year_min) & (self.year <= year_max))[0]
        else:
            iyear = np.arange(0, self.nstorms)

        # # Filter by vmax
        # if vmax_min and vmax_max:
        #     ivmax = np.where((self.year >= year_min) & (self.year <= year_max))[0]
        # else:
        #     iyear = np.arange(0, self.nstorms)

        # Filter by name
        if name:
            iname = np.array(
                [
                    i
                    for i in range(len(self.name))
                    if self.name[i].lower() == name.lower()
                ]
            )
        else:
            iname = np.arange(0, self.nstorms)

        # Filter by bounding box
        if isinstance(lon, list) and isinstance(lat, list):
            inear = np.where(
                (self.lon >= lon[0])
                & (self.lon <= lon[1])
                & (self.lat >= lat[0])
                & (self.lat <= lat[1])
            )
            ibbox = np.unique(inear[0])
        else:
            ibbox = np.arange(0, self.nstorms)

        # Filter by distance
        if distance:
            # Compute distance of all tracks to point
            d = CycloneTrackDatabase.compute_distance(lon, lat, self.lon, self.lat)
            dmin = np.nanmin(d, axis=1)
            idist = np.where(dmin < distance)[0]
        else:
            idist = np.arange(0, self.nstorms)

        # Intersect all filters
        index = reduce(np.intersect1d, (ibasin, iyear, iname, idist, ibbox))

        return index

    @staticmethod
    def compute_distance(lon1, lat1, lon2, lat2):
        """
        Compute the distance between two points on a sphere.

        Parameters
        ----------
        lon1 : numpy.ndarray
            Longitude of the first point.
        lat1 : numpy.ndarray
            Latitude of the first point.
        lon2 : numpy.ndarray
            Longitude of the second point.
        lat2 : numpy.ndarray
            Latitude of the second point.

        Returns
        -------
        d : numpy.ndarray
            Distance between the points.
        """

        # Convert to radians
        R = 6373.0
        lon1 = lon1 * np.pi / 180
        lat1 = lat1 * np.pi / 180
        lon2 = lon2 * np.pi / 180
        lat2 = lat2 * np.pi / 180
        dlon = lon2 - lon1

        # Wrap around
        dlon[np.where(dlon < -np.pi)] = dlon[np.where(dlon < -np.pi)] + 2 * np.pi
        dlon[np.where(dlon > np.pi)] = dlon[np.where(dlon > np.pi)] - 2 * np.pi
        dlat = lat2 - lat1

        # Compute distance
        a = (np.sin(dlat / 2)) ** 2 + np.cos(lat1) * np.cos(lat2) * (
            np.sin(dlon / 2)
        ) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
        return R * c
