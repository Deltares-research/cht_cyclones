import datetime
import os
from functools import reduce

import boto3
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import toml
import xarray as xr
from botocore import UNSIGNED
from botocore.client import Config
from cht_utils import fileops as fo

from .tropical_cyclone_refactored import TropicalCyclone


class CycloneTrackDataset:
    """
    CycloneTrackDataset is a class that reads in a dataset of tropical cyclone
    tracks and provides methods to filter the dataset and return a geopandas
    dataframe of the tracks.

    Parameters
    ----------
    name : str
        Name of the dataset. Currently only "ibtracs" is supported.
    file_name : str
        File name of the dataset. Currently only "ibtracs" is supported.

    Attributes
    ----------
    name : str
        Name of the dataset. Currently only "ibtracs" is supported.
    ds : xarray.Dataset
        Dataset of the dataset.
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
        Number of storms in the dataset.
    ntimes : int
        Number of time steps in the dataset.
    """

    def __init__(self, name, path):
        self.name = name
        self.path = path
        self.ds = None
        self.lon = None
        self.lat = None
        self.basin = None
        self.storm_name = None
        self.year = None
        self.nstorms = None
        self.ntimes = None

        self.read_metadata()

    def read(self):
        if self.ds is not None:
            # Already read in
            return
        filename = os.path.join(self.path, self.files[0])
        self._read_ibtracs(filename)

    def read_metadata(self):
        # Read metadata file
        tml_file = os.path.join(self.path, "metadata.tml")
        tml = toml.load(tml_file)
        for key in tml:
            setattr(self, key, tml[key])
        # Long name for backwards compatibility
        if "longname" in tml:
            self.long_name = tml["longname"]
        # Make sure there is always a long_name
        if self.long_name == "":
            self.long_name = self.name        

    def download(self):
        if self.s3_bucket is None:
            return
        # Check if download is needed
        for file in self.files:
            if not os.path.exists(os.path.join(self.path, file)):
                s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))
                break
        # Get all files defined in the toml file
        for file in self.files:
            if not os.path.exists(os.path.join(self.path, file)):
                print(f"Downloading {file} from track dataset {self.name} ...")
                s3_client.download_file(self.s3_bucket, f"{self.s3_key}/{file}", os.path.join(self.path, file))

    def _read_ibtracs(self, file_name):
        """
        Read in IBTrACS dataset

        Parameters
        ----------
        file_name : str
            File name of the IBTrACS dataset.

        Returns
        -------
        None
        """

        # Read in dataset
        self.ds = xr.open_dataset(file_name)

        # Convert to numpy arrays
        self.lon = self.ds["lon"].values[:]
        self.lat = self.ds["lat"].values[:]
        self.basin = self.ds["basin"].values[:, 0].astype(str).tolist()
        self.storm_name = self.ds["name"].values[:].astype(str).tolist()
        self.year = self.ds["season"].values[:].astype(int)
        self.nstorms = np.shape(self.lon)[0]
        self.ntimes = np.shape(self.lon)[1]

    def get_track(self, index):
        """
        Get a single track from the dataset.

        Parameters
        ----------
        index : int
            Index of the track in the dataset.

        Returns
        -------
        tc : TropicalCyclone
            TropicalCyclone object of the track.
        """

        # Create a TropicalCyclone object
        tc = TropicalCyclone(name=self.storm_name[index])

        gdf = gpd.GeoDataFrame()

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
            gdf_point = gpd.GeoDataFrame(
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
        gdf = gdf.set_crs(crs=4326, inplace=True)

        tc.track.gdf = gdf

        return tc

    def list_names(self, index=None):
        if index is not None:
            return [self.storm_name[i] for i in index]
        else:
            return self.storm_name

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
                description.append(self.storm_name[ind] + " (" + str(self.year[ind]) + ")")

        # Create geopandas dataframe
        gdf = gpd.GeoDataFrame(crs=4326, geometry=geom)
        names = [self.storm_name[i] for i in iok]
        years = self.year[iok]
        gdf["name"] = names
        gdf["year"] = years
        gdf["description"] = description
        gdf["dataset_index"] = iok
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
        Filter the dataset and return the indices of the filtered tracks.

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
                    for i in range(len(self.storm_name))
                    if self.storm_name[i].lower() == name.lower()
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
            d = CycloneTrackDataset.compute_distance(lon, lat, self.lon, self.lat)
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
