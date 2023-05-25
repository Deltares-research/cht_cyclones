import datetime
import os
from functools import reduce

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr

from .tropical_cyclone import TropicalCyclone


class CycloneTrackDatabase:
    def __init__(self, name, file_name=None):
        if name == "ibtracs":
            self.name = name
            self.ds = xr.open_dataset(file_name)
            # Read in database
            self.lon = self.ds["lon"].values[:]
            self.lat = self.ds["lat"].values[:]
            self.basin = self.ds["basin"].values[:, 0].astype(str).tolist()
            self.name = self.ds["name"].values[:].astype(str).tolist()
            self.year = self.ds["season"].values[:].astype(int)
            self.nstorms = np.shape(self.lon)[0]
            self.ntimes = np.shape(self.lon)[1]

        else:
            pass

    def get_track(self, index):
        tc = TropicalCyclone(name=self.name[index])
        for it in range(self.ntimes):
            if not np.isfinite(self.lon[index, it]):
                break
            point = shapely.Point(self.lon[index, it], self.lat[index, it])
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
            R100_NE = -999.0
            R100_SE = -999.0
            R100_SW = -999.0
            R100_NW = -999.0
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
            gdf.set_crs(epsg=4326, inplace=True)

            # Append self
            tc.track = pd.concat([tc.track, gdf])

        # Done with this
        tc.track = tc.track.reset_index(drop=True)
        tc.track = tc.track.drop([0])  # remove the dummy
        tc.track = tc.track.reset_index(drop=True)

        return tc

    def to_gdf(self, index=None):
        geom = []
        iok = []
        if index is None:
            index = range(self.nstorms)
        description = []
        # Returns a gdf with tracks
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
    ):
        if year:
            year_min = year
            year_max = year
        else:
            if not year_min:
                year_min = 0
            if not year_max:
                year_max = 9999

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

        if year_min and year_max:
            iyear = np.where((self.year >= year_min) & (self.year <= year_max))[0]
        else:
            iyear = np.arange(0, self.nstorms)

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

        if distance:
            # Compute distance of all tracks to point
            d = compute_distance(lon, lat, self.lon, self.lat)
            dmin = np.nanmin(d, axis=1)
            idist = np.where(dmin < distance)[0]

        else:
            idist = np.arange(0, self.nstorms)

        index = reduce(np.intersect1d, (ibasin, iyear, iname, idist, ibbox))

        return index

    def track_selector(
        self, app, lon=0.0, lat=0.0, distance=1000.0, year_min=1850, year_max=2030
    ):
        app.gui.setvar("cyclone_track_selector", "lon", lon)
        app.gui.setvar("cyclone_track_selector", "lat", lat)
        app.gui.setvar("cyclone_track_selector", "distance", distance)
        app.gui.setvar("cyclone_track_selector", "year0", year_min)
        app.gui.setvar("cyclone_track_selector", "year1", year_max)
        app.gui.setvar("cyclone_track_selector", "name", "")

        data = {}
        data["track_database"] = self

        # Read GUI config file
        config = app.gui.read_gui_config(
            os.path.dirname(__file__), "cyclone_track_selector.yml"
        )
        okay, data = app.gui.popup(config, data)

        track = None
        if okay:
            # Get the track from the database
            track = self.get_track(data["database_index"])

        return track, okay


def compute_distance(lon1, lat1, lon2, lat2):
    R = 6373.0
    lon1 = lon1 * np.pi / 180
    lat1 = lat1 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    lat2 = lat2 * np.pi / 180
    dlon = lon2 - lon1
    dlon[np.where(dlon < -np.pi)] = dlon[np.where(dlon < -np.pi)] + 2 * np.pi
    dlon[np.where(dlon > np.pi)] = dlon[np.where(dlon > np.pi)] - 2 * np.pi
    dlat = lat2 - lat1
    a = (np.sin(dlat / 2)) ** 2 + np.cos(lat1) * np.cos(lat2) * (np.sin(dlon / 2)) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))
    return R * c


### Callback functions for track selector ###


def map_ready():
    try:
        from app import app
    except Exception:
        # This should be removed. ddb in Delft Dashboard should be renamed to app!
        from ddb import ddb as app

    mp = app.gui.popup_window.find_element_by_id("track_selector_map").widget
    mp.jump_to(0.0, 0.0, 1)
    data = app.gui.popup_data
    # Container layers
    data["main_layer"] = mp.add_layer("track_selector")
    # Tracks layers
    data["track_layer"] = data["main_layer"].add_layer(
        "tracks",
        type="line_selector",
        file_name="tracks.geojson",
        select=select_track,
        selection_type="single",
        line_color="dodgerblue",
        line_width=2,
        line_color_selected="red",
        line_width_selected=3,
        hover_param="description",
    )

    # Update data in tracks layer
    update_tracks()


def update_tracks():
    try:
        from app import app
    except Exception:
        # This should be removed. ddb in Delft Dashboard should be renamed to app!
        from ddb import ddb as app

    data = app.gui.popup_data

    tdb = data["track_database"]
    tracks_layer = data["track_layer"]

    # Get filter data
    distance = app.gui.getvar("cyclone_track_selector", "distance")
    year_min = app.gui.getvar("cyclone_track_selector", "year0")
    year_max = app.gui.getvar("cyclone_track_selector", "year1")
    lon = app.gui.getvar("cyclone_track_selector", "lon")
    lat = app.gui.getvar("cyclone_track_selector", "lat")

    # Get indices based on filter
    index = tdb.filter(
        lon=lon, lat=lat, distance=distance, year_min=year_min, year_max=year_max
    )

    # Get GeoDataFrame of tracks
    gdf = tdb.to_gdf(index=index)

    tracks_layer.set_data(gdf, 0)


def map_moved(coords):
    pass


def select_track(feature):
    try:
        from app import app
    except Exception:
        # This should be removed. ddb in Delft Dashboard should be renamed to app!
        from ddb import ddb as app

    app.gui.popup_data["database_index"] = feature["properties"]["database_index"]


def edit_filter(*args):
    update_tracks()
