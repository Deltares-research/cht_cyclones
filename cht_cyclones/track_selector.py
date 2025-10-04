import os
from geopandas import GeoDataFrame
from shapely.geometry import Point
from pyproj import CRS, Transformer

# Make gui and map globally available in this module
gui = None
map = None

def track_selector(
    dataset, app, lon=0.0, lat=0.0, distance=1000.0, year_min=1850, year_max=2030, vmax_min=0.0, vmax_max=1000.0
):
    global gui

    gui = app.gui

    app.gui.setvar("cyclone_track_selector", "lon", round(lon, 2))
    app.gui.setvar("cyclone_track_selector", "lat", round(lat, 2))
    gui.setvar("cyclone_track_selector", "distance", distance)
    gui.setvar("cyclone_track_selector", "year0", year_min)
    gui.setvar("cyclone_track_selector", "year1", year_max)
    gui.setvar("cyclone_track_selector", "vmax0", vmax_min)
    gui.setvar("cyclone_track_selector", "vmax1", vmax_max)
    gui.setvar("cyclone_track_selector", "name", "")
    gui.setvar("cyclone_track_selector", "storm_names", [""])
    gui.setvar("cyclone_track_selector", "storm_indices", [0])
    gui.setvar("cyclone_track_selector", "storm_index", 0)

    data = {}
    data["track_dataset"] = dataset

    # Read GUI config file
    config_file = os.path.join(os.path.dirname(__file__), "cyclone_track_selector.yml")
    okay, data = gui.popup(config_file, id="track_selector", data=data)

    track = None
    if okay:
        # Get the track from the dataset
        track = dataset.get_track(data["dataset_index"])

    return track, okay


def map_ready(*args):

    global map

    map = gui.popup_window["track_selector"].find_element_by_id("track_selector_map").widget

    map.jump_to(0.0, 0.0, 1)

    # Container layer
    layer = map.add_layer("track_selector")

    # And add a layer with a circle for the search radius
    layer.add_layer(
        "search_radius",
        type="line",
        line_color="white",
        line_style="--",
    )

    layer.add_layer("tracks",
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
    update_circle()


def update_tracks():

    data = gui.popup_data

    tdb = data["track_selector"]["track_dataset"]

    # Get filter data
    distance = gui.getvar("cyclone_track_selector", "distance")
    year_min = gui.getvar("cyclone_track_selector", "year0")
    year_max = gui.getvar("cyclone_track_selector", "year1")
    vmax_min = gui.getvar("cyclone_track_selector", "vmax0")
    vmax_max = gui.getvar("cyclone_track_selector", "vmax1")
    lon = gui.getvar("cyclone_track_selector", "lon")
    lat = gui.getvar("cyclone_track_selector", "lat")

    # Get indices based on filter
    storm_indices = tdb.filter(
        lon=lon, lat=lat, distance=distance, year_min=year_min, year_max=year_max,
        vmax_min=vmax_min, vmax_max=vmax_max
    )

    if len(storm_indices) == 0:
        # No storms found, clear the layer
        map.layer["track_selector"].layer["tracks"].set_data(GeoDataFrame(), 0)
        gui.setvar("cyclone_track_selector", "storm_names", [""])
        gui.setvar("cyclone_track_selector", "storm_indices", [0])
        gui.setvar("cyclone_track_selector", "storm_index", 0)

    else:

        # Get GeoDataFrame of tracks
        gdf = tdb.to_gdf(index=storm_indices)

        # Also get list with names of datasets
        storm_names = tdb.list_names(index=storm_indices)

        # storm_indices is now a numpy array, convert to list
        storm_indices = list(storm_indices)

        # We want the names in alphabetical order, so sort storm_names AND storm_indices accordingly, using zip
        storm_names, storm_indices = zip(*sorted(zip(storm_names, storm_indices)))

        gui.setvar("cyclone_track_selector", "storm_names", storm_names)
        gui.setvar("cyclone_track_selector", "storm_indices", storm_indices)

        # Get the original storm name
        storm_index = gui.getvar("cyclone_track_selector", "storm_index")

        # If storm_index is in the list, set it as selected. Need to find the index in the list.
        if storm_index in storm_indices:
            selected_index = storm_indices.index(storm_index)
        else:
            selected_index = 0

        gui.setvar("cyclone_track_selector", "storm_index", storm_indices[selected_index])        

        map.layer["track_selector"].layer["tracks"].set_data(gdf, selected_index)

    gui.popup_window["track_selector"].update()

def update_circle():
    lon = gui.getvar("cyclone_track_selector", "lon")
    lat = gui.getvar("cyclone_track_selector", "lat")
    distance = gui.getvar("cyclone_track_selector", "distance")
    # Create a gdf with a circle around the lon, lat with the given radius, using exact geospatial operations
    # Make circle in web mercator coordinates
    transformer = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
    x, y = transformer.transform(lon, lat)
    center = Point(x, y)
    circle = center.buffer(distance * 1000)  # Convert km to meters
    # Convert circle to GeoDataFrame in web mercator coordinates
    gdf = GeoDataFrame(geometry=[circle], crs=3857).to_crs(4326)
    map.layer["track_selector"].layer["search_radius"].set_data(gdf)


def map_moved(coords, widget):
    pass

def select_track(feature, widget):
    storm_index = feature["properties"]["dataset_index"]
    gui.popup_data["track_selector"]["dataset_index"] = storm_index
    gui.setvar("cyclone_track_selector", "storm_index", storm_index)
    gui.popup_window["track_selector"].update()

def select_track_from_list(*args):
    storm_index = gui.getvar("cyclone_track_selector", "storm_index")
    storm_indices = gui.getvar("cyclone_track_selector", "storm_indices")
    index = storm_indices.index(storm_index)
    map.layer["track_selector"].layer["tracks"].set_selected_index(index)

def edit_filter(val, widget):
    update_tracks()
    update_circle()

def reset_search_location(*args):
    crds = map.map_center
    # Round lon and lat to 2 decimals
    lon = round(crds[0], 2)
    lat = round(crds[1], 2)
    gui.setvar("cyclone_track_selector", "lon", lon)
    gui.setvar("cyclone_track_selector", "lat", lat)
    update_tracks()
    update_circle()
