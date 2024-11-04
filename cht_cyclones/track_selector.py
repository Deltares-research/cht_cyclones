import os


def track_selector(dataset, app, lon=0.0, lat=0.0, distance=1000.0, year_min=1850, year_max=2030):

    app.gui.setvar("cyclone_track_selector", "lon", lon)
    app.gui.setvar("cyclone_track_selector", "lat", lat)
    app.gui.setvar("cyclone_track_selector", "distance", distance)
    app.gui.setvar("cyclone_track_selector", "year0", year_min)
    app.gui.setvar("cyclone_track_selector", "year1", year_max)
    app.gui.setvar("cyclone_track_selector", "name", "")

    data = {}
    data["track_dataset"] = dataset

    # Read GUI config file
    config_file = os.path.join(os.path.dirname(__file__), "cyclone_track_selector.yml")
    okay, data = app.gui.popup(config_file, id="track_selector", data=data)

    track = None
    if okay:
        # Get the track from the dataset
        track = dataset.get_track(data["dataset_index"])

    return track, okay

def map_ready(widget):

    print("Selector map is ready !")

    gui = widget.element.gui

    mp = gui.popup_window["track_selector"].find_element_by_id("track_selector_map").widget
    mp.jump_to(0.0, 0.0, 1)
    data = gui.popup_data
    # Container layers
    data["track_selector"]["main_layer"] = mp.add_layer("track_selector")
    # Tracks layers
    data["track_selector"]["track_layer"] = data["track_selector"]["main_layer"].add_layer(
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
    update_tracks(gui)


def update_tracks(gui):

    data = gui.popup_data

    tdb = data["track_selector"]["track_dataset"]
    tracks_layer = data["track_selector"]["track_layer"]

    # Get filter data
    distance = gui.getvar("cyclone_track_selector", "distance")
    year_min = gui.getvar("cyclone_track_selector", "year0")
    year_max = gui.getvar("cyclone_track_selector", "year1")
    lon      = gui.getvar("cyclone_track_selector", "lon")
    lat      = gui.getvar("cyclone_track_selector", "lat")

    # Get indices based on filter
    index = tdb.filter(
        lon=lon, lat=lat, distance=distance, year_min=year_min, year_max=year_max
    )

    # Get GeoDataFrame of tracks
    gdf = tdb.to_gdf(index=index)

    tracks_layer.set_data(gdf, 0)


def map_moved(coords, widget):
    pass


def select_track(feature, widget):
    widget.element.gui.popup_data["track_selector"]["dataset_index"] = feature["properties"]["dataset_index"]


def edit_filter(val, widget):
    update_tracks(widget.element.gui)
