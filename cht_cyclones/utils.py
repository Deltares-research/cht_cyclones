"""
Miscellaneous utility helpers for the cht_cyclones package.

Provides functions for serialising GeoDataFrames to JavaScript-wrapped GeoJSON
files and to Delft3D polyline (*.pli) files.
"""

import json


def gdf_to_geojson_js(gdf, filename: str, varname: str = "track_ensemble") -> None:
    """
    Write a GeoDataFrame to a JavaScript-wrapped GeoJSON file.

    The output file contains a single variable assignment of the form
    ``var <varname> = { ... }``.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Source data.
    filename : str
        Path to the output ``.js`` file.
    varname : str, optional
        JavaScript variable name used in the assignment statement.
    """
    text = f"var {varname} ="
    write_json_js(filename, json.loads(gdf.to_json()), text)


def write_json_js(file_name: str, jsn, first_line: str) -> None:
    """
    Write a JSON object or list to a JavaScript file with a leading assignment line.

    Parameters
    ----------
    file_name : str
        Path to the output file.
    jsn : dict or list
        JSON-serialisable data to write.
    first_line : str
        Text to write on the first line (typically a ``var x =`` assignment).
    """
    if isinstance(jsn, list):
        with open(file_name, "w") as f:
            f.write(f"{first_line}\n")
            f.write("[\n")
            for ix, x in enumerate(jsn):
                json_string = json.dumps(x)
                if ix < len(jsn) - 1:
                    f.write(json_string + ",")
                else:
                    f.write(json_string)
                f.write("\n")
            f.write("]\n")
    else:
        with open(file_name, "w") as f:
            f.write(f"{first_line}\n")
            json_string = json.dumps(jsn)
            f.write(f"{json_string}\n")


def gdf_to_pli(gdf, filename: str) -> None:
    """
    Write all LineString and Polygon geometries from a GeoDataFrame to a ``.pli`` file.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Source data; only ``LineString`` and ``Polygon`` geometries are written.
    filename : str
        Path to the output ``.pli`` file (opened in append mode).
    """
    for _irow, row in gdf.iterrows():
        if row["geometry"].geom_type == "LineString":
            nrp = len(row["geometry"].coords)
            with open(filename, "a") as f:
                f.write("BL01\n")
                f.write(f"{nrp} 2\n")
                for ip in range(nrp):
                    x = row["geometry"].coords[ip][0]
                    y = row["geometry"].coords[ip][1]
                    f.write(f"{x} {y}\n")
        elif row["geometry"].geom_type == "Polygon":
            nrp = len(row["geometry"].exterior.coords)
            with open(filename, "a") as f:
                f.write("BL01\n")
                f.write(f"{nrp} 2\n")
                for ip in range(nrp):
                    x = row["geometry"].exterior.coords[ip][0]
                    y = row["geometry"].exterior.coords[ip][1]
                    f.write(f"{x} {y}\n")
