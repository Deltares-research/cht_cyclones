import json


def gdf_to_geojson_js(gdf, filename, varname="track_ensemble"):
    """Write out the ensemble to a geojson or js file."""
    text = "var " + varname + " ="
    write_json_js(filename, json.loads(gdf.to_json()), text)

def write_json_js(file_name, jsn, first_line):
     # Writes json javascript file
     if type(jsn) == list:
         f = open(file_name, "w")
         f.write(first_line + "\n")
         f.write("[\n")
         for ix, x in enumerate(jsn):
             json_string = json.dumps(x)     
             if ix<len(jsn) - 1:
                 f.write(json_string + ",")
             else:
                 f.write(json_string)
             f.write('\n')
         f.write("]\n")
         f.close()
     else:
         f = open(file_name, "w")
         f.write(first_line + "\n")
         json_string = json.dumps(jsn)     
         f.write(json_string + "\n")
         f.close()         

def gdf_to_pli(gdf, filename):
    # Write all linstrings and polygons to pli file
    # Loop through all rows in gdf
    for irow, row in gdf.iterrows():
        # If geometry is linestring or polygon, extract coordinates x and y into list and write to text file
        if row["geometry"].geom_type == "LineString":
            # Write to file
            nrp = len(row["geometry"].coords)
            with open(filename, "a") as f:
                f.write("BL01\n")
                f.write(f"{nrp} 2\n")
                for ip in range(nrp):
                    f.write(str(row["geometry"].coords[ip][0]) + " " + str(row["geometry"].coords[ip][1]) + "\n")
        elif row["geometry"].geom_type == "Polygon":
            # Write to file
            nrp = len(row["geometry"].exterior.coords)
            with open(filename, "a") as f:
                f.write("BL01\n")
                f.write(f"{nrp} 2\n")
                for ip in range(nrp):
                    f.write(str(row["geometry"].exterior.coords[ip][0]) + " " + str(row["geometry"].exterior.coords[ip][1]) + "\n")
