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
