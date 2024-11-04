# Load modules
import os

from cht_cyclones.tropical_cyclone import TropicalCyclone, TropicalCycloneEnsemble

# Define a track
tc                  = TropicalCyclone(name="Idai")

# Read a track
current_directory   = os.path.dirname(os.path.abspath(__file__))
file_path           = os.path.join(current_directory, 'best_track_idai.cyc')
tc.from_ddb_cyc(file_path)

# Define settings
# Uses typical Wind Enhance Scheme (WES) formulations
tc.include_rainfall     = True      # using empirical formulation to compute rainfall (IPET default)
tc.rainfall_factor      = 2.0       # factor to calibrate rainfall 2.0 means 2x as much as relationship

# create (regular) ASCII spiderweb
file_path           = os.path.join(current_directory, 'best_track_idai.spw')
tc.to_spiderweb(file_path)

# create (netcdf) spiderweb
file_path           = os.path.join(current_directory, 'best_track_idai.nc')
tc.to_spiderweb(file_path, format_type='netcdf')

# Write out as shapefile
file_path           = os.path.join(current_directory, 'best_track_idai.shp')
tc.track.to_file(file_path)

# Write out as geojson
file_path           = os.path.join(current_directory, 'best_track_idai.json')
tc.track.to_file(file_path, driver="GeoJSON")

# Write as cyc again
tc.convert_units_metric_imperial()
file_path           = os.path.join(current_directory, 'best_track_idai_cht.cyc')
tc.write_track(file_path, 'ddb_cyc')