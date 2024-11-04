# Load modules
import os

from cht_cyclones.tropical_cyclone import TropicalCyclone, TropicalCycloneEnsemble

# Define a track
tc                  = TropicalCyclone(name="Maria_forecast")

# Read a track
current_directory   = os.path.dirname(os.path.abspath(__file__))
file_path           = os.path.join(current_directory, 'TRK_COAMPS_CTCX_3_2017092006_15L')
tc.read_track(file_path,'trk')

# Define settings
# Uses typical Wind Enhance Scheme (WES) formulations
tc.include_rainfall     = True      # using empirical formulation to compute rainfall (IPET default)
tc.rainfall_factor      = 2.0       # factor to calibrate rainfall 2.0 means 2x as much as relationship

# Write as cyc
tc.convert_units_metric_imperial()
file_path           = os.path.join(current_directory, 'forecasted_track_Maria.cyc')
tc.write_track(file_path, 'ddb_cyc')

# create (regular) ASCII spiderweb
file_path           = os.path.join(current_directory, 'forecasted_track_Maria.spw')
tc.to_spiderweb(file_path)

# create (netcdf) spiderweb
file_path           = os.path.join(current_directory, 'forecasted_track_Maria.nc')
tc.to_spiderweb(file_path, format_type='netcdf')

# Write out as shapefile
file_path           = os.path.join(current_directory, 'forecasted_track_Maria.shp')
tc.track.to_file(file_path)

# Write out as geojson
file_path           = os.path.join(current_directory, 'forecasted_track_Maria.json')
tc.track.to_file(file_path, driver="GeoJSON")

# TC-FF ensemble members
tc2                     = TropicalCycloneEnsemble(name="Maria_forecast_ensembles", TropicalCyclone=tc)
tc2.dt                  = 3                             # 3 hourly
tc2.include_best_track  = 0                             # not included best-track
tc2.compute_ensemble(100)
tc2.to_shapefile(current_directory)