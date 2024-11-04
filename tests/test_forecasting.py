# Load modules
import os
from datetime import datetime

from cht_cyclones.tropical_cyclone import TropicalCyclone, TropicalCycloneEnsemble

# Define folder and track
tc                  = TropicalCyclone(name="Idai")
current_directory   = os.path.dirname(os.path.abspath(__file__))

# Read a track
file_path           = os.path.join(current_directory, 'best_track_idai.cyc')
tc.from_ddb_cyc(file_path)
tc.account_for_forward_speed()
tc.estimate_missing_values()
tc.include_rainfall = True

# Forecasting
# Uses DeMaria et al. 2009 formulations
tc2                     = TropicalCycloneEnsemble(name="Idai_forecastng", TropicalCyclone=tc)
tc2.dt                  = 3                             # 3 hourly
tc2.include_best_track  = 0                             # not included best-track
tc2.tstart              = datetime(2019,3,9,0,0,0)      # March 5, 2019 at midnight UTC
tc2.tend                = datetime(2019,3,16,6,0,0)     # last timestamp in cyc
nm_to_m                 = float(1.852) * 1000           # nm to m
tc2.mean_abs_cte24      = 19.0397*nm_to_m               # mean absolute error in cross-track error (CTE) in meter
tc2.sc_cte              = 1.3253                        # auto-regression CTE; typically >1
tc2.compute_ensemble(10)

# Write out
tc2.to_shapefile(current_directory)
tc2.to_spiderweb(current_directory)
tc2.make_figures(current_directory)