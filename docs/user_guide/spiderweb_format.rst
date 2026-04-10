Spiderweb format
================

A *spiderweb* is the primary output of ``cht_cyclones``: a time-varying 2-D
wind-field on a polar (range × azimuth) grid centred on the storm eye.
Hydrodynamic models such as Delft3D-FM, SFINCS, and HurryWave can read
spiderweb files directly as meteorological forcing.

Grid structure
--------------

The grid has three dimensions:

* **time** — one snapshot per track point
* **range** — radial distance from the eye; equally spaced from the first bin
  centre to ``spiderweb_radius`` (default 400 km)
* **azimuth** — direction in nautical convention (clockwise from North);
  equally spaced in steps of ``360 / nr_directional_bins`` degrees

Variables stored at each (time, range, azimuth) cell:

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - Variable
     - Units
     - Description
   * - ``wind_x``
     - m/s
     - East component of wind velocity
   * - ``wind_y``
     - m/s
     - North component of wind velocity
   * - ``pressure``
     - Pa
     - Atmospheric pressure
   * - ``precipitation``
     - mm/h
     - Rainfall rate (if ``include_rainfall = True``)
   * - ``longitude_eye``
     - degrees
     - Longitude of the storm eye at each time step
   * - ``latitude_eye``
     - degrees
     - Latitude of the storm eye at each time step

ASCII (``.spw``) format
-----------------------

The Delft3D ASCII format is the most widely supported output.  The file begins
with a metadata header followed by one block per time step:

.. code-block:: text

   FileVersion                  = 1.03
   n_cols                       = 36
   n_rows                       = 100
   grid_unit                    = m
   spw_radius                   = 400000.0
   spw_rad_unit                 = m
   n_quantity                   = 3
   quantity1                    = wind_speed
   unit1                        = m/s
   ...
   NODATA_value                 = -999.0

   TIME = 2005-08-26 00:00:00 +00:00
   Longitude of centre          = -75.1
   Latitude of centre           = 23.4
   <data rows: wind_x, wind_y, pressure, [precipitation]>

NetCDF format
-------------

The NetCDF output stores all variables in a single self-describing file and is
the recommended format for archival and post-processing:

.. code-block:: python

   import xarray as xr

   ds = xr.open_dataset("katrina.nc")
   print(ds)

   # Extract wind speed at a single time step
   import numpy as np
   wx = ds["wind_x"].isel(time=0).values
   wy = ds["wind_y"].isel(time=0).values
   speed = np.hypot(wx, wy)

Reading a spiderweb
-------------------

Use :class:`~cht_cyclones.spiderweb.TropicalCycloneSpiderweb` to read either
format; the format is inferred from the file extension (``.spw`` or ``.nc``):

.. code-block:: python

   from cht_cyclones.spiderweb import TropicalCycloneSpiderweb

   spw = TropicalCycloneSpiderweb()
   spw.read("katrina.spw")
   print(spw.ds["wind_x"])

Model compatibility
-------------------

Both formats are accepted by:

* **Delft3D-FM** (``*.spw``)
* **SFINCS** (``*.spw`` or ``*.nc``)
* **HurryWave** (``*.spw``)
* **Delft3D 4** (``*.spw`` — classic structured grid version)

Refer to each model's reference manual for the keyword to specify the spiderweb
file path in the model input.
