Getting started
===============

Installation
------------

Install the latest release from PyPI:

.. code-block:: bash

   pip install cht_cyclones

For development, clone the repository and install in editable mode:

.. code-block:: bash

   git clone https://github.com/Deltares-research/cht_cyclones.git
   cd cht_cyclones
   pip install -e ".[tests]"

Dependencies
~~~~~~~~~~~~

``cht_cyclones`` requires Python 3.10 or later.  The following packages are
installed automatically: ``numpy``, ``scipy``, ``pandas``, ``geopandas``,
``matplotlib``, ``shapely``, ``xarray``, ``fiona``, ``toml``, ``boto3``,
``python-dateutil``, ``cht_utils``, ``geojson``, ``netCDF4``, and
``feedparser``.

Quick example
-------------

The workflow below shows the most common use case: loading a historical storm
from the IBTrACS database and generating a spiderweb wind field for use in a
hydrodynamic model.

1. Load a historical track from IBTrACS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The :class:`~cht_cyclones.track_database.CycloneTrackDatabase` class manages
one or more local track datasets.  A minimal database directory contains a
``cyclone_track_datasets.tml`` catalogue file and one or more dataset
sub-directories, each with a ``metadata.tml`` file.

.. code-block:: python

   from cht_cyclones import CycloneTrackDatabase

   # Point to the local directory that holds your IBTrACS data
   db = CycloneTrackDatabase("/path/to/cyclone_database")

   # Access the first (or only) dataset
   dataset = db.dataset[0]

   # Download data files from S3 if they are not yet available locally
   dataset.download()

   # Filter by proximity to a point, year range, and minimum intensity
   indices = dataset.filter(
       lon=90.0,
       lat=20.0,
       distance=1500.0,   # km
       year_min=2000,
       year_max=2023,
       vmax_min=64.0,     # knots – Category 1 threshold
   )

   # Retrieve the track as a TropicalCyclone object
   tc = dataset.get_track(indices[0])
   print(tc.name)

2. Inspect the track
~~~~~~~~~~~~~~~~~~~~

The track is stored as a :class:`geopandas.GeoDataFrame` at
``tc.track.gdf``.  Each row represents one 6-hourly (or finer) track point.

.. code-block:: python

   print(tc.track.gdf[["datetime", "vmax", "pc", "rmw"]].head())

Columns include:

* ``datetime`` — timestamp string ``"YYYYMMDD HHMMSS"``
* ``vmax`` — maximum sustained wind speed (knots)
* ``pc`` — central pressure (hPa)
* ``rmw`` — radius of maximum winds (nautical miles)
* ``r35_ne/se/sw/nw`` — 35-knot wind radii by quadrant (nautical miles)
* ``r50_*``, ``r65_*``, ``r100_*`` — 50-, 65-, and 100-knot wind radii

3. Generate a spiderweb wind field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:meth:`~cht_cyclones.tropical_cyclone_refactored.TropicalCyclone.compute_wind_field`
converts the track into a time-varying polar wind-field grid:

.. code-block:: python

   from cht_cyclones import TropicalCyclone

   tc = TropicalCyclone(track_file="katrina.cyc")

   # Optionally adjust configuration
   tc.config["spiderweb_radius"] = 500.0       # km
   tc.config["nr_radial_bins"]   = 100
   tc.config["nr_directional_bins"] = 36
   tc.config["wind_profile"]     = "holland2010"
   tc.config["include_rainfall"] = True

   # Compute the wind field (populates tc.spiderweb)
   tc.compute_wind_field()

   # Write Delft3D ASCII spiderweb
   tc.write_spiderweb("katrina.spw")

   # Or write NetCDF
   tc.write_spiderweb("katrina.nc", format="netcdf")

4. Save and reload a track file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Write the track in the native .cyc format
   tc.write_track("katrina.cyc")

   # Read it back
   tc2 = TropicalCyclone(track_file="katrina.cyc")
