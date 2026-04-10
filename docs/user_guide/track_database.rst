Track database
==============

``cht_cyclones`` organises historical storm tracks using a two-level hierarchy:

* :class:`~cht_cyclones.track_database.CycloneTrackDatabase` — top-level
  container that manages *multiple* datasets and can synchronise them from a
  remote S3 bucket.
* :class:`~cht_cyclones.track_dataset.CycloneTrackDataset` — a single labelled
  dataset (e.g. the IBTrACS global best-track archive) backed by a NetCDF file.

Directory layout
----------------

A database root directory must contain a TOML catalogue file
``cyclone_track_datasets.tml``.  Each entry points to a sub-directory that
holds a ``metadata.tml`` file and the actual data file(s):

.. code-block:: text

   /cyclone_database/
   ├── cyclone_track_datasets.tml   # catalogue listing all datasets
   └── ibtracs_global/
       ├── metadata.tml             # dataset metadata
       └── IBTrACS.ALL.v04r00.nc   # IBTrACS NetCDF

The catalogue file lists datasets by name and (optionally) explicit path:

.. code-block:: toml

   [[dataset]]
   name = "ibtracs_global"

The metadata file describes the dataset format and S3 source:

.. code-block:: toml

   format     = "ibtracs"
   long_name  = "IBTrACS v4 – Global"
   files      = ["IBTrACS.ALL.v04r00.nc"]
   s3_bucket  = "my-bucket"
   s3_key     = "cyclones/ibtracs_global"

Opening a database
------------------

.. code-block:: python

   from cht_cyclones import CycloneTrackDatabase

   db = CycloneTrackDatabase(
       path="/path/to/cyclone_database",
       s3_bucket="my-bucket",
       s3_key="cyclones",
       check_online=True,   # check S3 for new datasets on startup
   )

   # All loaded datasets
   for ds in db.dataset:
       print(ds.name, ds.long_name)

When ``check_online=True`` the database queries S3 for new datasets and
downloads their metadata files.  The actual track data is *not* downloaded
until :meth:`~cht_cyclones.track_dataset.CycloneTrackDataset.download` or
:meth:`~cht_cyclones.track_dataset.CycloneTrackDataset.read` is called.

Downloading data
----------------

.. code-block:: python

   dataset = db.dataset[0]
   dataset.download()   # fetches missing files from S3

Reading the track data into memory is deferred until the first filtering or
retrieval call, so loading a database is fast even when the NetCDF files are
large.

Supported formats
-----------------

IBTrACS
~~~~~~~

`IBTrACS <https://www.ncei.noaa.gov/products/international-best-track-archive>`_
(International Best Track Archive for Climate Stewardship) provides 6-hourly
best-track positions and intensities for all recorded tropical cyclones
worldwide from 1842 onwards.  The ``usa_wind`` and ``usa_pres`` variables are
used for intensity, and the ``usa_r34`` / ``usa_r50`` / ``usa_r64`` variables
provide quadrant wind radii where available.

Specify ``format = "ibtracs"`` in ``metadata.tml``.
