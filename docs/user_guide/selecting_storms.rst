Selecting storms
================

Once a :class:`~cht_cyclones.track_dataset.CycloneTrackDataset` is loaded you
can filter the full catalogue and retrieve individual tracks.

Filtering the dataset
---------------------

:meth:`~cht_cyclones.track_dataset.CycloneTrackDataset.filter` returns an
array of integer indices that satisfy *all* of the supplied criteria:

.. code-block:: python

   dataset = db.dataset[0]

   indices = dataset.filter(
       lon=90.5,          # degrees East (search centre)
       lat=22.0,          # degrees North
       distance=1000.0,   # km – only tracks passing within this radius
       year_min=1990,
       year_max=2023,
       vmax_min=64.0,     # knots – exclude weak tropical depressions/storms
       vmax_max=1000.0,
   )
   print(f"{len(indices)} storms match the filter")

All parameters are optional; omitting them includes all storms on that
criterion.

Listing matching storms
-----------------------

:meth:`~cht_cyclones.track_dataset.CycloneTrackDataset.list_names` returns a
list of storm names for a given set of indices:

.. code-block:: python

   names = dataset.list_names(index=indices)
   for i, name in zip(indices, names):
       print(i, name)

Retrieving a single track
-------------------------

:meth:`~cht_cyclones.track_dataset.CycloneTrackDataset.get_track` converts one
IBTrACS record into a
:class:`~cht_cyclones.tropical_cyclone_refactored.TropicalCyclone` object:

.. code-block:: python

   tc = dataset.get_track(indices[0])
   print(tc.name)
   print(tc.track.gdf.head())

The returned :class:`~cht_cyclones.tropical_cyclone_refactored.TropicalCyclone`
can then be used directly for wind-field computation or written to a file.

Interactive GUI selector (DelftDashboard)
-----------------------------------------

When running inside DelftDashboard, the
:func:`~cht_cyclones.track_selector.track_selector` helper opens a popup window
with an interactive map, filter controls, and a name list:

.. code-block:: python

   from cht_cyclones import track_selector

   tc, ok = track_selector(
       dataset=db.dataset[0],
       app=app,                # DelftDashboard app object
       lon=90.0,
       lat=20.0,
       distance=1000.0,
       year_min=2000,
       year_max=2023,
   )
   if ok:
       tc.compute_wind_field()

This is intended for interactive use only and requires the ``guitars``
framework to be available.

Reading individual track files
------------------------------

Tracks can also be loaded from stand-alone files:

.. code-block:: python

   from cht_cyclones import TropicalCyclone

   # Auto-detect format (.cyc / ddb_cyc / COAMPS-TC .trk / PAGASA CSV)
   tc = TropicalCyclone(track_file="my_storm.cyc")

   # Explicit format
   tc = TropicalCyclone()
   tc.read_track("advisory.txt", format="jtwc_advisory_text")

Supported formats
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Format identifier
     - Description
   * - ``cyc``
     - Native ``.cyc`` format (auto-detects legacy ``ddb_cyc`` and
       COAMPS-TC ``.trk`` variants)
   * - ``ddb_cyc``
     - Old DelftDashboard ``.cyc`` format with embedded config header
   * - ``jmv30``
     - JMV 3.0 advisory format (``.tcw``)
   * - ``jtwc_advisory_text``
     - JTWC operational advisory text file

Writing tracks
--------------

.. code-block:: python

   # Native .cyc
   tc.write_track("storm.cyc")

   # DelftDashboard format with embedded config
   tc.write_track("storm.cyc", format="ddb_cyc")

Merging track files
-------------------

Forecast products are often split across multiple hourly files.
Pass a list of paths to merge them chronologically:

.. code-block:: python

   tc = TropicalCyclone(
       track_file=["track_00h.trk", "track_06h.trk", "track_12h.trk"],
       tau=6,   # discard the first 6 h of each file after the first
   )
