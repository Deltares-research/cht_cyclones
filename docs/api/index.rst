API reference
=============

This page documents all public classes and functions exposed by
``cht_cyclones``.  The primary public API is importable directly from the
package:

.. code-block:: python

   from cht_cyclones import (
       TropicalCyclone,
       TropicalCycloneTrack,
       CycloneTrackDatabase,
       CycloneTrackDataset,
       track_selector,
   )

----

TropicalCyclone
---------------

.. autoclass:: cht_cyclones.tropical_cyclone_refactored.TropicalCyclone
   :members:
   :undoc-members:
   :show-inheritance:

----

TropicalCycloneTrack
--------------------

.. autoclass:: cht_cyclones.track.TropicalCycloneTrack
   :members:
   :undoc-members:
   :show-inheritance:

----

TropicalCycloneSpiderweb
------------------------

.. autoclass:: cht_cyclones.spiderweb.TropicalCycloneSpiderweb
   :members:
   :undoc-members:
   :show-inheritance:

----

TropicalCycloneEnsemble
-----------------------

.. autoclass:: cht_cyclones.ensemble.TropicalCycloneEnsemble
   :members:
   :undoc-members:
   :show-inheritance:

----

CycloneTrackDatabase
--------------------

.. autoclass:: cht_cyclones.track_database.CycloneTrackDatabase
   :members:
   :undoc-members:
   :show-inheritance:

----

CycloneTrackDataset
-------------------

.. autoclass:: cht_cyclones.track_dataset.CycloneTrackDataset
   :members:
   :undoc-members:
   :show-inheritance:

----

track_selector
--------------

.. autofunction:: cht_cyclones.track_selector.track_selector

----

Wind-profile functions
----------------------

.. autofunction:: cht_cyclones.wind_profiles.holland2010

.. autofunction:: cht_cyclones.wind_profiles.wind_radii_nederhoff
