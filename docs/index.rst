cht_cyclones documentation
###########################

**cht_cyclones** is a Python toolkit for working with tropical cyclone data in
coastal hazard studies.  It is part of the Deltares `Coastal Hazards Toolkit
<https://github.com/Deltares-research>`_ (CHT) family of packages.

Key capabilities
----------------

* **Historical track databases** — load and filter IBTrACS global best-track
  data; filter by region, year, basin, and storm intensity.
* **Parametric wind-field generation** — compute 2-D wind and pressure fields
  using the Holland et al. (2010) profile fitted to observed wind radii.
* **Spiderweb output** — write time-varying wind fields on a polar (range ×
  azimuth) grid in Delft3D ASCII (``.spw``) or NetCDF format, ready for use in
  SFINCS, Delft3D-FM, HurryWave, and other hydrodynamic models.
* **Probabilistic ensembles** — generate perturbed track and intensity
  realisations based on NHC forecast-error statistics (DeMaria et al., 2009).
* **Multiple track formats** — read and write ``.cyc``, ``ddb_cyc``,
  COAMPS-TC ``.trk``, PAGASA CSV, JTWC advisory text, and JMV 3.0 formats.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   getting_started
   user_guide/index
   api/index
   changelog
