Creating wind fields
====================

``cht_cyclones`` generates parametric 2-D wind and pressure fields from the
track using the Holland et al. (2010) radial profile.  The output is a
*spiderweb* — a time-varying field on a polar (range × azimuth) grid centred
on the storm eye at each time step.

Wind-field configuration
------------------------

All settings are stored in ``tc.config`` (a plain Python dict).  The defaults
are suitable for most studies; change only what you need:

.. code-block:: python

   from cht_cyclones import TropicalCyclone

   tc = TropicalCyclone(track_file="storm.cyc")

   # Grid geometry
   tc.config["spiderweb_radius"]    = 400.0   # outer radius (km)
   tc.config["nr_radial_bins"]      = 100     # radial resolution
   tc.config["nr_directional_bins"] = 36      # 10° azimuthal bins

   # Physics
   tc.config["wind_profile"]         = "holland2010"   # only supported option
   tc.config["wind_pressure_relation"] = "holland2008"
   tc.config["rmw_relation"]         = "nederhoff2019"
   tc.config["background_pressure"]  = 1012.0   # hPa
   tc.config["wind_conversion_factor"] = 0.93   # 1-min to 10-min conversion

   # Asymmetry
   tc.config["asymmetry_option"] = "schwerdt1979"  # or "mvo" or "none"
   tc.config["phi_spiral"]       = 20.0   # inflow angle (degrees)
   tc.config["phi_trans"]        = 45.0   # translation-vector angle offset

   # Rainfall
   tc.config["include_rainfall"] = True
   tc.config["rainfall_relationship"] = "ipet"

   # Reference time for the spiderweb NetCDF output
   tc.config["tref"] = "20000101 000000"

Configuration can also be loaded from a TOML file:

.. code-block:: python

   tc.read_config("my_settings.toml")
   tc.write_config("output_settings.toml")   # save current settings

Holland et al. (2010) profile
------------------------------

The wind speed at radius *r* from the storm centre is:

.. math::

   V(r) = V_{\max}
          \left( \frac{R_{\max}}{r} \right)^B
          \exp\!\left[1 - \left(\frac{R_{\max}}{r}\right)^B \right]^x

where *B* is the Holland shape parameter derived from the pressure–wind
relationship (Holland, 2008) and *x* varies between 0.5 (inside *R*\ :sub:`max`)
and a fitted value outside it.

When quadrant wind radii (R35, R50, R65, R100) are present in the track data
the profile exponent *x* and the asymmetry amplitude are fitted to match the
observations at each time step.  When no radii are available, *R*\ :sub:`max`
is estimated from the Nederhoff et al. (2019) statistical relationship.

Computing and writing the wind field
-------------------------------------

.. code-block:: python

   # Compute — prints one line per track point
   tc.compute_wind_field()

   # Write Delft3D ASCII spiderweb (.spw)
   tc.write_spiderweb("storm.spw")

   # Write NetCDF spiderweb
   tc.write_spiderweb("storm.nc", format="netcdf")

The resulting ``tc.spiderweb`` object is a
:class:`~cht_cyclones.spiderweb.TropicalCycloneSpiderweb` whose underlying
data live in the :class:`xarray.Dataset` at ``tc.spiderweb.ds``.  Variables
are described in :doc:`spiderweb_format`.

Reading an existing spiderweb
------------------------------

.. code-block:: python

   from cht_cyclones.spiderweb import TropicalCycloneSpiderweb

   spw = TropicalCycloneSpiderweb()
   spw.read("storm.spw")        # auto-detects .spw vs .nc by extension
   print(spw.ds)                # xarray.Dataset

Wind-profile functions
----------------------

The underlying one-dimensional profile functions are available as standalone
utilities for scripting and research:

.. code-block:: python

   import numpy as np
   from cht_cyclones.wind_profiles import holland2010, wind_radii_nederhoff

   # Radial wind speed profile
   r = np.linspace(1, 500, 500)   # km
   vr, pr = holland2010(
       r=r,
       vmax=50.0,    # m/s
       pc=945.0,     # hPa
       pn=1012.0,    # hPa
       rmax=40.0e3,  # m  (40 km)
       dpdt=0.0,
       lat=20.0,
       vt=5.0,       # m/s forward speed
       xn=0.5,
   )

   # Estimate RMW and R35 from Nederhoff et al. (2019)
   rmax_stats, dr35_stats = wind_radii_nederhoff(
       vmax=50.0,     # m/s
       lat=20.0,
       region=7,      # 7 = global
       probability=0, # 0 = mode only
   )
   print(rmax_stats["mode"])  # km
