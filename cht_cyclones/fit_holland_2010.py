"""
Optimisation routines that fit Holland (2010) wind-field parameters to observed
wind-radii quadrant data.
"""

import numpy as np
from scipy.interpolate import interp1d

from cht_cyclones.wind_profiles import holland2010


def fit_wind_field_holland2010(
    vmax: float,
    rmax: float,
    pc: float,
    vtreal: float,
    phit: float,
    pn: float,
    phi_spiral: float,
    lat: float,
    dpdt: float,
    obs: dict,
) -> list:
    """
    Fit Holland (2010) wind-field parameters to observed quadrant wind radii.

    Iteratively adjusts the profile exponent *xn*, the asymmetry magnitude *vt*,
    and the asymmetry direction *phia* to minimise the RMS error between the
    modelled wind field and the observed quadrant radii.

    Parameters
    ----------
    vmax : float
        Maximum sustained wind speed (m/s).
    rmax : float
        Radius of maximum winds (km).
    pc : float
        Central pressure (hPa).
    vtreal : float
        Actual storm translation speed (m/s).
    phit : float
        Track heading angle (degrees, Cartesian convention).
    pn : float
        Background pressure (hPa).
    phi_spiral : float
        Inflow angle (degrees).
    lat : float
        Storm latitude (degrees).
    dpdt : float
        Rate of pressure change (hPa/h).
    obs : dict
        Observed data with keys ``"quadrants_radii"`` (ndarray, shape (4, 4))
        and ``"quadrants_speed"`` (ndarray, length 4).

    Returns
    -------
    xn : float
        Fitted Holland profile exponent.
    vt : float
        Fitted asymmetry-correction wind speed (m/s).
    phia : float
        Fitted asymmetry direction (degrees).
    """
    # function to fit wind field based on Holland 2010
    size_factor = 1
    phi = np.arange(90, -270 - 10, -10)  # radial angles (cartesian, degrees)
    rmax = rmax * 1000  # convert from km to m
    r = np.arange(5000, 500000 + 5000, 5000)

    # first estimates
    xn = 0.5
    vt = 0.6 * vtreal

    if lat > 0:
        phia = 45  # angle with respect to track angle (cartesian degrees, i.e. counter-clockwise)
    else:
        phia = -45

    # More variables
    dxn = 0.01
    dvt = 0.5
    dphia = 5
    nrad = 2
    nobs = 0

    for irad in range(np.size(obs["quadrants_radii"], 0)):
        for iquad in range(np.size(obs["quadrants_radii"], 1)):
            if not np.isnan(obs["quadrants_radii"][irad, iquad]):
                nobs = nobs + 1
    wrad = obs["quadrants_speed"]

    # just for plotting
    xx = np.zeros((len(phi), len(r)))
    yy = np.zeros((len(phi), len(r)))
    # plt.pcolor(xx,yy,w)

    for j in range(len(phi)):
        for h in range(len(r)):
            xx[j, h] = 0.001 * r[h] * np.cos(phi[j] * np.pi / 180)
            yy[j, h] = 0.001 * r[h] * np.sin(phi[j] * np.pi / 180)

    if nobs > 0 and vmax > 21.0:
        # Do three fits
        for irep in range(3):
            # Default fit
            w = compute_wind_field(
                r,
                phi,
                vmax,
                pc,
                rmax,
                pn,
                vtreal,
                phit,
                lat,
                dpdt,
                phi_spiral,
                xn,
                vt,
                phia,
            )
            [mean_error, rms_error, err] = compute_mean_error(r, w, obs, wrad)
            nit = 0

            # 1) now first adjust size of storm with xn
            while True:
                nit = nit + 1
                wu = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn + dxn,
                    vt,
                    phia,
                )
                [mean_error_u, rms_error_u, err_u] = compute_mean_error(
                    r, wu, obs, wrad
                )
                wd = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn - dxn,
                    vt,
                    phia,
                )
                [mean_error_d, rms_error_d, err_d] = compute_mean_error(
                    r, wd, obs, wrad
                )

                # better with larger value of xn
                if rms_error_u < rms_error:
                    xn = xn + dxn
                    rms_error = rms_error_u

                # better with smaller value of xn
                elif rms_error_d < rms_error:
                    xn = xn - dxn
                    rms_error = rms_error_d

                # optimum reached
                else:
                    break

                # Other criteria
                if xn < 0.3 or xn > 1.0:
                    break
                if nit > 100:
                    break

            nit = 0

            # 2) now the asymmetry magnitude
            while True:
                nit = nit + 1
                wu = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn,
                    vt + dvt,
                    phia,
                )
                [mean_error_u, rms_error_u, err_u] = compute_mean_error(
                    r, wu, obs, wrad
                )
                wd = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn,
                    vt - dvt,
                    phia,
                )
                [mean_error_d, rms_error_d, err_d] = compute_mean_error(
                    r, wd, obs, wrad
                )

                # better with larger value of vt
                if rms_error_u < rms_error:
                    vt = vt + dvt
                    rms_error = rms_error_u

                # better with smaller value of vt
                elif rms_error_d < rms_error:
                    vt = vt - dvt
                    rms_error = rms_error_d

                # optimum reached
                else:
                    break

                # Other criteria
                if vt < 0.0 or vt > 1.5 * vtreal:
                    break
                if nit > 100:
                    break

            nit = 0

            # 3. now the asymmetry direction
            while True:
                nit = nit + 1

                wu = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn,
                    vt,
                    phia + dphia,
                )
                [mean_error_u, rms_error_u, err_u] = compute_mean_error(
                    r, wu, obs, wrad
                )
                wd = compute_wind_field(
                    r,
                    phi,
                    vmax,
                    pc,
                    rmax,
                    pn,
                    vtreal,
                    phit,
                    lat,
                    dpdt,
                    phi_spiral,
                    xn,
                    vt,
                    phia - dphia,
                )
                [mean_error_d, rms_error_d, err_d] = compute_mean_error(
                    r, wd, obs, wrad
                )

                if rms_error_u < rms_error:
                    # better with larger value of phia
                    phia = phia + dphia
                    rms_error = rms_error_u

                elif rms_error_d < rms_error:
                    # better with smaller value of phia
                    phia = phia - dphia
                    rms_error = rms_error_d
                else:
                    # optimum reached
                    break

                if phia < -180:
                    phia = phia + 360

                if phia > 180:
                    phia = phia - 360

                if nit > 100:
                    break

    return [xn, vt, phia]


def compute_mean_error(
    r: "np.ndarray",
    w: "np.ndarray",
    obs: dict,
    wrad: "np.ndarray",
) -> list:
    """
    Compute the mean and RMS error between modelled and observed wind radii.

    Parameters
    ----------
    r : numpy.ndarray
        Radii array (m).
    w : numpy.ndarray
        Modelled wind-speed array with shape ``(n_directions, n_radii)``.
    obs : dict
        Observed data; must contain ``"quadrants_radii"`` (ndarray shape (4, 4)).
    wrad : numpy.ndarray
        Target wind speeds for each radius threshold (m/s).

    Returns
    -------
    mean_error : float
        Mean signed error (m/s).
    rms_error : float
        Root-mean-square error (m/s).
    err : numpy.ndarray
        Per-quadrant/radius error array (NaN where no observation exists).
    """
    # variables
    nrad = np.size(obs["quadrants_radii"], 0)
    nrad = 3  # not we are only fitting R35, nothing else
    nq = np.size(obs["quadrants_radii"], 1)
    iq1 = [0, 9, 18, 27]
    iq2 = [9, 18, 27, 36]
    err = np.zeros((nq, nrad))

    # Go over the quadrants and radii
    for irad in range(nrad):
        for iquad in range(nq):
            vrad = 0
            for j in range(iq1[iquad], iq2[iquad] + 1):
                ww = w[j, :]
                if not np.isnan(obs["quadrants_radii"][irad, iquad]):
                    # compute wind speed vrad at required radius
                    wf = interp1d(r, ww, bounds_error=False, fill_value=0)
                    w0 = wf(obs["quadrants_radii"][irad, iquad] * 1000)
                    vrad = max(vrad, w0)
                else:
                    # maximum wind speed must be lower than wrad
                    w0 = max(ww)
                    vrad = max(vrad, w0)
            if not np.isnan(obs["quadrants_radii"][irad, iquad]):
                err[iquad, irad] = vrad - wrad[irad]
            else:
                err[iquad, irad] = np.nan

    # Get error values
    mask = ~np.isnan(err)  # Create the mask
    err = err[mask]
    mean_error = np.nanmean(err)
    rms_error = np.sqrt(np.mean(err**2))

    # Return
    return [mean_error, rms_error, err]


def compute_forward_speed_heading(
    t: "np.ndarray",
    x: "np.ndarray",
    y: "np.ndarray",
) -> list:
    """
    Compute forward speed and heading from a sequence of track positions.

    Uses forward differences at the first point, backward differences at the
    last point, and central differences elsewhere.

    Parameters
    ----------
    t : numpy.ndarray
        Time array (any unit; currently unused in the speed calculation).
    x : numpy.ndarray
        Longitude array (degrees).
    y : numpy.ndarray
        Latitude array (degrees).

    Returns
    -------
    forward_speed : numpy.ndarray
        Forward speed at each track point (m/s).
    heading : numpy.ndarray
        Track heading at each point (radians, Cartesian convention).
    """
    # variables
    forward_speed = np.zeros((len(x)))
    heading = np.zeros((len(x)))
    for it in range(len(x)):
        # Get basics
        geofacy = 111111
        geofacx = geofacy * np.cos(y[it] * np.pi / 180)

        if it == 0:
            # print('Forward')
            # datetime_forward    =
            # dt                  = datetime_forward - datetime_it
            # dt                  = dt.total_seconds()
            dx = (x[it + 1] - x[it]) * geofacx
            dy = (y[it + 1] - y[it]) * geofacy

        elif it == (len(x) - 1):
            # print('Backward')
            # datetime_backward   = datetime.strptime(self.track.datetime[it-1], dateformat_module)
            # coords_backward     = self.track.geometry[it-1]
            # dt                  = datetime_it - datetime_backward
            # dt                  = dt.total_seconds()
            dx = (x[it] - x[it - 1]) * geofacx
            dy = (y[it] - y[it - 1]) * geofacy

        else:
            # Backward
            dx1 = (x[it] - x[it - 1]) * geofacx
            dy1 = (y[it] - y[it - 1]) * geofacy

            # Forward
            dx2 = (x[it + 1] - x[it]) * geofacx
            dy2 = (y[it + 1] - y[it]) * geofacy

            # Combined yields central differences
            # print('Central')
            dx = np.mean([dx1, dx2])
            dy = np.mean([dy1, dy2])

        # Compute angle
        forward_speed[it] = 0.0
        heading[it] = np.arctan2(dy, dx)

    # Return
    return [forward_speed, heading]


def compute_wind_field(
    r: "np.ndarray",
    phi: "np.ndarray",
    vmax: float,
    pc: float,
    rmax: float,
    pn: float,
    vtreal: float,
    phit: float,
    lat: float,
    dpdt: float,
    phi_spiral: float,
    xn: float,
    vt: float,
    phia: float,
) -> "np.ndarray":
    """
    Compute a 2-D wind-speed field for a parametric tropical-cyclone model.

    Parameters
    ----------
    r : numpy.ndarray
        Radii (m).
    phi : numpy.ndarray
        Azimuth angles (degrees, Cartesian convention).
    vmax : float
        Maximum sustained wind speed (m/s).
    pc : float
        Central pressure (hPa).
    rmax : float
        Radius of maximum winds (m).
    pn : float
        Background pressure (hPa).
    vtreal : float
        Storm translation speed (m/s).
    phit : float
        Track heading angle (degrees).
    lat : float
        Storm latitude (degrees).
    dpdt : float
        Rate of pressure change (hPa/h).
    phi_spiral : float
        Inflow angle (degrees).
    xn : float
        Holland profile exponent.
    vt : float
        Asymmetry correction speed (m/s).
    phia : float
        Asymmetry direction (degrees).

    Returns
    -------
    wind_speed : numpy.ndarray
        2-D array of wind speed (m/s) with shape ``(n_phi, n_r)``.
    """
    vms = vmax - vt

    # compute wind profile (vr and pr)
    [vr, pr] = holland2010(r, vms, pc, pn, rmax, dpdt, lat, vtreal, xn)

    wind_speed = np.zeros((phi.shape[0], r.shape[0]))
    wind_to_direction_cart = np.zeros((phi.shape[0], r.shape[0]))
    for iphi in range(len(phi)):
        wind_speed[iphi, :] = vr
        if lat >= 0:
            # northern hemisphere
            dr = 90 + phi[iphi] + phi_spiral
        else:
            # southern hemisphere
            dr = -90 + phi[iphi] - phi_spiral
        wind_to_direction_cart[iphi, :] = dr

    vnorm = wind_speed / np.max(wind_speed)
    vnorm = np.zeros(np.shape(wind_speed)) + 1.0

    ux = vt * np.cos((phit + phia) * np.pi / 180)
    uy = vt * np.sin((phit + phia) * np.pi / 180)

    vx = wind_speed * np.cos(wind_to_direction_cart * np.pi / 180) + ux * vnorm
    vy = wind_speed * np.sin(wind_to_direction_cart * np.pi / 180) + uy * vnorm

    wind_speed = np.sqrt(vx**2 + vy**2)

    return wind_speed
