"""
Parametric wind-profile and wind-pressure-relation functions used by cht_cyclones.

References
----------
Holland, G. J., Belanger, J. I., & Fritz, A. (2010).
    A revised model for radial profiles of hurricane winds.
    Monthly Weather Review, 138(12), 4393-4401.

Holland, G. J. (2008).
    A revised hurricane pressure–wind model.
    Monthly Weather Review, 136(9), 3432-3445.

Nederhoff, K., et al. (2019).
    Estimates of tropical cyclone geometry parameters based on best-track data.
    Natural Hazards and Earth System Sciences, 19(11), 2359-2370.
"""

import numpy as np


def holland2010(
    r: "np.ndarray",
    vmax: float,
    pc: float,
    pn: float,
    rmax: float,
    dpdt: float,
    lat: float,
    vt: float,
    xn: float,
) -> list:
    """
    One-dimensional Holland et al. (2010) parametric wind profile.

    Parameters
    ----------
    r : numpy.ndarray
        Radii at which to evaluate the profile (km).
    vmax : float
        Maximum sustained wind speed (m/s).
    pc : float
        Central pressure (hPa).
    pn : float
        Background (ambient) pressure (hPa).
    rmax : float
        Radius of maximum winds (m).
    dpdt : float
        Rate of pressure change (hPa/h).
    lat : float
        Latitude of the storm centre (degrees).
    vt : float
        Storm translation speed (m/s).
    xn : float
        Holland profile exponent (dimensionless, typically 0.5).

    Returns
    -------
    vr : numpy.ndarray
        Radial wind speed profile (m/s).
    pr : numpy.ndarray
        Radial pressure profile (hPa).
    """
    # calculate Holland b parameter based pn Holland (2008) - assume Dvorak method
    vms = vmax
    dp = max(pn - pc, 1)
    x = 0.6 * (1 - dp / 215)
    if np.isnan(dpdt):
        dpdt = 0
    b = (
        -4.4 * 10**-5 * dp**2
        + 0.01 * dp
        + 0.03 * dpdt
        - 0.014 * abs(lat)
        + 0.15 * vt**x
        + 1
    )
    b = min(np.nanmax(np.append(b, 1)), 2)  # numerical limits

    # initialise
    x = np.zeros(np.size(r)) + xn
    pr = np.zeros((np.size(r)))
    vr = np.zeros((np.size(r)))

    # Compute
    index = r <= rmax
    x[index] = 0.5
    rp = (rmax / r) ** b
    pr = pc + dp * np.exp(-rp)
    vr = vms * (rp * np.exp(1.0 - rp)) ** x

    # output
    return [vr, pr]


def wind_radii_nederhoff(
    vmax: float,
    lat: float,
    region: int,
    probability: int,
) -> list:
    """
    Estimate the radius of maximum winds (RMW) and gale-force wind radius (R35).

    Based on the statistical relationships from Nederhoff et al. (2019).

    Parameters
    ----------
    vmax : float
        Maximum sustained wind speed (m/s, 1-minute average).
    lat : float
        Latitude of the storm centre (degrees).
    region : int
        Region index (0–7); use 7 for a globally applicable estimate.
    probability : int
        ``0`` returns only the mode; ``1`` also draws 1 000 random samples.

    Returns
    -------
    rmax : dict
        Estimated RMW statistics: keys ``"mode"``, and optionally ``"mean"``,
        ``"median"``, ``"lowest"``, ``"highest"``, ``"numbers"``.
    dr35 : dict
        Estimated *delta* R35 statistics (same keys as ``rmax``).
    """
    # radius of maximum winds (rmw or rmax)
    # 1. coefficients for A
    coefficients_a = np.array(
        [
            0.306982540000000,
            0.338409237000000,
            0.342791450000000,
            0.363546490000000,
            0.358572938000000,
            0.310729085000000,
            0.395431764000000,
            0.370190027000000,
        ]
    )

    # 2. coefficients for B
    coefficients_b = np.array(
        [
            [
                132.411906200000,
                14.5640379700000,
                -0.00259703300000000,
                20.3808036500000,
            ],
            [229.245844100000, 9.53865069100000, 0.00398810500000000, 28.4457367200000],
            [85.2576655100000, 30.6920872600000, 0.00243248000000000, 5.78116540600000],
            [127.833300700000, 11.8474757400000, 0.0159363120000000, 25.4682000500000],
            [153.733294700000, 11.4788885400000, 0.00747119300000000, 28.9489788700000],
            [261.528874200000, 7.01151785400000, 0.0261912560000000, 29.2022787100000],
            [19.0899242800000, 24.0885573100000, 0.106240340000000, 23.1802014600000],
            [44.8241743300000, 23.3717128800000, 0.0304690570000000, 22.4282036100000],
        ]
    )

    # 3. get the best guess for a and b given wind speed and latitude
    a_value = coefficients_a[region]
    b_value = (
        coefficients_b[region, 0]
        * np.exp(-vmax / coefficients_b[region, 1])
        * (1 + coefficients_b[region, 2] * abs(lat))
        + coefficients_b[region, 3]
    )

    rmax = {}
    rmax["mode"] = {}
    rmax["mean"] = {}
    rmax["median"] = {}
    rmax["lowest"] = {}
    rmax["highest"] = {}
    rmax["numbers"] = {}

    # 4. compute 1000 delta r35 values
    rmax["mode"] = np.exp(np.log(b_value) - a_value**2)
    if probability == 1:
        numbers = np.sort(
            np.exp(np.random.normal(size=(1000, 1)) * a_value + np.log(b_value))
        )
        rmax["mean"] = np.mean(numbers)
        rmax["median"] = np.median(numbers)
        rmax["lowest"] = numbers[int(0.05 * len(numbers))][0]
        rmax["highest"] = numbers[int(0.95 * len(numbers))][0]
        rmax["numbers"] = np.sort(numbers)

    # delta radius of 35 knots (r35)
    dr35 = {}
    dr35["mode"] = {}
    dr35["mean"] = {}
    dr35["median"] = {}
    dr35["lowest"] = {}
    dr35["highest"] = {}
    dr35["numbers"] = {}

    # Only if wind speed is more than 20 m/s
    if vmax > 20:
        # 1. coefficients for a
        coefficients_a = np.array(
            [
                [0.121563729, -0.052184289, 0.032953813],
                [0.131188105, -0.044389473, 0.002253258],
                [0.122286754, -0.045355772, 0.013286154],
                [0.120490659, -0.035029431, -0.005249445],
                [0.156059522, -0.041685377, 0.004952978],
                [-0.251333213, -0.009072243, -0.00506365],
                [0.131903526, -0.042096876, 0.012443195],
                [0.190044585, -0.044602083, 0.006117124],
            ]
        )

        # 2. coefficients for b
        coefficients_b = np.array(
            [
                [30.92867473, 0.530681714, -0.012001645],
                [30.21210133, 0.414897465, 0.021689596],
                [26.58686237, 0.425916004, 0.028547278],
                [23.88007085, 0.43109144, 0.038119083],
                [33.26829485, 0.42859578, 0.017209431],
                [18.11013691, 0.486399912, 0.02955688],
                [16.9973011, 0.453713419, 0.054643743],
                [29.61141102, 0.4132484, 0.024418947],
            ]
        )

        # 3. get the best guess for a and b given wind speed and latitude
        a_value = coefficients_a[region, 0] + np.exp(
            vmax * coefficients_a[region, 1]
        ) * (1 + coefficients_a[region, 2] * abs(lat))
        b_value = (
            coefficients_b[region, 0]
            * (vmax - 18) ** coefficients_b[region, 1]
            * (1 + coefficients_b[region, 2] * abs(lat))
        )

        # 4. compute 1000 delta r35 values
        dr35["mode"] = np.exp(np.log(b_value) - a_value**2)
        if probability == 1:
            numbers = np.sort(
                np.exp(np.random.normal(size=(1000, 1)) * a_value + np.log(b_value))
            )
            dr35["mean"] = np.mean(numbers)
            dr35["median"] = np.median(numbers)
            dr35["lowest"] = numbers[int(0.05 * len(numbers))][0]
            dr35["highest"] = numbers[int(0.95 * len(numbers))][0]
            dr35["numbers"] = np.sort(numbers)

    # output
    return [rmax, dr35]


def wpr_holland2008(
    pc: float | None = None,
    pn: float | None = None,
    phi: float | None = None,
    vt: float | None = None,
    dpcdt: float | None = None,
    rhoa: float | None = None,
    SST: float | None = None,
    vmax: float | None = None,
) -> float:
    """
    Holland (2008) wind-pressure relation.

    Computes either the minimum central pressure (when ``vmax`` is given) or the
    maximum sustained wind speed (when ``vmax`` is ``None``).

    Parameters
    ----------
    pc : float, optional
        Central pressure (hPa).  Required when estimating ``vmax``.
    pn : float, optional
        Background pressure (hPa).
    phi : float, optional
        Latitude (degrees).
    vt : float, optional
        Storm translation speed (m/s).
    dpcdt : float, optional
        Rate of pressure change (hPa/h).
    rhoa : float, optional
        Air density (kg/m³).  Estimated from SST when not provided.
    SST : float, optional
        Sea surface temperature (°C).  Used to derive ``rhoa`` when ``rhoa``
        is not supplied.
    vmax : float, optional
        Maximum sustained wind speed (m/s).  When provided the function returns
        the central pressure; otherwise it returns ``vmax``.

    Returns
    -------
    float
        Central pressure (hPa) if ``vmax`` is given, otherwise maximum
        sustained wind speed (m/s).
    """
    # used when pc needs to be determined
    if not rhoa:
        if vmax:
            dp1 = np.arange(1, 151 + 5, 5)
            pc1 = pn - dp1
        else:
            dp1 = pn - pc
            pc1 = pc

        if not SST:
            Ts = 28.0 - 3 * (phi - 10) / 20  # surface temperature
        else:
            Ts = SST - 1

        prmw = pc1 + dp1 / 3.7
        qm = 0.9 * (3.802 / prmw) * np.exp(17.67 * Ts / (243.5 + Ts))  # vapor pressure
        Tvs = (Ts + 273.15) * (1.0 + 0.81 * qm)  # virtual surface air temperature
        Rspecific = 287.058
        rhoa = 100 * pc1 / (Rspecific * Tvs)

    # vmax to be determined
    if not vmax:
        pc = min(pc, pn - 1.0)
        dp = pn - pc
        x = 0.6 * (1 - dp / 215)
        bs = (
            -4.4e-5 * dp**2
            + 0.01 * dp
            + 0.03 * dpcdt
            - 0.014 * phi
            + 0.15 * vt**x
            + 1.0
        )
        vmax = np.sqrt(100 * bs * dp / (rhoa * np.e))
        output = vmax
    else:
        vtkmh = vt * 3.6
        vmaxkmh = vmax * 3.6 * 0.88
        dp = (
            0.00592
            * (1 - 0.0687 * vtkmh**0.33)
            * (1 + 0.00285 * abs(phi) ** 1.35)
            * vmaxkmh**1.81
        )
        output = pn - dp

    return output
