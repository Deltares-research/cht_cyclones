import numpy as np
import matplotlib.pyplot as plt

def phram_rainfall(Vmax, Rm, Re, shear_mag, shear_dir_deg, grid_size_km=200, dx_km=1):
    """
    Generate 2D rainfall field (mm/h) from simplified PHRaM parameters.

    Parameters:
    - Vmax: Maximum wind speed (m/s)
    - Rm: Radius of maximum rainfall (km)
    - Re: Exponential decay scale (km)
    - shear_mag: Vertical wind shear magnitude (m/s)
    - shear_dir_deg: Shear direction (degrees clockwise from North)
    - grid_size_km: Half-size of square grid (km)
    - dx_km: Grid resolution (km)

    Returns:
    - X, Y: 2D coordinate grids (km)
    - rainfall: 2D rainfall rate field (mm/h)
    """

    # Grid coordinates (centered on storm eye)
    x = np.arange(-grid_size_km, grid_size_km + dx_km, dx_km)
    y = np.arange(-grid_size_km, grid_size_km + dx_km, dx_km)
    X, Y = np.meshgrid(x, y)

    # Convert shear direction to radians (meteorological: degrees clockwise from North)
    shear_dir_rad = np.deg2rad(shear_dir_deg)

    # Calculate radius from storm center
    R = np.sqrt(X**2 + Y**2)

    # Azimuth relative to shear direction
    theta = np.arctan2(Y, X)  # radians CCW from East
    # Adjust theta so 0 corresponds to shear direction (clockwise from North)
    # Convert arctan2 range (-pi, pi) CCW from East to degrees CW from North:
    # East=0 rad CCW → 90 deg CW from North
    theta_deg = (450 - np.rad2deg(theta)) % 360
    azimuth = np.deg2rad((theta_deg - shear_dir_deg) % 360)  # relative azimuth to shear

    # Empirical regression for rainfall rates (simplified)
    # These coefficients are example placeholders from Lonfat et al. (2007)
    # You can replace them with refined regression if you want.
    Tm = 20 + 0.5 * (Vmax - 33)  # Peak rainfall rate at Rm (mm/h)
    T0 = 5                        # Rainfall rate at center (mm/h)

    # Symmetric rainfall profile
    P_sym = np.where(
        R < Rm,
        T0 + (Tm - T0) * (R / Rm),
        Tm * np.exp(-(R - Rm) / Re)
    )

    # Asymmetric rainfall components (Fourier n=1 and n=2)
    # Coefficients scale with shear magnitude (example values)
    a1 = 0.3 * shear_mag * np.exp(-R / (2 * Rm))
    b1 = 0.2 * shear_mag * np.exp(-R / (2 * Rm))
    a2 = 0.1 * shear_mag * np.exp(-R / (3 * Rm))
    b2 = 0.05 * shear_mag * np.exp(-R / (3 * Rm))

    P_asym = a1 * np.cos(azimuth) + b1 * np.sin(azimuth) + a2 * np.cos(2 * azimuth) + b2 * np.sin(2 * azimuth)

    rainfall = P_sym + P_asym

    # Ensure no negative rainfall
    rainfall = np.maximum(rainfall, 0)

    return X, Y, rainfall

# Example usage
if __name__ == "__main__":
    # Parameters for a typical hurricane
    Vmax = 60           # m/s (~100 knots)
    Rm = 20             # km
    Re = 60             # km
    shear_mag = 5       # m/s
    shear_dir_deg = 45  # degrees clockwise from North

    X, Y, rain = phram_rainfall(Vmax, Rm, Re, shear_mag, shear_dir_deg)

    plt.figure(figsize=(8, 6))
    cs = plt.contourf(X, Y, rain, levels=np.linspace(0, np.max(rain), 20), cmap="Blues")
    plt.colorbar(cs, label="Rainfall rate (mm/h)")
    plt.title("Simplified PHRaM Rainfall Field")
    plt.xlabel("X (km)")
    plt.ylabel("Y (km)")
    plt.axis("equal")
    plt.show()
