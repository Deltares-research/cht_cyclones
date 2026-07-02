"""
NHC Forecast/Advisory downloader and organiser.

Downloads the latest National Hurricane Center Forecast/Advisory (TCM) products
for all active storms, organises them by year/basin/storm, and maintains a set
of merged per-storm ``.cyc`` track files.  Also provides
:func:`find_nhc_track_file` for locating a storm track that intersects a given
model forecast area and time window.

This is the NHC counterpart to :mod:`cht_cyclones.jtwc.jtwc`; the two share the
same on-disk archive layout and the same public helpers.
"""

import os
import shutil
from datetime import datetime

import requests
from bs4 import BeautifulSoup

from cht_cyclones import TropicalCyclone, TropicalCycloneTrack

# Machine-readable index of all currently active NHC storms.
CURRENT_STORMS_URL = "https://www.nhc.noaa.gov/CurrentStorms.json"
# ATCF official-forecast files (fst/<basin><cy><year>.fst).
ATCF_FST_URL = "https://ftp.nhc.noaa.gov/atcf/fst/"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/115.0.0.0 Safari/537.36"
    )
}


def download(path: str, source: str = "tcm") -> None:
    """
    Download the latest NHC forecasts and organise them under ``path``.

    Parameters
    ----------
    path : str
        Root directory for the NHC track archive.
    source : str, optional
        Which NHC product to download:

        * ``"tcm"`` (default) — the Forecast/Advisory (TCM) text product.
        * ``"atcf"`` — the ATCF official-forecast file (``.fst``), a
          comma-delimited table with 34/50/64 kt quadrant wind radii.
    """
    if source == "tcm":
        download_forecast_advisories(os.path.join(path, "_downloads"))
    elif source == "atcf":
        download_fst_files(os.path.join(path, "_downloads"))
    else:
        raise ValueError(f"Unknown source '{source}' (expected 'tcm' or 'atcf')")
    organize(path)


def _extract_advisory_text(url: str) -> str:
    """
    Fetch a Forecast/Advisory ``.shtml`` page and return the raw product text.

    The NHC serves the plain-text product wrapped in a ``<pre>`` element inside
    an HTML page; this strips the HTML and returns the enclosed text.

    Parameters
    ----------
    url : str
        URL of the Forecast/Advisory ``.shtml`` page.

    Returns
    -------
    str
        The raw Forecast/Advisory text.

    Raises
    ------
    requests.HTTPError
        If the HTTP request returns a non-200 status code.
    RuntimeError
        If no ``<pre>`` block is found in the page.
    """
    r = requests.get(url, headers=HEADERS)
    r.raise_for_status()
    soup = BeautifulSoup(r.content, "html.parser")
    pre = soup.find("pre")
    if pre is None:
        raise RuntimeError(f"No advisory text (<pre>) found at {url}")
    return pre.get_text()


def download_forecast_advisories(path: str) -> None:
    """
    Download the Forecast/Advisory of every active NHC storm into ``path``.

    Each product is saved as ``<storm_id>_adv<NNN>.txt`` (e.g.
    ``al092024_adv015.txt``), where the storm id encodes basin, storm number and
    year.

    Parameters
    ----------
    path : str
        Destination directory for downloaded ``.txt`` products.

    Raises
    ------
    requests.HTTPError
        If the current-storms index cannot be retrieved.
    """
    # If path does not exist, create it
    if not os.path.exists(path):
        os.makedirs(path)

    # Delete all files in path
    for f in os.listdir(path):
        os.remove(os.path.join(path, f))

    r = requests.get(CURRENT_STORMS_URL, headers=HEADERS)
    r.raise_for_status()
    data = r.json()

    storms = data.get("activeStorms") or []
    if not storms:
        print("No active NHC storms found.")
        return

    for storm in storms:
        storm_id = storm.get("id")
        fcst = storm.get("forecastAdvisory")
        if not storm_id or not fcst or not fcst.get("url"):
            continue

        adv_num = fcst.get("advNum", "xxx")
        url = fcst["url"]
        print(f"Found NHC Forecast/Advisory: {storm_id} adv {adv_num} ({url})")

        try:
            text = _extract_advisory_text(url)
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            continue

        filename = os.path.join(path, f"{storm_id}_adv{adv_num}.txt")
        with open(filename, "w") as f:
            f.write(text)
        print(f"Downloaded: {filename}")


def download_fst_files(path: str) -> None:
    """
    Download the ATCF official-forecast ``.fst`` file of every active storm.

    Each file is saved as ``<storm_id>.fst`` (e.g. ``al092024.fst``), where the
    storm id encodes basin, storm number and year.

    Parameters
    ----------
    path : str
        Destination directory for downloaded ``.fst`` files.

    Raises
    ------
    requests.HTTPError
        If the current-storms index cannot be retrieved.
    """
    # If path does not exist, create it
    if not os.path.exists(path):
        os.makedirs(path)

    # Delete all files in path
    for f in os.listdir(path):
        os.remove(os.path.join(path, f))

    r = requests.get(CURRENT_STORMS_URL, headers=HEADERS)
    r.raise_for_status()
    data = r.json()

    storms = data.get("activeStorms") or []
    if not storms:
        print("No active NHC storms found.")
        return

    for storm in storms:
        storm_id = storm.get("id")
        if not storm_id:
            continue

        url = f"{ATCF_FST_URL}{storm_id}.fst"
        print(f"Found NHC ATCF forecast: {storm_id} ({url})")

        try:
            resp = requests.get(url, headers=HEADERS)
            resp.raise_for_status()
        except Exception as e:
            print(f"Error downloading {url}: {e}")
            continue

        filename = os.path.join(path, f"{storm_id}.fst")
        with open(filename, "wb") as f:
            f.write(resp.content)
        print(f"Downloaded: {filename}")


def _merge_storm_tracks(pth: str, basin: str, year: str, storm_num: str) -> None:
    """
    Rebuild the merged per-storm ``.cyc`` file from all individual cycle tracks.

    Parameters
    ----------
    pth : str
        Storm directory containing the individual ``.cyc`` files.
    basin, year, storm_num : str
        Storm identifiers used to name the merged file.
    """
    # Merge all individual tracks (chronological order). Exclude the previously
    # written merged file to avoid feeding it back in.
    cyc_files = sorted(
        os.path.join(pth, f)
        for f in os.listdir(pth)
        if f.endswith(".cyc") and not f.endswith("_merged.cyc")
    )
    tc = TropicalCyclone(track_file=cyc_files)
    tc.track.write(os.path.join(pth, f"{basin}_{year}_{storm_num}_merged.cyc"))


def organize(nhc_path: str) -> None:
    """
    Organise downloaded NHC forecast files into the archive structure.

    Handles both Forecast/Advisory (TCM) text files (``.txt``) and ATCF
    official-forecast files (``.fst``).  For each file: copies the raw file to a
    format-specific sub-folder, writes a ``.cyc`` track file tagged by advisory
    number (TCM) or synoptic init time (ATCF), and rebuilds the merged per-storm
    track file.

    Parameters
    ----------
    nhc_path : str
        Root directory of the NHC track archive.
    """
    download_path = os.path.join(nhc_path, "_downloads")

    # Loop through all downloaded files
    for filename in os.listdir(download_path):
        # Storm id is the leading token of the filename (e.g. "al092024").
        storm_id = filename.split("_")[0].split(".")[0]
        if len(storm_id) < 8:
            continue
        basin = storm_id[0:2].upper()
        storm_num = storm_id[2:4]
        year = storm_id[4:8]

        if filename.endswith(".txt"):
            # Forecast/Advisory (TCM) text product.
            track = TropicalCycloneTrack()
            config, name, advisory = track.read(
                os.path.join(download_path, filename), format="nhc_forecast_advisory"
            )
            tag = f"adv{advisory:03d}" if advisory is not None else "advxxx"
            raw_subfolder = "forecast_advisory"
            raw_name = f"{storm_id}_{tag}.txt"

        elif filename.endswith(".fst"):
            # ATCF official-forecast file. Tag by synoptic init time (= tau 0).
            track = TropicalCycloneTrack()
            config, name, advisory = track.read(
                os.path.join(download_path, filename), format="atcf"
            )
            init = datetime.strptime(track.gdf.datetime[0], "%Y%m%d %H%M%S")
            tag = init.strftime("%Y%m%d%H")
            raw_subfolder = "atcf"
            raw_name = f"{storm_id}_{tag}.fst"

        else:
            continue

        # Copy raw data file to the format-specific sub-folder.
        pth = os.path.join(nhc_path, year, basin, storm_num, raw_subfolder)
        os.makedirs(pth, exist_ok=True)
        shutil.copy(os.path.join(download_path, filename), os.path.join(pth, raw_name))

        # Write the *.cyc file for this cycle.
        pth = os.path.join(nhc_path, year, basin, storm_num)
        os.makedirs(pth, exist_ok=True)
        track.write(os.path.join(pth, f"{basin}_{year}_{storm_num}_{tag}.cyc"))

        # Rebuild the merged per-storm track.
        _merge_storm_tracks(pth, basin, year, storm_num)


def find_nhc_track_file(
    nhc_path: str,
    t0,
    t1,
    forecast_area,
) -> tuple:
    """
    Search the NHC archive for a storm that overlaps with a forecast domain.

    Walks the archive for the year of ``t0``, checks each merged track file,
    and returns the first one whose track points fall within ``forecast_area``
    during the ``[t0, t1]`` window.

    Parameters
    ----------
    nhc_path : str
        Root directory of the NHC track archive.
    t0 : datetime
        Start of the forecast window.
    t1 : datetime
        End of the forecast window.
    forecast_area : shapely.geometry.Polygon or None
        Spatial domain to test against; if ``None`` only the time window is used.

    Returns
    -------
    track_file_name : str or None
        Path to the merged ``.cyc`` file, or ``None`` if no match was found.
    storm_name : str or None
        Storm identifier string, or ``None`` if no match was found.
    """
    storm_name = None
    track_file_name = None

    # Get the year of the cycle
    year = t0.strftime("%Y")
    nhc_yr_path = os.path.join(nhc_path, year)
    # Read in every merged track for this year, limit it to the start/stop time
    # of the cycle, and use it if it has data in the window and overlaps extents.
    for root, dirs, files in os.walk(nhc_yr_path):
        for file in files:
            if file.endswith("_merged.cyc"):
                # Read the track file
                tc = TropicalCyclone(track_file=os.path.join(root, file))
                tc.track.shorten(tstart=t0, tend=t1)
                gdf = tc.track.gdf
                # Check if the track has data in the time window
                if len(gdf) > 0 and forecast_area is not None:
                    # Check if any track point falls within the forecast area
                    for idx, row in gdf.iterrows():
                        if row["geometry"].within(forecast_area):
                            track_file_name = os.path.join(root, file)
                            # storm name is file without path and without _merged.cyc
                            storm_name = os.path.splitext(os.path.basename(file))[
                                0
                            ].replace("_merged", "")
                            break
    return track_file_name, storm_name
