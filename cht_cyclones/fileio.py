"""
Alternative file I/O implementation for tropical-cyclone track files.

Provides :class:`TropicalCycloneTrack`, a class-based reader/writer that
supports the ``ddb_cyc`` and ``trk`` (COAMPS-TC) formats using explicit header
and data parsing routines.  This module co-exists with the functions in
:mod:`cht_cyclones.track`.
"""

from datetime import datetime
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point


class TropicalCycloneTrack:
    """
    Reader/writer for tropical-cyclone track files (``ddb_cyc`` and ``trk`` formats).

    The track data are stored in ``self.track``, a
    :class:`geopandas.GeoDataFrame` with one row per time step.

    Parameters
    ----------
    epsg : int, optional
        EPSG code for the track CRS (default 4326).
    debug : int, optional
        Verbosity level; ``1`` prints progress messages.
    """

    wind_profile: str
    wind_pressure_relation: str
    rmw_relation: str
    background_pressure: float
    phi_spiral: float
    wind_conversion_factor: float
    spiderweb_radius: float
    nr_radial_bins: int
    nr_directional_bins: int
    EPSG: int
    debug: int

    # Mapping of the header lines to the corresponding attribute names
    LINE_PREFIX_MAPPING_DDB_CYC = {
        "WindProfile": "wind_profile",
        "WindPressureRelation": "wind_pressure_relation",
        "RMaxRelation": "rmw_relation",
        "Backgroundpressure": "background_pressure",
        "PhiSpiral": "phi_spiral",
        "WindConversionFactor": "wind_conversion_factor",
        "SpiderwebRadius": "spiderweb_radius",
        "NrRadialBins": "nr_radial_bins",
        "NrDirectionalBins": "nr_directional_bins",
    }

    # List of field names corresponding to the extracted data values
    DDB_CYC_FIELD_NAMES = [
        "vmax",
        "pc",
        "RMW",
        "R35_NE",
        "R35_SE",
        "R35_SW",
        "R35_NW",
        "R50_NE",
        "R50_SE",
        "R50_SW",
        "R50_NW",
        "R65_NE",
        "R65_SE",
        "R65_SW",
        "R65_NW",
        "R100_NE",
        "R100_SE",
        "R100_SW",
        "R100_NW",
    ]

    def __init__(self, epsg: int = 4326, debug: int = 0) -> None:
        self.EPSG = epsg
        self.debug = debug

    def read(self, filename: Path, fmt: str) -> None:
        """
        Read a track file and populate ``self.track``.

        Parameters
        ----------
        filename : Path
            Path to the track file.
        fmt : str
            Format identifier: ``"ddb_cyc"`` or ``"trk"``.

        Raises
        ------
        ValueError
            If ``fmt`` is not one of the supported formats.
        """
        if fmt == "ddb_cyc":
            self._read_ddb_cyc(filename)
        elif fmt == "trk":
            self._read_trk(filename)
        else:
            raise ValueError("This file format is not supported as read track!")

    def write(self, filename: Path, fmt: str) -> None:
        """
        Write the track to a file.

        Parameters
        ----------
        filename : Path
            Output file path.
        fmt : str
            Format identifier.  Only ``"ddb_cyc"`` is currently supported.

        Raises
        ------
        ValueError
            If ``fmt`` is not a supported output format.
        """
        if fmt == "ddb_cyc":
            self._write_ddb_cyc(filename)
        else:
            raise ValueError(
                "Unsupported file format for writing the track! For other methods of writing the track; please used the 'tc.track.to_file' option"
            )

    def _read_ddb_cyc(self, filename: Path) -> None:
        # Read all the lines first
        with open(filename, "r") as f:
            lines = f.readlines()
        data_start = self._read_data_start(lines)

        # Read the name
        self.name = None
        for line in lines[:data_start]:
            if line.startswith("Name"):
                string_value = line[5:]
                string_value = "".join(ch for ch in string_value if ch.isalnum())
                self.name = string_value
                break
        if self.name is None:
            raise ValueError("Could not find the name of the tropical cyclone")

        self._read_header_ddb_cyc(lines[:data_start])
        self._read_data_ddb_cyc(lines[data_start:])

    def _read_trk(self, filename: Path) -> None:
        # Initialize variables
        lasttime = np.datetime64(datetime.strptime("2000010118", "%Y%m%d%H"))
        newtime2 = lasttime.astype("O")
        tc_time_string = newtime2.strftime("%Y%m%d %H%M%S")
        it = 0
        point = Point(0, 0)
        gdf = gpd.GeoDataFrame(
            {
                "datetime": [tc_time_string],
                "geometry": [point],
                "vmax": [-999],
                "pc": [-999],
                "RMW": [-999],
                "R35_NE": [-999],
                "R35_SE": [-999],
                "R35_SW": [-999],
                "R35_NW": [-999],
                "R50_NE": [-999],
                "R50_SE": [-999],
                "R50_SW": [-999],
                "R50_NW": [-999],
                "R65_NE": [-999],
                "R65_SE": [-999],
                "R65_SW": [-999],
                "R65_NW": [-999],
                "R100_NE": [-999],
                "R100_SE": [-999],
                "R100_SW": [-999],
                "R100_NW": [-999],
            }
        )
        gdf.set_crs(epsg=self.EPSG, inplace=True)
        self.track = gdf  # Always start clean

        # Open the file
        with open(filename, "r") as fid:
            # Read line by line
            for s0 in fid:
                s0 = s0.strip()
                if not s0:
                    continue

                # Start
                s = s0.split(",")
                tc = {"name": "not_known"}
                tc["basin"] = s[0]
                tc["storm_number"] = int(s[1])

                # Read time
                tstr = str(s[2])
                tstr = tstr.replace(" ", "")
                time = datetime.strptime(tstr, "%Y%m%d%H")
                newtime = np.datetime64(time)

                # Add forecasted time
                hrs = float(s[5])
                newtime = newtime + np.timedelta64(int(hrs * 3600), "s")
                if newtime > lasttime + np.timedelta64(1, "s"):
                    # Update all variables in gdf (when it is not the first one)
                    if it > 0:
                        gdf["geometry"] = [Point(x, y)]  # Assign the new geometry
                        gdf.vmax = vmax
                        gdf.pc = pc
                        gdf.RMW = rmax
                        if R35_NE > 0:
                            gdf.R35_NE = R35_NE
                        if R35_SE > 0:
                            gdf.R35_SE = R35_SE
                        if R35_SW > 0:
                            gdf.R35_SW = R35_SW
                        if R35_NW > 0:
                            gdf.R35_NW = R35_NW

                        if R50_NE > 0:
                            gdf.R50_NE = R50_NE
                        if R50_SE > 0:
                            gdf.R50_SE = R50_SE
                        if R50_SW > 0:
                            gdf.R50_SW = R50_SW
                        if R50_NW > 0:
                            gdf.R50_NW = R50_NW

                        if R65_NE > 0:
                            gdf.R65_NE = R65_NE
                        if R65_SE > 0:
                            gdf.R65_SE = R65_SE
                        if R65_SW > 0:
                            gdf.R65_SW = R65_SW
                        if R65_NW > 0:
                            gdf.R65_NW = R65_NW

                        if R100_NE > 0:
                            gdf.R100_NE = R100_NE
                        if R100_SE > 0:
                            gdf.R100_SE = R100_SE
                        if R100_SW > 0:
                            gdf.R100_SW = R100_SW
                        if R100_NW > 0:
                            gdf.R100_NW = R100_NW

                    # New time point found
                    it += 1
                    lasttime = newtime

                    # Append self
                    self.track = pd.concat([self.track, gdf])

                    # Make TC string
                    newtime2 = newtime.astype("O")
                    tc_time_string = newtime2.strftime("%Y%m%d %H%M%S")

                    # Start fresh
                    gdf = gpd.GeoDataFrame(
                        {
                            "datetime": [tc_time_string],
                            "geometry": [point],
                            "vmax": [-999],
                            "pc": [-999],
                            "RMW": [-999],
                            "R35_NE": [-999],
                            "R35_SE": [-999],
                            "R35_SW": [-999],
                            "R35_NW": [-999],
                            "R50_NE": [-999],
                            "R50_SE": [-999],
                            "R50_SW": [-999],
                            "R50_NW": [-999],
                            "R65_NE": [-999],
                            "R65_SE": [-999],
                            "R65_SW": [-999],
                            "R65_NW": [-999],
                            "R100_NE": [-999],
                            "R100_SE": [-999],
                            "R100_SW": [-999],
                            "R100_NW": [-999],
                        }
                    )
                    gdf.set_crs(epsg=self.EPSG, inplace=True)

                    # Reset all values
                    x = 0
                    y = 0
                    vmax = -999
                    pc = -999
                    RMW = -999
                    R35_NE = -999
                    R35_SE = -999
                    R35_SW = -999
                    R35_NW = -999

                    R50_NE = -999
                    R50_SE = -999
                    R50_SW = -999
                    R50_NW = -999

                    R65_NE = -999
                    R65_SE = -999
                    R65_SW = -999
                    R65_NW = -999

                    R100_NE = -999
                    R100_SE = -999
                    R100_SW = -999
                    R100_NW = -999

                # Latitude
                if s[6][-1] == "N":
                    y = 0.1 * float(s[6][:-1])
                else:
                    y = -0.1 * float(s[6][:-1])
                # Longitude
                if s[7][-1] == "E":
                    x = 0.1 * float(s[7][:-1])
                else:
                    x = -0.1 * float(s[7][:-1])

                # vmax
                if float(s[8]) > 0:
                    vmax = float(s[8])

                # Pressure and raddi
                if len(s) > 9:
                    # Pressure
                    if float(s[9]) > 0:
                        pc = float(s[9])

                    # Radii
                    r = float(s[11])
                    if r in [34, 35]:
                        R35_NE = float(s[13])
                        R35_SE = float(s[14])
                        R35_SW = float(s[15])
                        R35_NW = float(s[16])
                    elif r == 50:
                        R50_NE = float(s[13])
                        R50_SE = float(s[14])
                        R50_SW = float(s[15])
                        R50_NW = float(s[16])
                    elif r in [64, 65]:
                        R65_NE = float(s[13])
                        R65_SE = float(s[14])
                        R65_SW = float(s[15])
                        R65_NW = float(s[16])
                    elif r == 100:
                        R100_NE = float(s[13])
                        R100_SE = float(s[14])
                        R100_SW = float(s[15])
                        R100_NW = float(s[16])

                    # Other things
                    if len(s) > 17:
                        if s[17]:
                            pressure_last_closed_isobar = float(s[17])
                        if s[18]:
                            radius_last_closed_isobar = float(s[18])
                        if s[19]:
                            try:
                                if float(s[19]) > 0:
                                    rmax = float(s[19])
                            except ValueError:
                                print("Error: Unable to convert to float for RMW")
                                rmax = -999
                        if len(s) >= 28:
                            if s[27]:
                                name = s[27]

        # Done with this => rest so track looks good
        self.track = self.track.reset_index(drop=True)
        self.track = self.track.drop([0])  # remove the dummy
        self.track = self.track.reset_index(drop=True)
        self.track = self.track.drop([0])  # remove the dummy
        self.track = self.track.reset_index(drop=True)
        if self.debug == 1:
            print("Successfully read track - trk")

    def _read_data_start(self, lines: list[str]) -> int:
        data_start = None
        for i, line in enumerate(lines):
            if line.startswith("##    Datetime ") or line.startswith("#   Date   Time"):
                data_start = i
                break

        if data_start is None:
            raise ValueError("Could not find the start of the data in the file")
        return data_start

    def _read_header_ddb_cyc(self, header_lines: list[str]) -> None:
        for line in header_lines:
            for prefix, attribute in self.LINE_PREFIX_MAPPING_DDB_CYC.items():
                if line.startswith(prefix):
                    string_value = self._extract_header_value(
                        line.strip(), len(prefix) + 1
                    )

                    # Check if the attribute requires a float or integer conversion
                    if attribute in [
                        "background_pressure",
                        "phi_spiral",
                        "wind_conversion_factor",
                        "spiderweb_radius",
                    ]:
                        setattr(self, attribute, float(string_value))
                    elif attribute in ["nr_radial_bins", "nr_directional_bins"]:
                        setattr(self, attribute, int(string_value))
                    else:
                        setattr(self, attribute, string_value)
                    break  # Exit the loop once a match is found

    def _read_data_ddb_cyc(self, data_lines: list[str]) -> None:
        # Read data
        for line in data_lines:
            if line.strip().startswith("#"):
                continue
            # Extract the date and values
            tc_time, values = self._extract_data_values_from_line(line)
            tc_time_string = tc_time.strftime("%Y%m%d %H%M%S")

            # Create a dictionary for the fields
            if len(values[2:]) != len(self.DDB_CYC_FIELD_NAMES):
                raise ValueError(
                    f"Expected {len(self.DDB_CYC_FIELD_NAMES)} values, but got {len(values[2:])}"
                )
            track_data = dict(
                zip(
                    ["datetime"] + self.DDB_CYC_FIELD_NAMES,
                    [tc_time_string] + values[2:],
                    strict=True,
                )
            )

            # Create the GeoDataFrame
            point = Point(
                values[0], values[1]
            )  # Assuming x = values[0] and y = values[1]
            gdf = gpd.GeoDataFrame({**track_data, "geometry": [point]}, index=[0])
            gdf.set_crs(epsg=self.EPSG, inplace=True)

            # Append the GeoDataFrame to the track
            if hasattr(self, "track"):
                self.track = pd.concat([self.track, gdf], ignore_index=True)
            else:
                self.track = gdf

        # Done with this
        self.track = self.track.reset_index(drop=True)
        self.track = self.track.drop([0])  # remove the dummy
        self.track = self.track.reset_index(drop=True)
        if self.debug == 1:
            print("Successfully read track - ddb_cyc")

    def _write_ddb_cyc(self, filename: Path) -> None:
        """
        Write the track to a Delft Dashboard ``ddb_cyc`` format file.

        Parameters
        ----------
        filename : Path
            Output file path.
        """
        with open(filename, "wt") as f:
            # Print header
            f.writelines(
                f"# Tropical Cyclone Toolbox - Coastal Hazards Toolkit - "
                f"{self.creation_time.strftime(dateformat_module)}\n"
            )

            # Print rest
            f.writelines(f'Name                   "{self.name}"\n')
            f.writelines(f"WindProfile            {self.wind_profile}\n")
            f.writelines(f"WindPressureRelation   {self.wind_pressure_relation}\n")
            f.writelines(f"RMaxRelation           {self.rmw_relation}\n")
            f.writelines(f"Backgroundpressure     {self.background_pressure}\n")
            f.writelines(f"PhiSpiral              {self.phi_spiral}\n")
            f.writelines(f"WindConversionFactor   {self.wind_conversion_factor}\n")
            f.writelines(f"SpiderwebRadius        {self.spiderweb_radius}\n")
            f.writelines(f"NrRadialBins           {self.nr_radial_bins}\n")
            f.writelines(f"NrDirectionalBins      {self.nr_directional_bins}\n")
            epsg = self.track.crs.name
            f.writelines(f"EPSG                   {epsg}\n")
            f.writelines(f"UnitIntensity          {self.unit_intensity}\n")
            f.writelines(f"UnitWindRadii          {self.unit_radii}\n")

            # Print header for the track
            f.writelines("#  \n")
            f.writelines(
                "#   Date   Time               Lat        Lon         Vmax       Pc          Rmax         R35(NE)      R35(SE)     R35(SW)     R35(NW)     R50(NE)     R50(SE)    R50(SW)    R50(NW)     R65(NE)     R65(SE)     R65(SW)     R65(NW)    R100(NE)    R100(SE)    R100(SW)    R100(NE)  \n"
            )
            f.writelines("#  \n")

            # Print the actual track
            for i in range(len(self.track)):
                f.writelines(self.track.datetime[i].rjust(20))
                coords = self.track.geometry[i]
                f.writelines(str(round(coords.y, 2)).rjust(12))
                f.writelines(str(round(coords.x, 2)).rjust(12))

                f.writelines(str(round(self.track.vmax[i], 1)).rjust(12))
                f.writelines(str(round(self.track.pc[i], 1)).rjust(12))
                f.writelines(str(round(self.track.RMW[i], 1)).rjust(12))

                f.writelines(str(round(self.track.R35_NE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R35_SE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R35_SW[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R35_NW[i], 1)).rjust(12))

                f.writelines(str(round(self.track.R50_NE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R50_SE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R50_SW[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R50_NW[i], 1)).rjust(12))

                f.writelines(str(round(self.track.R65_NE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R65_SE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R65_SW[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R65_NW[i], 1)).rjust(12))

                f.writelines(str(round(self.track.R100_NE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R100_SE[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R100_SW[i], 1)).rjust(12))
                f.writelines(str(round(self.track.R100_NW[i], 1)).rjust(12))

                f.writelines("\n")

        if self.debug == 1:
            print("Successfully written track - ddb_cyc")

    @staticmethod
    def _extract_header_value(line: str, start_index: int) -> str:
        string_value = line[start_index:].strip()
        return "".join(ch for ch in string_value if ch.isalnum())

    @staticmethod
    def _extract_data_values_from_line(line: str) -> tuple[datetime, list[float]]:
        values = line.strip().split()
        date_format = "%Y%m%d %H%M%S"
        date_string = f"{values[0]} {values[1]}"
        tc_time = datetime.strptime(date_string, date_format)
        values_float = [
            float(v) for v in values[2:]
        ]  # Extract float values from the line
        return tc_time, values_float
