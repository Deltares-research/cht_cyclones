from datetime import datetime

import geopandas as gpd
from shapely.geometry import Point


class TropicalCycloneTrack:
    def __init__(self):
        pass

    def read(self, filename, fmt):
        
            # If ddb_cyc
        if fmt == "ddb_cyc":
            # Read all the lines first
            with open(filename, "r") as f:
                lines = f.readlines()

            # Define the name first
            for line in lines:
                if line[0:4] == "Name":
                    string_value = line[5:]
                    string_value = "".join(ch for ch in string_value if ch.isalnum())
                    self.name = string_value

            # Define other variables names (if they exist)
            for i in range(len(lines)):
                line = lines[i]
                if line[0:11] == "WindProfile":
                    string_value = line[23:]
                    string_value = "".join(ch for ch in string_value if ch.isalnum())
                    self.wind_profile = string_value
                if line[0:20] == "WindPressureRelation":
                    string_value = line[23:]
                    string_value = "".join(ch for ch in string_value if ch.isalnum())
                    self.wind_pressure_relation = string_value
                if line[0:12] == "RMaxRelation":
                    string_value = line[23:]
                    string_value = "".join(ch for ch in string_value if ch.isalnum())
                    self.rmw_relation = string_value
                if line[0:18] == "Backgroundpressure":
                    string_value = line[23:]
                    self.background_pressure = float(string_value)
                if line[0:9] == "PhiSpiral":
                    string_value = line[23:]
                    self.phi_spiral = float(string_value)
                if line[0:20] == "WindConversionFactor":
                    string_value = line[23:]
                    self.wind_conversion_factor = float(string_value)
                if line[0:15] == "SpiderwebRadius":
                    string_value = line[23:]
                    self.spiderweb_radius = float(string_value)
                if line[0:12] == "NrRadialBins":
                    string_value = line[23:]
                    self.nr_radial_bins = int(string_value)
                if line[0:17] == "NrDirectionalBins":
                    string_value = line[23:]
                    self.nr_directional_bins = int(string_value)

            # Read the track
            for i in range(len(lines)):
                line = lines[i]
                if line[0:15] == "##    Datetime " or line[0:15] == "#   Date   Time":
                    break

            # Place coordinates in Tropical Cyclone Track
            for j in range(i + 2, len(lines)):
                # Get values
                line = lines[j]
                line = line.split()
                date_format = "%Y%m%d %H%M%S"
                date_string = line[0] + " " + line[1]
                tc_time = datetime.strptime(date_string, date_format)
                tc_time_string = tc_time.strftime(date_format)
                y = float(line[2])
                x = float(line[3])
                vmax = float(line[4])
                pc = float(line[5])
                RMW = float(line[6])
                R35_NE = float(line[7])
                R35_SE = float(line[8])
                R35_SW = float(line[9])
                R35_NW = float(line[10])
                R50_NE = float(line[11])
                R50_SE = float(line[12])
                R50_SW = float(line[13])
                R50_NW = float(line[14])
                R65_NE = float(line[15])
                R65_SE = float(line[16])
                R65_SW = float(line[17])
                R65_NW = float(line[18])
                R100_NE = float(line[19])
                R100_SE = float(line[20])
                R100_SW = float(line[21])
                R100_NW = float(line[22])

                # Make GeoDataFrame
                point = Point(x, y)
                gdf = gpd.GeoDataFrame(
                    {
                        "datetime": [tc_time_string],
                        "geometry": [point],
                        "vmax": [vmax],
                        "pc": [pc],
                        "RMW": [RMW],
                        "R35_NE": [R35_NE],
                        "R35_SE": [R35_SE],
                        "R35_SW": [R35_SW],
                        "R35_NW": [R35_NW],
                        "R50_NE": [R50_NE],
                        "R50_SE": [R50_SE],
                        "R50_SW": [R50_SW],
                        "R50_NW": [R50_NW],
                        "R65_NE": [R65_NE],
                        "R65_SE": [R65_SE],
                        "R65_SW": [R65_SW],
                        "R65_NW": [R65_NW],
                        "R100_NE": [R100_NE],
                        "R100_SE": [R100_SE],
                        "R100_SW": [R100_SW],
                        "R100_NW": [R100_NW],
                    }
                )
                gdf.set_crs(epsg=self.EPSG, inplace=True)

                # Append self
                self.track = pd.concat([self.track, gdf])

            # Done with this
            self.track = self.track.reset_index(drop=True)
            self.track = self.track.drop([0])  # remove the dummy
            self.track = self.track.reset_index(drop=True)
            if self.debug == 1:
                print("Successfully read track - ddb_cyc")

        elif fmt == "trk":
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

        else:
            raise Exception("This file format is not supported as read track!")


    def write(self, filename, fmt):
        # If ddb_cyc
        if fmt == "ddb_cyc":
            # Open file
            with open(filename, "wt") as f:
                # Print header
                f.writelines(
                    "# Tropical Cyclone Toolbox - Coastal Hazards Toolkit - "
                    + self.creation_time.strftime(dateformat_module)
                    + "\n"
                )

                # Print rest
                f.writelines('Name                   "' + self.name + '"\n')
                f.writelines("WindProfile            " + self.wind_profile + "\n")
                f.writelines(
                    "WindPressureRelation   " + self.wind_pressure_relation + "\n"
                )
                f.writelines("RMaxRelation           " + self.rmw_relation + "\n")
                f.writelines(
                    "Backgroundpressure     " + str(self.background_pressure) + "\n"
                )
                f.writelines("PhiSpiral              " + str(self.phi_spiral) + "\n")
                f.writelines(
                    "WindConversionFactor   " + str(self.wind_conversion_factor) + "\n"
                )
                f.writelines(
                    "SpiderwebRadius        " + str(self.spiderweb_radius) + "\n"
                )
                f.writelines(
                    "NrRadialBins           " + str(self.nr_radial_bins) + "\n"
                )
                f.writelines(
                    "NrDirectionalBins      " + str(self.nr_directional_bins) + "\n"
                )
                epsg = self.track.crs.name
                f.writelines("EPSG                   " + epsg + "\n")
                f.writelines(
                    "UnitIntensity          " + str(self.unit_intensity) + "\n"
                )
                f.writelines("UnitWindRadii          " + str(self.unit_radii) + "\n")

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
        else:
            print(
                'For other methods of writing the track; please used the "tc.track.to_file" option'
            )
