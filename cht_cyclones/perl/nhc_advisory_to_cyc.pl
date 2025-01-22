#!/usr/bin/perl
#********************************************************************************
#  read_fcst_advisory.pl
#
#  This Perl script reads a NOAA National Weather Service (NWS)/Tropical Prediction
#  Center (TPC) marine advisory (AKA "marine warning" or "technical advisory") and
#  obtains date, time, location, central pressure (CP) and radius of maximum winds (RMW)
#  for use in the Holland hurricane wind model.  The current version (as of 18 Apr 2006)
#  reads from a local, user-supplied input file.  The output is written to a user-
#  supplied output file name, and optionally can be appended to an existing output
#  file.  The output format is consistent with that of a cyclone.dat file used by the
#  Holland hurricane wind model.
#
#  Portions of this software are taken and/or derived from the program 'hurricane'
#  which, as of Apr 2006, is available from http://www.solar.ifa.hawaii.edu/Tropical/Bin/
#  and is:
#
#   Copyright (C) 1996 Thomas R. Metcalf
#
#   This software is provided "as is" and is subject to change without
#   notice.  No warranty of any kind is made with regard to this software,
#   including, but not limited to, the implied warranties of
#   merchantability and fitness for a particular purpose.  The author shall
#   not be liable for any errors or for direct, indirect, special,
#   incidental or consequential damages in connection with the furnishing,
#   performance, or use of this software:  use it at your own risk.
#
#    This program is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public
#    License as published by the Free Software Foundation; either
#    version 2 of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#    You should have received a copy of the GNU General Public
#    License along with this library; if not, write to the Free
#    Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
#  Syntax: read_fcst_advisory.pl InFile num_rad_bins num_dir_bins tc_radius [btfile]
#  where:  InFile is a local input file (forecast file),
#          num_rad_bins is the number of Delft3D Spiderweb file radial bins,
#          num_dir_bins is the number of Delft3D Spiderweb file directional bins, and
#          tc_radius is the tropical cyclone radius (km)
#          btfile is the optional best track downloaded with the forecast
#
#  Revision History:
#  10 Apr 2006  Initial coding.  R.S. Linzell, PSI/Neptune Sciences Div.
#  14 Apr 2006  Added code from 'hurricane' by T.R. Metcalf (see Copyright above);
#               modified same to output needed data for Holland wind model.  (RSL)
#  18 Apr 2006  Put most of the processing into new sub. process_data; added
#               a valid data record counter in sub. process_data; added printing
#               of the number of valid data records to the main program; added
#               sub. get_user_input to obtain the user input; added output file
#               ($outfile) and append flag ($append) to command line argument
#               list; added sub. print_data to write the output data to the output
#               file; updated comments at beginning of file; modified sub. get_radius
#               to ignore zero radius values when calculating the average radius,
#               and to estimate the radius where only 1 forecast record has radii
#               instead of ignoring such records; modified def. of $FMT in sub.
#               print_data to print zero-padded hours in the time field; modified
#               sub. DUMP to estimate forecast central pressures by multiplying the
#               present CP by the ratio of the present max. wind speed to the max.
#               forecast wind speed.  (RSL)
#  27 Apr 2006  Added printing of header information to sub. print_data (if not appending
#               to the output file); added calculation of Holland "B" parameter to sub.
#               process_data, including header parameters, conversion factors and other
#               definitions required for calculating B; added "_O" to output record
#               to indicate present observation and "_F" to indicate forecast; added
#               code to sub. process_data to calculate an observed RMW (if not found)
#               from wind speed ratio and RMW of 1st valid forecast.  (RSL)
#  27 Jun 2006  Added code to check for invalid CP values (in subs. process_data & DUMP)
#               and defined minimum CP values for Atlantic & Pacific basins in sub.
#               process_data (cf. http://en.wikipedia.org/wiki/List_of_notable_tropical_cyclones);
#               added code to check for invalid RMW values (in sub. get_radius) &
#               defined max. RMW value in sub. process_data (cf. http://www.srh.weather.gov/
#               srh/jetstream/tropics/tc_structure.htm); fixed bug in sub. print_data
#               so that header info. is printed to an appended file if the file size
#               is zero (i.e., for the 1st advisory's output).  (RSL)
#  21 Jul 2006  Fixed bug in 3rd date algorithm in sub. process_data.  (RSL)
#  12 Sep 2006  Defined $MeanRMW in sub. DUMP; modified sub. print_data to set any
#               zero RMW value to the mean; changed ambient air pressure, $pn, in sub.
#               process_data to 1013.25 (std. atmosphere).  (RSL)
#  28 Jul 2011  Added sub. wind_pressure_relationship to better estimate minimum
#               central pressure given wind speed (Dvorak technique); replaced the def.
#               of $for_cp in sub. DUMP with a call to sub. wind_pressure_relationship.  (RSL)
#  01 Nov 2011  Began adding support for the Delft3D program, WES.EXE (i.e., added
#               creation of the *.inp and *.trk files for that program); moved the 'use' and
#               'require' statements to the top of the source code; .  (RSL)
#  02 Nov 2011  Upgraded sub. get_user_input to support Delft3D input parameters; updated
#               the help content in several places; upgraded sub. process_data to
#               calculate the storm movement direction & speed; upgraded sub. get_radius
#               to perform the conversion from R34/50/64 to R35/50/65/100 wind radii
#               if Delft3D output is generated.  (RSL)
#  03 Nov 2011  Fixed some problems in sub. get_radius; added an IF block to the main
#               program to call sub. print_data if NOT writing Delft3D files, or to
#               call subs. print_wes_inp_file, print_track_file if doing so.  (RSL)
#  07 Nov 2011  Added IF blocks to sub. print_data, main program to change longitude
#               from [0,360] to [-180,180] degrees if needed.  (RSL)
#  16 Dec 2011  Added code to the main program to build the TC storm name ($tc_name), and
#               added the storm name to the input parameter list of the call to sub.
#               print_track_file.  (RSL)
#  17 Feb 2012  Added a line to convert radius of maximum winds from km to nautical
#               miles & a line to replace a zero RMW with the $MeanRMW value prior to
#               calling print_track_file for WES.EXE output; fixed a bug in sub.
#               get_radius; added an IF block to sub. DUMP to include the analysis
#               wind radius in the calculation of $MeanRMW.  (RSL)
#  12 Mar 2012  Modified sub. print_data to ensure that longitudes are on [0,360]
#               as required by PCTides.  (RSL)
#  09 May 2012  Upgraded sub. process_data: added an array to store the 3-D R34/R50/R64
#               wind radii; added an array to store the date numbers for the analysis
#               and forecast periods, and added calls to Date::Calc::Mktime to load the
#               date number array.  (RSL)
#  10 May 2012  Moved the storm name code from the main program to sub. process_data;
#               added '$storm_id' and '$warning_number' parambers to sub. process_data;
#               upgraded sub. get_radius to add parameters & code to load the output 3-D
#               wind radii array; updated the call to sub. print_track_file to include
#               the environmental pressure ($pn) and to use the new 3-D wind radii
#               array (@r100_65_50_35_3D).  (RSL)
#  18 May 2012  Added '$warn_id' for the storm warning ID (based on the input file
#               name without the extension).  (RSL)
#  20 Aug 2012  Upgraded sub. print_data: added writing of the storm name to a file
#               named 'cyclone_name.txt' for use by an external plotting script.  (RSL)
#  05 Sep 2012  Changed all debug 'print' statements to 'warn'; upgraded sub. get_radius
#               to always use the ratio method to calculate the RMW value, and to set
#               the RMW to the (average) minimum of R34/R50/R64 - 5 nmi if it is greater
#               than or equal to that value, and changed the algorithm to use the
#               square of the ratio; updated a warning message in sub. DUMP; added
#               'last' statements in the IF block where appropriate in sub. wind_
#               pressure_relationship.  (RSL)
#  13 Sep 2012  Commented-out the call in sub. get_radius to sub. interp_rads and
#               the last FOR loop inside the last IF block for Delft3D output, so
#               that interpolation of the R34/... wind radii is not performed.  (RSL)
#  14 Sep 2012  Modified the IF statements in sub. get_radius within the FOR $k
#               block to include any wind radii of zero, and added a WARN statement
#               to notify the user that zero radii are included in the output.  (RSL)
#  18 Sep 2012  Updated comments & help content; enclosed the zero radii warning
#               in sub. get_radius in an IF block.  (RSL)
#  01 Nov 2012  Commented-out the IF block where the eye diameter, if present,
#               was used for the RMW, in sub. process_data.  (RSL)
#  11 Jan 2013  Fixed an issue in sub. process_data to prevent the storm name from
#               being mangled by text after the file header.  (RSL)
#  01 Feb 2013  Added code to support output in extended best track format; modified
#               the IF block in the main program to define the $wes_output flag
#               according to the name of the output file (.trk for WES.EXE output,
#               .bst for best track output); added IF blocks to the main program
#               to call newly added sub. write_best_track.pl instead of the
#               calls to print_track_file & print_wes_inp_file; added code to sub.
#               process_data to redefine $MISS for best track output; added code
#               to sub. get_radius to correctly initialize & load the @r100_65_50_35_3D
#               array for best track output.  (RSL)
#  04 Feb 2013  Added storage of the eye diameter in $eye_diameter to sub. process_data;
#               added $eye_diameter to the call to sub. write_best_track; changed the
#               'require' command to use 'best_track_routines.pl' vice 'write_best_track.pl'
#               to use the latest routine; changed the best track output routine call to
#               'write_ext_best_track()' to use the latest routine; added def. of
#               $EBT_MISS (missing value code), @rad_outer_iso (radius of the
#               outermost closed isobar) for the extended best track output option;
#               added @pres_outer_iso to store the background pressure in an array for
#               the best track output option; replaced $pn with @pres_outer_iso in the
#               call to sub. write_ext_best_track.  (RSL)
#  07 Feb 2013  Added new array @storm_id_strs to store the storm ID string & pass
#               to sub. write_ext_best_track for the extended best track option.
#               (RSL)
#  13 Feb 2013  Added new array @seye_diam to store the eye diameter value & pass
#               to sub. write_ext_best_track for the extended best track option.
#               (RSL)
#  22 Feb 2013  Added a nested IF block at the end of the main WHILE loop in sub.
#               process_data to check for & correct invalid (i.e., zero) forecast
#               pressure values; added an IF block to sub. get_radius to check for
#               & handle the case of an empty @spd array; added support for reading
#               and input best track file and prepending it to the extended best
#               track output.  (RSL)
#  08 Mar 2013  Added defs. $time_best_track (date number of current best track
#               record), $start_time_data (date number of 1st TC warning record),
#               & added an IF block to use only those best track records that
#               precede the 1st TC warning record (i.e., the analysis); added
#               array index ranges to the code that concatenates best track &
#               forecast/analysis arrays to correctly omit later best track
#               records.  (RSL)
#  12 Mar 2013  Added IF blocks within the best track code of the main program to
#               check the best track file type according to the file extension,
#               to call new sub. read_ext_best_track_file (in best_track_routines.pl),
#               and to provide support for NRL extended best track files.  (RSL)
#  15 May 2013  Added code to the main program to provide support for PCTides
#               output using a best track file, either with or without a forecast
#               /analysis file; upgraded sub. print_data to avoid printing if
#               no data are present.  (RSL)
#  17 May 2013  Fixed printing issues in sub. print_data to correctly handle the
#               case of only best track file input; added an IF block to the main
#               program to check whether the call to sub. read_std_best_track was
#               successful, and if not, to call sub. read_ext_best_track_file with
#               a flag of 1 to handle "standard" extended best track files.  (RSL)
#  05 Sep 2013  Added sub. storm_velocity, and upgraded sub. process_data to call
#               this newly added sub.  (RSL)
#  06 Sep 2013  Added code to main program to support best track file input for
#               Delft3D WES.EXE output format.  (RSL)
#
#********************************************************************************

#  Include file(s):
use File::Basename;
my $BaseDir = dirname($0);		     # Get name of directory in which this program resides
require "$BaseDir/interpl.pl";		 # External linear interpolation routine
require "$BaseDir/cyc_routines.pl";  # Extended best track format output routine
use POSIX qw(mktime);

our $DEBUG = 0;   # set to 1 to enable debug print statements
our $DEBUG1 = 0;   # set to 1 to enable debug print statements (in sub. get_radius)
$| = 1;		# set autoflush (of buffers)
my $numrecs = 0;

#  Define a flag for whether Delft3D WES.EXE files will be generated.
my $EBT_MISS = -99;      # Missing value code for optional extended best track output
my @rad_outer_iso = ();  # Radius of the outermost closed isobar (ext. best track option)
my @pres_outer_iso = (); # Pressure of the outermost closed isobar (ext. best track option)
my @storm_id_strs = ();  # Array of storm ID strings (ext. best track option)
my @eye_diam = ();       # Array of eye diameter values (ext. best track option)
my @dateb = ();          # Array of date strings for WES.EXE output with best track option
my $bt_input = 0;        # Flag for input best track file
my $READ_ALL = 0;        # Flag to read only complete input records

#  Similar parameters for the input best track file:
our (@latt, @long, @RMW, @MSLP, @datestr, @max_wind_kt, @storm_name, $envpress, @r100_65_50_35_3Db, @eye_diameter);
our (@R34, @R34_Rad, @R50, @R50_Rad, @R64, @R64_Rad, @data_type, @bt_datenum, @pres_outer_isobar, @rad_outer_isobar);
our (@for_yearb, @for_nmonthb, @for_dayb, @for_timeb, $st_name, $i34, $i50, $i64, $time_best_track, $start_time_data);
my ($sum34, $sum50, $sum64, $istop, $iebt, $ismissing) = (0) x 6;

#  Check the command line syntax and obtain user input.
my ($infile, $num_rad_bins, $num_dir_bins, $tc_radius, $btfile) = &get_user_input;
my $nm2km = 1.852; 		# Conversion for nmi to km

#  Check whether an input best track file was specified.
if (length($btfile) > 0) {
   #  This is the case, so set the flag to 1.
   $bt_input = 1;
   if ($DEBUG) {
      warn " BEST TRACK INPUT....\n";
   }
}

# Process the input data.
process_data($infile);

# Creating the cyclone file name; it will be put in the same directory as the forecast
# The format is 2 digit basin, 2 digit storm number, 2 digit year “_” storm name _ forecast time.cyc
my @tmp = split(" ", $tc_name); #Get just the storm name (i.e. remove "Hurricane" "Typoon" etc)
$sname = uc($tmp[$#tmp]); #Change name to upper case
$sid = substr(${storm_id},0,4); # This is the 2-dgit basin + storm number for NHC forecasts
$yr  = substr($year,2,2);
my $DIR = dirname(${infile}); # Get the directory where the forecast file resides
my $trkfile = "${DIR}/${sid}${yr}_${sname}_${year}${nmonth}${day}${time}00.cyc";

if ($nmonth < 10) {
   warn "Month# < 10; append 0 to the month number\n";
   $trkfile = "${DIR}/${sid}${yr}_${sname}_${year}0${nmonth}${day}${time}00.cyc";
}

if ($DEBUG) { warn "In MAIN: InFile: $infile  \n";
   warn " Trkfile = $trkfile, NumRadBins = $num_rad_bins, NumDirBins = $num_dir_bins, TC_Rad = $tc_radius.\n";
}

if ($DEBUG) {
   warn "In MAIN: TC_Dir = @tc_dir; TC_Speed = @tc_speed; For_Day = @for_day.\n";
   warn "TC_Dir has ", scalar @tc_dir, " elems.\n";
   warn "R100 array (1st group): $r100_65_50_35[0][0] $r100_65_50_35[0][1] $r100_65_50_35[0][2] $r100_65_50_35[0][3].\n";
   warn "Storm name is ", uc($tc_name), "; Storm ID is $storm_id; Warning ID is $warn_id\n";
   warn "Hour: @for_time\n";
}

#  Print any valid data records to the output file.
#  Prepend the initial values to the respective arrays.  Actually, the
#  1st value (i.e., index 0) is replaced, since for these arrays, the 1st valid
#  value is stored at index 1 (cf. sub. print_data).
$for_year[0] = $year;
$for_nmonth[0] = $nmonth;
$for_day[0] = $day;
$for_time[0] = $time;
$for_lat[0] = $latitude;
$for_lon[0] = $longitude;
$for_rad[0] = $radius;
$for_cp[0] = $cp;
$for_winds[0] = $winds;

#  Calculate the date number of the analysis record from the forecast/analysis file.
$start_time_data = mktime(0,0,$time,$day,$nmonth-1,$year-1900);

#warn "Start_Time_Date = $start_time_data for $year $nmonth $day $time:00:00.\n";


#  Combine the year, month, day for the date array.
my @date = ();
for my $i (0..$#for_year) {
   $date[$i] = sprintf("%4d%2.2d%2.2d", $for_year[$i], $for_nmonth[$i], $for_day[$i]);
   #  Append "00" minutes to time
   $for_time[$i] = sprintf("%2.2d00", $for_time[$i]);
   #  Change longitude from [0,360] to [-180,180] if needed.
   if ($for_lon[$i] > 180) { $for_lon[$i] -= 360; }
   if ($for_rad[$i] == 0) { $for_rad[$i] = $MeanRMW; }
   #  Convert RMW from km to nmi for non-PCTides output.
   $for_rad[$i] /= $nm2km;
   $pres_outer_iso[$i] = $pn;
   $storm_id_strs[$i] = $storm_id;
   if ($i == 0) {
      $eye_diam[$i] = $eye_diameter;
   } else {
      $eye_diam[$i] = $EBT_MISS;
   }
}

#  Check whether the @rad_outer_iso array has any values.
if (!@rad_outer_iso) {
   #  It has none, so load it with missing values.
   @rad_outer_iso = ($EBT_MISS) x scalar(@for_rad);   # Load this array
}

#  Check whether the @pres_outer_iso array has any values.
if (!@pres_outer_iso) {
   #  It has none, so load it with background pressure values.
   @pres_outer_iso = ($pn) x scalar(@for_rad);   # Load this array
}

#  Check the input best track file flag.
if ($bt_input == 1) {
   #  Read the input best track file, using the routine that depends on the file type.
   #  Standard JTWC or NHC ATCF best track file
   #  Check whether output is extended best track format; if so, then read all input best track records.
   $READ_ALL = 1;
   ($envpress,*latt,*long,*RMW,*MSLP,*datestr,*max_wind_kt,*storm_name,*R34,*R34_Rad,*R50,*R50_Rad,
      *R64,*R64_Rad,*bt_datenum,*pres_outer_isobar,*rad_outer_isobar,*storm_id,*eye_diameter) = read_std_best_track($READ_ALL,$btfile);

   #  Check whether the data read was successful.
   if (scalar(@datestr) < 1) {
      die "ERROR: Input best track file $btfile apparently is an unknown format; stopped";
   }
}
my @name = split ' ', $storm_name[$#storm_name];
$st_name = uc($name[$#name]);

#  Check for valid start_time_data, and correct if invalid (for PCTides format).
if ($start_time_data <= 0) {
   $ismissing = 1;
   #  Calculate the date number for this best track record, using the date/time parameters
   $start_time_data = mktime(0,0,substr($datestr[$#datestr], 8, 2),
      substr($datestr[$#datestr], 6, 2),substr($datestr[$#datestr], 4, 2)-1,substr($datestr[$#datestr], 0, 4)-1900);
}

#  Loop over best track records....
for my $i (0..$#latt) {
   #  Calculate the date number for this best track record, using the date/time parameters
   $time_best_track = mktime(0,0,substr($datestr[$i], 8, 2),
      substr($datestr[$i], 6, 2),substr($datestr[$i], 4, 2)-1,substr($datestr[$i], 0, 4)-1900);
   #  Check whether the current best track record date number equals or exceeds
   #  that of the 1st forecast/analysis date number (for combined data sets).
   $istop = $i;
   if ($time_best_track >= $start_time_data && $ismissing == 0) {
      #  It does equal or exceed, so store the previous
      #  loop index and drop out of the loop.
      $istop = $i - 1;
      #warn "  HERE000\n";
      last;
   }

   #  Extract the year, month, day, hour for standard best track file.
   $for_yearb[$i] = substr($datestr[$i], 0, 4);
   $for_nmonthb[$i] = substr($datestr[$i], 4, 2);
   $for_dayb[$i] = substr($datestr[$i], 6, 2);
   $for_timeb[$i] = substr($datestr[$i], 8, 2);

   #  Determine which of the R34/50/64 wind radii are nonzero by summing each.
   $sum34 = $R34_Rad[$i][0] + $R34_Rad[$i][1] + $R34_Rad[$i][2] + $R34_Rad[$i][3];
   $sum50 = $R50_Rad[$i][0] + $R50_Rad[$i][1] + $R50_Rad[$i][2] + $R50_Rad[$i][3];
   $sum64 = $R64_Rad[$i][0] + $R64_Rad[$i][1] + $R64_Rad[$i][2] + $R64_Rad[$i][3];

   #  Define the indices of the R34/50/64 wind radii, depending on which values are nonzero.
   if ($sum64 == 0 && $sum50 == 0 && $sum34 >= 0) {
      #  R34 values may be present, R50 & R64 are not.
      $i34 = 3;
      $i50 = 2;
      $i64 = 1;
   } elsif ($sum64 == 0 && $sum50 > 0 && $sum34 > 0) {
      #  R34 & R50 values are present, R64 are not.
      $i34 = 2;
      $i50 = 3;
      $i64 = 1;
   } elsif ($sum64 > 0 && $sum50 > 0 && $sum34 > 0) {
      #  All 3 radii are present.
      $i34 = 1;
      $i50 = 2;
      $i64 = 3;
   }
   #  Loop over the radii to load the output 3D matrix.
   for my $j (0..$#{$R34_Rad[$i]}) {
      $r100_65_50_35_3Db[$i][0][$j] = 0;
      $r100_65_50_35_3Db[$i][$i34][$j] = $R34_Rad[$i][$j];
      $r100_65_50_35_3Db[$i][$i50][$j] = $R50_Rad[$i][$j];
      $r100_65_50_35_3Db[$i][$i64][$j] = $R64_Rad[$i][$j];
   }
}

#  If the storm name from the F/A file is empty, replace it w/that of the best track file.
if (length($tc_name) < 1 || lc($tc_name) !~ /[a-z]+/) { $tc_name = $st_name; }

#  Prepend the output arrays w/those of the best track file.  For those best track
#  arrays that weren't built in the previous FOR loop(s), use the range of valid
#  indices to omit later best track records.
@storm_id_strs    = (@storm_id[0..$istop],          @storm_id_strs);
@eye_diam         = (@eye_diameter[0..$istop],      @eye_diam);
@pres_outer_iso   = (@pres_outer_isobar[0..$istop], @pres_outer_iso);
@rad_outer_iso    = (@rad_outer_isobar[0..$istop],  @rad_outer_iso);
#  Check the best track file type to append the date, time parameters.
#  Standard best track file; these were loaded above, so use as-is.
@for_year         = (@for_yearb,                  @for_year);
@for_nmonth       = (@for_nmonthb,                @for_nmonth);
@for_day          = (@for_dayb,                   @for_day);
@for_time         = (@for_timeb,                  @for_time);

@for_lat          = (@latt[0..$istop],              @for_lat);
@for_lon          = (@long[0..$istop],              @for_lon);
@for_winds        = (@max_wind_kt[0..$istop],       @for_winds);
@for_rad          = (@RMW[0..$istop],               @for_rad);
@for_cp           = (@MSLP[0..$istop],              @for_cp);
@r100_65_50_35_3D = (@r100_65_50_35_3Db,            @r100_65_50_35_3D);

#  Sanity check on Mean RMW:
if ($MeanRMW <= 0) {
   #  It is invalid, so calculate a valid mean from the data.
   $sum34 = 0;
   $i34 = 0;
   for $i (0..$#for_rad) {
      if ($for_rad[$i] > 0) {
         $i34++;
         $sum34 += $for_rad[$i];
      }
   }
   if ($i34 > 0) { $MeanRMW = $sum34/$i34; }
   #warn "MeanRMW = $MeanRMW\n";
}

#  Check the output format.
#  Flag is -1, so print extended best track file.
#write_ext_best_track($trkfile,\@storm_id_strs,$tc_name,\@eye_diam,\@pres_outer_iso,\@rad_outer_iso,\@for_year,\@for_nmonth,\@for_day,\@for_time,\@for_lat,\@for_lon,\@for_winds,\@for_rad,\@r100_65_50_35_3D,\@for_cp);
write_cyc($trkfile,$tc_name,\@eye_diam,\@pres_outer_iso,\@rad_outer_iso,
   \@for_year,\@for_nmonth,\@for_day,\@for_time,\@for_lat,\@for_lon,\@for_winds,\@for_rad,
   \@r100_65_50_35_3D,\@for_cp, $tc_radius, $num_rad_bins, $num_dir_bins);

#  Calculate the number of records.
$numrecs = scalar @for_year;

#  Print the number of valid data records to standard error.
if ($numrecs == 0) {
   warn " NOTICE: No valid data records were found in $infile.\n";
} else {
   if ($appendflag == 1) {
      $done = "added";
   } else {
      $done = "written";
   }
   warn " NOTE: $numrecs valid data records were $done to $trkfile.\n";
}

print " Cyclone parameters written to: ${trkfile}";

exit;

#*******************************************************************************
# For debugging
sub print_matrix {
    my $array = shift;
    for my $i ( 0 .. $#{ $array } ) {
        my $row = $array->[$i];
        for my $j ( 0 .. $#{ $row } ) {
            print $array->[$i][$j], ' ';
        }
        #print "\n";
    }
}


#********************************************************************************
sub get_user_input {
    #
    #  This subroutine checks the command line syntax and obtains user input.
    #
    #  Calls:  [No external routines used]
    #  Called by:     main
    #
    #  Revision History:
    #  18 Apr 2006  Initial coding.  R.S. Linzell, PSI/Neptune Sciences Div.
    #  01 Nov 2011  Began adding support for DelftD WES.EXE output; modified the usage
    #               message; added initialization of new input parameters.  (RSL)
    #  02 Nov 2011  Modified the IF blocks to add support of Delft3D parameters.  (RSL)
    #  22 Feb 2013  Added support for an input best track file name.  (RSL)
    #  29 May 2022  WES.EXE is redundant; pared down to output only the cyclone
    #               file
    #
    #********************************************************************************

    #  Check the command line arguments.
    my $usage = " Usage: $0 InFile num_rad_bins num_dir_bins tc_radius [best_track_file]\n for Extended Best Track or Delft3D output.\n";

    #  Initialization of Delft3D parameters.
    my ($append, $ia, $nradbins, $ndirbins, $tcrad) = (0) x 5;
    my $btfile = '';

    #  If an incorrect number of arguments (<2 or >3) was entered, exit with a usage message.
    if ( ($#ARGV+1) != 4 && ($#ARGV+1) != 5) {
        die $usage;

    } else {

        #  Delft3D output:
        #  Input file is the 1st arg., and output wes.inp file is the 2nd..
        #  is ignored since it is not an append flag.
        $infile = $ARGV[0];
        #  The optional inputs were entered, so parse them out.
        $nradbins = $ARGV[1];  # Number of radial bins
        $ndirbins = $ARGV[2];  # Number of directional bins
        $tcrad = $ARGV[3];     # TC radius (km)
        if (($#ARGV+1) == 5) {
            $btfile = $ARGV[4];  # Input best track file name
        }
    }

    if ($DEBUG) { warn "In sub get_user_input: Infile: $infile  Outfile: $outfile Append Flag: $append.\n";
        warn " NumRadBins = $nradbins, NumDirBins = $ndirbins, TC_Rad = $tcrad, BT_File = $btfile.\n";
    }

    #  Warn user if an existing (nonzero sized) output file is going to be overwritten.
    if (-s $outfile && $append == 0) { warn " WARNING: Existing output file $outfile will be overwritten!\n"; }

    #  Return the input, output file names & append flag to calling program.
    return($infile, $nradbins, $ndirbins, $tcrad, $btfile);
}


#********************************************************************************
sub process_data {
    #
    #  This subroutine processes the input data.  The output values are stored in global
    #  variables used in various subroutines.
    #
    #  Input(s):      $infile     -- Name of input data file
    #
    #  Output(s):     $year       -- Year of advisory
    #                 $nmonth     -- Month number of advisory
    #                 $day        -- Day number (1 to 31) of advisory
    #                 $time       -- Hour of advisory
    #                 $minute     -- Minute of advisory
    #                 $latitude   -- Latitude of storm
    #                 $longitude  -- Longitude of storm
    #                 $cp         -- Central pressure of storm
    #                 $radius     -- Radius of maximum winds (cf. sub. get_radius)
    #                 @for_year   -- Array of years for each forecast time
    #                 @for_nmonth -- Array of months for each forecast time
    #                 @for_day    -- Array of days for each forecast time
    #                 @for_time   -- Array of hours for each forecast time
    #                 @for_minute -- Array of minutes for each forecast time
    #                 @for_lat    -- Array of latitudes for each forecast time
    #                 @for_lon    -- Array of longitudes for each forecast time
    #                 @for_cp     -- Array of central pressures for each forecast time
    #                 @for_rad    -- Array of radii of maximum winds for each forecast time
    #
    #  Called by:     main
    #
    #
    #  Revision History:
    #  18 Apr 2006  Initial coding.  R.S. Linzell, PSI/Neptune Sciences Div.
    #  01 Aug 2011  Added radius of maximum winds (RMW) algorithms based on latitude and
    #               maximum sustained wind speed.  (RSL)
    #  02 Nov 2011  Added code to support Delft3D output: added 'use'-ing of the 3rd-party
    #               Perl modules, Math::Trig (for storm direction of travel) and Date::Calc
    #               (for storm speed of travel); added code to calculate speed & direction
    #               where necessary; added new parameters needed for Delft3D code.  (RSL)
    #  09 May 2012  Added 'Mktime' to the 'use Date::Calc' pragma; added the new array
    #               '@r100_65_50_35_3D' to store the 3-D R34/R50/R64 wind radii; added
    #               the new array '@datenum' to store the date numbers for the analysis
    #               and forecast periods, and added calls to Mktime to load the @datenum
    #               array.  (RSL)
    #  10 May 2012  Added '$warning_number' & '$storm_id' to store those respective
    #               parameters.  (RSL)
    #  18 May 2012  Added '$warn_id' as the warning ID string based on the input file
    #               name without the extension.  (RSL)
    #  05 Sep 2012  Changed all debug 'print' statements to 'warn' statements; put a
    #               debug warn statement inside an IF block (it wasn't in an IF block
    #               previously); added "\n" to the Holland B warnings.  (RSL)
    #  01 Nov 2012  Commented-out the IF block where the eye diameter, if present,
    #               was used for the RMW.  (RSL)
    #  11 Jan 2013  Added an IF block about the code in which $tc_name is created, so
    #               that it isn't clobbered by invalid text after the header section.  (RSL)
    #  01 Feb 2013  Added an IF block to redefine $MISS for best track output.  (RSL)
    #  04 Feb 2013  Added def. of $eye_diameter to store any present EYE DIAMETER value
    #               for output to extended best track format.  (RSL)
    #  22 Feb 2013  Added a nested IF block at the end of the main WHILE loop to check
    #               for & correct invalid (i.e., zero) forecast pressure values.  (RSL)
    #
    #********************************************************************************

    my ($infile) = @_;

    #  Import the date calculation function.
    use Date::Calc qw(Delta_DHMS Mktime);

    if ($DEBUG) { warn " -- in sub. process_data InFile: $infile.\n"; }

    &RESET();				# Initialize global variables.

    #  Local definitions and initializations:
    $nm2km = 1.852; 		# Conversion for nmi to km
    $fi = 1;
    my ($dtime, $idump, $junkmonth, $numrec) = (0) x 4;
    my $ppp = "...";
    my @line = ();
    my $INIT = -9999;       # Initialization value for some parameters.
    $BINDEX = 1;            #  Atlantic basin index for sub. wind_pressure_relationship
    #  Wind speed values needed for the corresponding WES.EXE radii, in KT:
    @spd_needed = (35,50,65,100);
    $MISS = 1.0e30;         # Missing value code for WES.EXE files
    $MISS = 0;              # For optional extended best track output
    @tc_dir = ();           # TC direction of travel, for Delft3D output
    @tc_speed = ();         # TC speed of travel, for Delft3D output
    @r100_65_50_35 = ();    # Array of R100, R65, R50, R35 radii for Delft3D WES.EXE
    @r100_65_50_35_3D = (); # 3-D Array of R100, R65, R50, R35 radii for Delft3D WES.EXE
    @datenum = ();          # Array of date numbers (sec. since 1970-01-01 00:00:00)
    $warning_number = 0;    # Warning number (sequential through life of a storm)
    $storm_id = '';         # Storm ID (e.g., AL142011 = Atlantic Storm #14 for year 2011)
    $tc_name = '';          # Storm name (e.g., TROPICAL STORM FERNANDA)
    $warn_id = '';          # Warning ID string (input file name without extension)

    #  Constants for output file header:
    $rho = 1.15;            # Air density, kg/m^3
    $e = exp(1);            # Base of natural log
    $pn = 1013.25;          # Ambient air pressure (far from edge of storm), mb - Originally 1005
    $mb2npm2 = 100;         # Conversion factor, mb to N/m^2
    $kn2mps = 0.514444;     # Conversion factor, knots to m/s
    $BMin = 1;              # Minimum valid Holland "B" value
    $BMax = 2.5;            # Maximum valid Holland "B" value
    $B = $BMin;             # Holland "B" parameter (for output file header)
    $tstep = 1;             # Output time step, hr
    $CP_MIN_WTN = 882;      # Lowest known CP of Atlantic Basin (H. Wilma, 2005), mb
    $CP_MIN_WTP = 870;      # Lowest known CP of Pacific Basin (T. Tip, 1979), mb
    $CP_MIN = $CP_MIN_WTN;  # Minimum CP, mb
    $RMW_MAX = 500;         # Maximum radius of max. winds (typ. TS-force wind radius is 483 km), km
    $RMW_MAX /= $nm2km;     # Convert Max. RMW to nmi

    #  Parameters from input file:
    $cp = 0;
    @for_cp = ();
    @for_rad = ();
    @for_winds = ();
    $winds = 0;
    $radius = 0;
    $eye_diameter = -99;    # For optional extended best track output
    %nummonth = ('JAN',1,'FEB',2,'MAR',3,'APR',4,'MAY',5,'JUN',6,'JUL',7,'AUG',8,
        'SEP',9,'OCT',10,'NOV',11,'DEC',12);
    @daymonth = (0,31,28,31,30,31,30,31,31,30,31,30,31);
    @daymonthleap = (0,31,29,31,30,31,30,31,31,30,31,30,31);
    @dayyear = (0,31,59,90,120,151,181,212,243,273,304,334);
    #Original:  $lastline = "";

    #  Get current system time, and define associated time & date variables.
    ($nowsec,$nowmin,$nowhour,$nowmday,$nowmon,$nowyear,$nowwday,$nowyday,$nowisdst) = gmtime;
    if ($nowyear>=94) {$nowyear+=1900;} else {$nowyear+=2000;}
    if (isleap($nowyear)) {
        $now = $nowyear + ($nowyday+1 + ($nowhour/24.0) + ($nowmin/1440.))/366.0;
    } else {
        $now = $nowyear + ($nowyday+1 + ($nowhour/24.0) + ($nowmin/1440.))/365.0;
    }
    @badnames = ("RELOCATED","TEST","TD","TS","NONAME");


    #  Open file & loop over input records.
    open(INPUT, "<$infile") or die "ERROR opening input file $infile \($!\); stopped";

    #  Extract a warning ID from the file name, without the extension.
    ($warn_id, undef, undef) = fileparse($infile, ('.txt',".txt_[0-9]{14}"));

    while(<INPUT>) {
        chop;
        s/\cM//g;   # Get rid of ctrl-m's
        s/\cC/NNNN/g;
        s/\cA/NNNN/g;
        @line = split ' ', $_;	# Split line into fields.

        #  Check for end of message or input file.
        if (/NEXT\s+(ADVISORY|WARNINGS?)\s+AT/i && $idump == 0) {   # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/\$\$/i && $idump == 0) {   # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/\*\*[^\*]+\*\*/i && $idump == 0) {   # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/NNNN/i && $idump == 0) {   # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/^\.\.\.NOTICE\.\.\./i && $idump == 0) {   # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/^\s*REMARKS\:\s*$/i && $idump == 0) {  # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/^G\.\s+/i && $abcmessage && $idump == 0) {  # End of message for abc type messages
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/^\cC/ || /^\cA/ && $idump == 0) {  # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (/unsubscribe/i && $idump == 0) {  # End of message
            &DUMP();
            $idump++;
            #         &RESET();
        }
        if (eof && $idump == 0) {  # End of file
            &DUMP();
            $idump++;
            #         &RESET();
        }

        if (/^[\*\s]*((WT|TP|WH|AB)[A-Z][A-Z]\d{1,2})\s+([A-Z]{4})\s+(\d{6}|DDHHMM)(\s+[A-Z]{3}[\*\s]*|[\*\s]*)$/i) { # Start of Message
            #Original:         &DUMP();
            #Original:         &RESET();
            if (length($1) == 5) {  # Something like WTPN1 --> WTPN01
                $wmoheader = substr($1,0,4)."0".substr($1,4,1);
            } else {
                $wmoheader = $1;
            }
            if ($DEBUG) {warn "WMO: $wmoheader\n";}

            #  Define parameters for Pacific basin storms.
            if (uc($wmoheader) =~ /^WTP/) {
                $BINDEX = 2;             # Basin index for sub. wind_pressure_relationship
                $CP_MIN = $CP_MIN_WTP;	# Define min CP for Pacific
            }
        }


        #  Determine storm name.
        if (/(TROPICAL\s+DEPRESSION|TROPICAL\s+STORM|TROPICAL\s+CYCLONE|HURRICANE|TYPHOON)\s+(.*)(\s+|\s*\.\.\.\s*)(SPECIAL\s+|)(FORECAST\/ADVISORY|ADVISORY|DISCUSSION|WARNING|POSITION\s+ESTIMATE)(\s+NUMBER|\s+NR|)/i) {   # start of message

            $nlen = @namarr = split(/\s+/,$1);
            if ($nlen==1) {
                #Original:            $type = substr($namarr[0],0,3);
                $type = $namarr[0];
            } elsif ($nlen==2) {
                #Original:            $type = substr($namarr[1],0,3);
                $type = "$namarr[0] $namarr[1]";
            }
            $nlen = @namarr = split(/\s+/,$2);
            $name = $namarr[0];
            if ($nlen==2) {   # Check for a secondary name in parentheses (PGTW)
                if ($namarr[1]=~/\((.*)\)/) {
                    $secondname = $1;
                    if (&ISDEPRESSION($name) && !&ISDEPRESSION($secondname)) {
                        $temp = $name;         # Use the secondary name if its better
                        $name = $secondname;
                        $secondname = $temp;
                    }
                    #Original:               if ($DEBUG) {print "$name $secondname\n"}
                    if ($DEBUG) {warn " Data for $name $secondname:\n";}
                }
            }
            $name =~ s/\(//;
            $name =~ s/\)//;
            #  Store this forecast/advisory number.
            $warning_number = $line[$#line];
            #Original:         if ($DEBUG) {print "$type $name\n";}
            if ($DEBUG) {warn " Data from Forecast/Advisory Number $warning_number for $type $name:\n";}

            #  Build the storm name if it hasn't been done.
            if (length($tc_name) < 1) {
                #  First, check for a multi-word storm type.
                if ($type =~ /\s+/) {
                    #  This is a multi-part storm type, so split the type into separate words.
                    my @part = split ' ', $type;
                    #  Loop over words....
                    foreach $word (@part) {
                        #  Append the initalcap word, then a space.
                        $tc_name .= ucfirst(lc($word));
                        $tc_name .= ' ';
                    }
                    #  Append the initialcap storm name.
                    $tc_name .= ucfirst(lc($name));
                } else {
                    #  Single word storm type, so simply concatenate the initialcap storm type and name.
                    $tc_name = ucfirst(lc($type)) . ' ' . ucfirst(lc($name));
                }
            }
            if ($DEBUG) { warn "TC_NAME = $tc_name.\n"; }
        }

        #  Determine the storm ID (e.g., AL142011, EP062011)
        if (/^NWS/ && /HURRICANE CENTER/) {
            #  Storm ID record
            $storm_id = $line[$#line];
            if ($DEBUG) { warn "  Storm ID: $storm_id\n"; }
        }

        #  Obtain the present movement of the storm system for Delft3D support.
        if ($_ =~ /PRESENT MOVEMENT/) {
            $tc_dir[0] = $line[$#line-4]; # was [$fi]
            $tc_speed[0] = $line[$#line-1]; # was [$fi]
            if ($DEBUG) { warn "  Present Movement, FI = $fi, TC_Dir = $tc_dir[0]; TC_Speed = $tc_speed[0]\n"; }
        }

        #  Determine central pressure, and perform gross error check.
        #  if ($_ =~ /(ESTIMATED\s+MINIMUM\s+CENTRAL\s+PRESSURE)/) {
        if ($_ =~ /(CENTRAL\s+PRESSURE)/) {
            $cp = $line[4];
            if ($cp < $CP_MIN) {
                warn "  WARNING: Central Pressure \($cp\) < Minimum \($CP_MIN mb\); resetting to Minimum.\n";
                $cp = $CP_MIN;
            } elsif ($cp > $pn) {
                warn "  WARNING: Central Pressure \($cp\) > Ambient \($pn mb\); resetting to Ambient.\n";
                $cp = $pn;
            }
        }

        #  Find eye diameter, if present, and store for output to extended best track file
        #  if that option was selected.
        if (/(EYE\s+DIAMETER)/i) {
            $eye_diameter = $line[2];
            #  $radius = $line[2] / 2;
            #  if ($line[$#line] =~ /NM/) { $radius *= $nm2km; }	 # Convert nmi to km.
        }

        #  Find the max. wind speed, and radius of maximum winds if not already found.
        if (/(MAX\s+SUSTAINED\s+WINDS)/) {
            $winds = $line[3];
            if ($radius == 0) { $radius = get_radius($winds) * $nm2km; }
        }  # end if($_

        # forecasts are considered active when the word FORECAST is found.
        # This is somewhat restrictive, but is necessary as garbled
        # observations would otherwise look into the forecast data as if
        # it was an observation.  This is avoided since no forecasts are
        # output in DUMP if there was no valid observation.
        if (/AT\s+(\d{6,6})\s+UTC/i && $day eq "???" && $time eq "???" && !$forecast_active) {  # time
            $day = substr($1,0,2);
            $time = substr($1,2,2);
            $minute = substr($1,4,2);
            ($nowsec,$nowmin,$nowhour,$nowmday,$nowmon,$nowyear,$nowwday,$nowyday,$nowisdst) = gmtime(time);
            $nmonth = $nowmon+1;
            if ($nowyear > 93) {$year = $nowyear+1900;}
            else {$year = $nowyear+2000;}
            ($year,$nmonth) = &GETNMONTH($year,$nmonth,$nowmday,$day);
            $year2 = substr($year,2,2);
        }
        if (/AT\s+(\d{8,8})Z/i && $day eq "???" && $time eq "???" && !$forecast_active) {  # time
            $nmonth = substr($1,0,2);
            $day = substr($1,2,2);
            $time = substr($1,4,2);
            $minute = substr($1,6,2);
            ($nowsec,$nowmin,$nowhour,$nowmday,$nowmon,$nowyear,$nowwday,$nowyday,$nowisdst) = gmtime(time);
            if ($nowyear > 93) {$year = $nowyear+1900;}
            else {$year = $nowyear+2000;}
            ($year,$junkmonth) = &GETNMONTH($year,$nmonth,$nowmday,$day);
            $year2 = substr($year,2,2);
        }

        # Get time of format "B. 082030Z" or "B. 14/2226Z"
        if (/(\d{1,4}|MIDNIGHT|NOON)\s*(Z| AM| PM| UTC)\s*(HST |PST |PDT |MST |MDT |CST |CDT |EST |EDT |AST |ADT |)(SUN|MON|TUE|WED|THU|FRI|SAT)\s+(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)\s+([^ ]*)\s+(\d\d\d\d)/i) {    # Time
            if ("\U$1" eq "MIDNIGHT") {
                $time = 0;
            } elsif ("\U$1" eq "NOON") {
                $time = 12;
            } else {
                $timelen = length($1);
                if ($timelen < 3) {$time = $1; $minute="00";}
                elsif ($timelen == 3) {
                    $time=substr($1,0,1);
                    $minute=substr($1,1,2);
                } else {$time=substr($1,0,2); $minute=substr($1,2,2);}

            }
            $dow = $4;
            $month = "\U$5";
            $nmonth = $nummonth{$month};
            $day = $6;
            $year = $7;
            $year2 = substr($year,2,2);
            if ("\U$2" eq " PM" && $1 < 12) {$time += 12;}
            if ("\U$2" eq " AM" && $1 == 12) {$time = "00";}
            if ("\U$3" eq "HST ") {$time += 10;}   # convert to UT
            if ("\U$3" eq "PST ") {$time += 8;}
            if ("\U$3" eq "MST ") {$time += 7;}
            if ("\U$3" eq "CST ") {$time += 6;}
            if ("\U$3" eq "EST ") {$time += 5;}
            if ("\U$3" eq "AST ") {$time += 4;}
            if ("\U$3" eq "PDT ") {$time += 7;}
            if ("\U$3" eq "MDT ") {$time += 6;}
            if ("\U$3" eq "CDT ") {$time += 5;}
            if ("\U$3" eq "EDT ") {$time += 4;}
            if ("\U$3" eq "ADT ") {$time += 3;}
            while ($time >= 24) {$time -= 24; $day += 1;}
            if ($year2 % 4 == 0) {
                if ($day > $daymonthleap[$nmonth]) {
                    $day = $day % $daymonthleap[$nmonth];
                    $nmonth++;
                }
            } else {
                if ($day > $daymonth[$nmonth]) {
                    $day = $day % $daymonth[$nmonth];
                    $nmonth++;
                }
            }
            if ($nmonth > 12) {$nmonth=1; ++$year; ++$year2;}
            if ($DEBUG) {print "Time: $time $dow $month $day $year\n"; }
            push @datenum, (Mktime($year, $nmonth, $day, $time, $minute, 0.0));
        }
        # Get observation of type "INITIAL     25/0300Z 21.6N 166.6W    25 KTS"
        if (m#INITIAL +(\D*)(\d*)/(\d*)Z(\D*)(\d*\.\d*)( +N| +S|N|S)(\D*)(\d*\.\d*)( E| W|E|W)(\D*)(\d*)(\D*)(KTS|KT|MPH)#i) {
            $latstr = $5.$6;
            $lonstr = $8.$9;
            $latstr =~ s/ORTH//gi;  # NORTH --> N
            $latstr =~ s/OUTH//gi;  # SOUTH --> S
            $lonstr =~ s/AST//gi;   # EAST  --> E
            $lonstr =~ s/EST//gi;   # WEST  --> W
            $latitude = substr($latstr,0,length($latstr)-1);
            $longitude = substr($lonstr,0,length($lonstr)-1);
            $latns = substr($latstr,length($latstr)-1,1);
            $lonew = substr($lonstr,length($lonstr)-1,1);
            $winds = $11;
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ("\U$13" eq "MPH") {$winds = int($winds*0.869+0.5);}
            if ($DEBUG) {warn "Initial $latitude $latns $longitude $lonew $winds\n";}
        }

        $rmw0 = 66.785 - (0.09102 * $winds) + (1.0619 * (abs($latitude) - 25));  # Knaff & Zehr (2007) in km
        $rmw1 = (35.37 - (0.111 * $winds) + (0.570 * (abs($latitude) - 25)))  * $nm2km;  # Gross et al. (2004) in nmi converted to km
        $rmw2 = 51.6 * exp((-0.0223 * $winds * $kn2mps) + (0.0281 * abs($latitude)));     # Willoughby & Rahn (2004) in km
        #warn "Analysis: Radius = $radius; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n";
        if ($DEBUG) {warn "Analysis: Radius = $radius; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n"; }

        # Get forecasts of type "12HR VT     25/1200Z 22.1N 167.6W    25 KTS"
        if (m#(\d*)HR VT(\D*)(\d*)/(\d*)Z(\D*)(\d*\.\d*)( N| S|N|S)(\D*)(\d*\.\d*)( E| W|E|W)(\D*)(\d*)(\D*)(KTS|KT|MPH)#i) {
            $forecast_active = 1;  # this is a forecast
            $latstr = $6.$7;
            $lonstr = $9.$10;
            $latstr =~ s/ORTH//gi;  # NORTH --> N
            $latstr =~ s/OUTH//gi;  # SOUTH --> S
            $lonstr =~ s/AST//gi;   # EAST  --> E
            $lonstr =~ s/EST//gi;   # WEST  --> W
            $nforecast = $#for_day+1;
            $for_name[$nforecast] = "\U$name";
            $for_day[$nforecast] = $3;
            $timelen = length($4);
            if ($timelen == 4) {
                $for_time[$nforecast] = substr($4,0,2);
                $for_minute[$nforecast] = substr($4,2,2)
            } elsif ($timelen == 3) {
                $for_time[$nforecast] = substr($4,0,1);
                $for_minute[$nforecast] = substr($4,1,2)
            } elsif ($timelen == 2) {
                $for_time[$nforecast] = substr($4,0,2);
                $for_minute[$nforecast] = "00";
            } elsif ($timelen == 1) {
                $for_time[$nforecast] = $4;
                $for_minute[$nforecast] = "00";
            }
            ($for_year[$nforecast],$for_nmonth[$nforecast]) = &GETNMONTH($year,$nmonth,$day,$for_day[$nforecast]);
            $for_year2[$nforecast] = substr($for_year[$nforecast],2,2);
            $for_lat[$nforecast] = substr($latstr,0,length($latstr)-1);
            $for_lon[$nforecast] = substr($lonstr,0,length($lonstr)-1);
            $for_latns[$nforecast] = substr($latstr,length($latstr)-1,1);
            $for_lonew[$nforecast] = substr($lonstr,length($lonstr)-1,1);
            $for_winds[$nforecast] = $12;
            if (uc($for_lonew[$nforecast]) eq "W") { $for_lon[$nforecast] = 360 - $for_lon[$nforecast]; }
            if ("\U$14" eq "MPH") {$for_winds[$nforecast] = int($for_winds[$nforecast]*0.869+0.5);}
            $for_gusts[$nforecast] = "???";

            #  Update the storm direction & speed.
            if ($nforecast >= 1) {
                if ($nforecast == 1) {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($longitude), deg2rad(90 - $latitude));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($year,$nmonth,$day,$time,0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$longitude,$latitude,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$year,$nmonth,$day,$time);
                } else {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($for_lon[$nforecast-1]), deg2rad(90 - $for_lat[$nforecast-1]));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1],0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$for_lon[$nforecast-1],$for_lat[$nforecast-1],$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1]);
                }
                #  Ending lon, lat (from the current record), converted to radians:
                #@end = (deg2rad($for_lon[$nforecast]), deg2rad(90 - $for_lat[$nforecast]));
                #  Great circle distance, in km:
                #$km = great_circle_distance(@start, @end, 6378.137);
                #  Bearing from start to end point (radians):
                #$dir = great_circle_bearing(@start, @end);
                #  Convert direction to degrees:
                #$tc_dir[$nforecast] = rad2deg($dir);
                #if ($Dd == 0) {
                #  $dtime = $Dh;
                #} else {
                #  $dtime = $Dh + ($Dd * 24);
                #}
                #  TC translation speed is distance / time:
                #$tc_speed[$nforecast] = ($km/$nm2km)/$dtime;

                if ($DEBUG) {warn " ^^^^^^vvv FI = $fi; TC_Dir1 = $tc_dir[$nforecast]; TC_Speed = $tc_speed[$nforecast]; Delta_Time = $dtime hr, NForecast = $nforecast.\n";}
                @start = (); @end = (); $km = $INIT; $dir = $INIT; $dtime = $INIT;
            }

            #Original:         if ($DEBUG) {print "FC: $latstr $lonstr\n";}
            if ($DEBUG) {warn "FC $nforecast: $for_lat[$nforecast] $for_latns[$nforecast] $for_lon[$nforecast] $for_lonew[$nforecast]\n";}

            #  Find radius of max. winds.
            $_ = <INPUT>;
            @line = split ' ', $_;

            #  Find the max. wind speed.
            if ($_ =~ /(MAX\s+WIND)/) {
                $for_winds[$fi] = $line[2];
                $for_rad[$fi] = get_radius($for_winds[$fi]) * $nm2km;
                #  Use both RMW algorithms to check whether either is better than the above.
                $rmw0 = 66.785 - (0.09102 * $for_winds[$fi]) + (1.0619 * (abs($for_lat[$nforecast]) - 25));  # Knaff & Zehr (2007) in km
                $rmw1 = (35.37 - (0.111 * $for_winds[$fi]) + (0.570 * (abs($for_lat[$nforecast]) - 25)))  * $nm2km;  # Gross et al. (2004) in nmi converted to km
                $rmw2 = 51.6 * exp((-0.0223 * $for_winds[$fi] * $kn2mps) + (0.0281 * abs($for_lat[$nforecast])));     # Willoughby & Rahn (2004) in km
                #warn "Forecast $fi: Radius = $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n";
                if ($DEBUG) { warn " Forecast $fi Wind, Radius = $for_winds[$fi], $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n"; }
            }  # end if($_
            $fi++;
        }
        # Get forecasts of type "FORECAST VALID  25/1200Z 22.1N 167.6W "
        #                       OUTLOOK VALID 24/1800Z 27.0N  69.5W
        if (m#[FORECAST|OUTLOOK] +VALID(\D*)(\d+)/(\d+)Z(\D*)(\d+\.\d*)( N| S|N|S)(\D*)(\d+\.\d*)( E| W|E|W)#i) {
            $forecast_active = 1;  # this is a forecast
            $latstr = $5.$6;
            $lonstr = $8.$9;
            $latstr =~ s/ORTH//gi;  # NORTH --> N
            $latstr =~ s/OUTH//gi;  # SOUTH --> S
            $lonstr =~ s/AST//gi;   # EAST  --> E
            $lonstr =~ s/EST//gi;   # WEST  --> W
            $nforecast = $#for_day+1;
            $for_name[$nforecast] = "\U$name";
            $for_day[$nforecast] = $2;
            $timelen = length($3);
            if ($timelen == 4) {
                $for_time[$nforecast] = substr($3,0,2);
                $for_minute[$nforecast] = substr($2,2,2)
            } elsif ($timelen == 3) {
                $for_time[$nforecast] = substr($3,0,1);
                $for_minute[$nforecast] = substr($3,1,2)
            } elsif ($timelen == 2) {
                $for_time[$nforecast] = substr($3,0,2);
                $for_minute[$nforecast] = "00";
            } elsif ($timelen == 1) {
                $for_time[$nforecast] = $3;
                $for_minute[$nforecast] = "00";
            }
            ($for_year[$nforecast],$for_nmonth[$nforecast]) = &GETNMONTH($year,$nmonth,$day,$for_day[$nforecast]);
            $for_year2[$nforecast] = substr($for_year[$nforecast],2,2);
            $for_lat[$nforecast] = substr($latstr,0,length($latstr)-1);
            $for_lon[$nforecast] = substr($lonstr,0,length($lonstr)-1);
            $for_latns[$nforecast] = substr($latstr,length($latstr)-1,1);
            $for_lonew[$nforecast] = substr($lonstr,length($lonstr)-1,1);
            if (uc($for_lonew[$nforecast]) eq "W") { $for_lon[$nforecast] = 360 - $for_lon[$nforecast]; }

            #  Update the storm direction & speed.
            if ($nforecast >= 1) {
                if ($nforecast == 1) {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($longitude), deg2rad(90 - $latitude));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($year,$nmonth,$day,$time,0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$longitude,$latitude,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$year,$nmonth,$day,$time);
                } else {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($for_lon[$nforecast-1]), deg2rad(90 - $for_lat[$nforecast-1]));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1],0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$for_lon[$nforecast-1],$for_lat[$nforecast-1],$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1]);
                }
                #  Ending lon, lat (from the current record), converted to radians:
                #@end = (deg2rad($for_lon[$nforecast]), deg2rad(90 - $for_lat[$nforecast]));
                #  Great circle distance, in km:
                #$km = great_circle_distance(@start, @end, 6378.137);
                #  Bearing from start to end point (radians):
                #$dir = great_circle_bearing(@start, @end);
                #  Convert direction to degrees:
                #$tc_dir[$nforecast] = rad2deg($dir);
                #if ($Dd == 0) {
                #  $dtime = $Dh;
                #} else {
                #  $dtime = $Dh + ($Dd * 24);
                #}
                #  TC translation speed is distance / time:
                #$tc_speed[$nforecast] = ($km/$nm2km)/$dtime;

                if ($DEBUG1) {warn " ^^^^^^vvv-- FI = $fi; TC_Dir2 = $tc_dir[$nforecast]; TC_Speed = $tc_speed[$nforecast]; Delta_Time = $dtime hr, NForecast = $nforecast.\n"; }
                @start = (); @end = (); $km = $INIT; $dir = $INIT; $dtime = $INIT;
            }

            #Original:         if ($DEBUG) {print "F/OA $latstr $lonstr\n";}
            if ($DEBUG) {warn "F/OA $nforecast: $for_lat[$nforecast] $for_latns[$nforecast] $for_lon[$nforecast] $for_lonew[$nforecast]\n";}

            #  Find radius of max. winds.
            $_ = <INPUT>;
            @line = split ' ', $_;

            #  Find the max. wind speed.
            if ($_ =~ /(MAX\s+WIND)/) {
                $for_winds[$fi] = $line[2];
                $for_rad[$fi] = get_radius($for_winds[$fi]) * $nm2km;
                #  Use both RMW algorithms to check whether either is better than the above.
                $rmw0 = 66.785 - (0.09102 * $for_winds[$fi]) + (1.0619 * (abs($for_lat[$nforecast]) - 25));  # Knaff & Zehr (2007) in km
                $rmw1 = (35.37 - (0.111 * $for_winds[$fi]) + (0.570 * (abs($for_lat[$nforecast]) - 25)))  * $nm2km;  # Gross et al. (2004) in nmi converted to km
                $rmw2 = 51.6 * exp((-0.0223 * $for_winds[$fi] * $kn2mps) + (0.0281 * abs($for_lat[$nforecast])));     # Willoughby & Rahn (2004) in km
                #warn "Forecast $fi: Radius = $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n";
                if ($DEBUG1) { warn " Forecast $fi Wind, Radius = $for_winds[$fi], $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n"; }
            }  # end if($_
            $fi++;
        }
        #Original:      if (/CENTER +(|RE)LOCATED +NEAR +(\D*)(\d*\.\d*)( +N| +S|N|S)(\D*)(\d*\.\d*)( +E| +W|E|W) +AT/i && !$forecast_active) {	# Position
        #2nd:      if (/(CENTER\s+(|RE)LOCATED\s+NEAR)(\s+\d*\.\d*)(\s*N|S)(\s*\d*\.\d*)(\s*E|W)(\s*AT)/i && !$forecast_active) {	# Position
        if (/(CENTER\s+LOCATED\s+NEAR)/i && !$forecast_active) {
            @line = split ' ', $_;
            #Original:         $latstr = $3.$4;
            #Original:         $lonstr = $6.$7;
            $latstr = $line[($#line - 3)];
            $lonstr = $line[($#line - 2)];
            $latstr =~ s/ORTH//gi;  # NORTH --> N
            $latstr =~ s/OUTH//gi;  # SOUTH --> S
            $lonstr =~ s/AST//gi;   # EAST  --> E
            $lonstr =~ s/EST//gi;   # WEST  --> W
            $latitude = substr($latstr,0,length($latstr)-1);
            $longitude = substr($lonstr,0,length($lonstr)-1);
            $latns = substr($latstr,length($latstr)-1,1);
            $lonew = substr($lonstr,length($lonstr)-1,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "Position: $latstr $latitude $latns $lonstr $longitude $lonew; FI = $fi; NForecast = $nforecast\n";}
        }
        #Original:      if (/\sPOSITION(\D*)(\d*\.\d*)( +N| +S|N|S)(\D*)(\d*\.\d*)( +E| +W|E|W)/i && !$forecast_active) {    # Position
        if (/\sPOSITION(\D*)(\d*\.\d*)(\s*N|S)(\s*\d*\.\d*)(\s*E|W)/i && !$forecast_active) {    # Position
            $latstr = $2.$3;
            $lonstr = $5.$6;
            $latstr =~ s/ORTH//gi;  # NORTH --> N
            $latstr =~ s/OUTH//gi;  # SOUTH --> S
            $lonstr =~ s/AST//gi;   # EAST  --> E
            $lonstr =~ s/EST//gi;   # WEST  --> W
            $latitude = substr($latstr,0,length($latstr)-1);
            $longitude = substr($lonstr,0,length($lonstr)-1);
            $latns = substr($latstr,length($latstr)-1,1);
            $lonew = substr($lonstr,length($lonstr)-1,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "Position2: $latitude $latns $longitude $lonew\n";}
        }
        if ($doubleline =~ /NEAR\s+(\d*\.\d*)\s+(NORTH|SOUTH)\s+(\d*\.\d*)\s+(EAST|WEST)/i && !$forecast_active) {
            $latitude = $1;
            $latns = substr($2,0,1);
            if ($DEBUG) {warn "*** $latitude $latns\n";}
            $longitude = $3;
            $lonew = substr($4,0,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "*** $longitude $lonew\n";}
        }
        if ($doubleline =~ /NEAR\s+(\d*\.\d*)\s+(NORTH|SOUTH)\s+(\d*\.\d*)\s+(EAST|WEST)\s+AT\s+(\d{4,6})\s+UTC/i && !$forecast_active) {
            $latitude = $1;
            $latns = substr($2,0,1);
            if ($DEBUG) {warn "*** $latitude $latns\n";}
            $longitude = $3;
            $lonew = substr($4,0,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "*** $longitude $lonew\n";}
            if ($nmonth ne "???" && $day ne "???" && $year ne "???") {
                if (length($5) == 4) {
                    ($year,$nmonth) = &GETNMONTH($year,$nmonth,$day,$day);
                    $year2[$nforecast] = substr($year,2,2);
                    $time = substr($5,0,2);
                    $minute = substr($5,2,2);
                } elsif (length($5) == 6) {
                    ($year,$nmonth) = &GETNMONTH($year,$nmonth,$day,substr($5,0,2));
                    $year2[$nforecast] = substr($year,2,2);
                    $day = substr($5,0,2);
                    $time = substr($5,2,2);
                    $minute = substr($5,4,2);
                }
            }
        }
        if ($doubleline =~ /AT\s+(\d{4,6})\s+UTC\D*NEAR\s+(\d*\.\d*)\s+(NORTH|SOUTH)\s+(\d*\.\d*)\s+(EAST|WEST)/i && !$forecast_active) {
            $latitude = $2;
            $latns = substr($3,0,1);
            if ($DEBUG) {warn "*** $latitude $latns\n";}
            $longitude = $4;
            $lonew = substr($5,0,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "*** $longitude $lonew\n";}
            if ($nmonth ne "???" && $day ne "???" && $year ne "???") {
                if (length($1) == 4) {
                    ($year,$nmonth) = &GETNMONTH($year,$nmonth,$day,$day);
                    $year2[$nforecast] = substr($year,2,2);
                    $time = substr($1,0,2);
                    $minute = substr($1,2,2);
                } elsif (length($1) == 6) {
                    ($year,$nmonth) = &GETNMONTH($year,$nmonth,$day,substr($1,0,2));
                    $year2[$nforecast] = substr($year,2,2);
                    $day = substr($1,0,2);
                    $time = substr($1,2,2);
                    $minute = substr($1,4,2);
                }
            }
        }
        if (/FORECAST\s+POSITION\s+AT\s+(\d{6,6})\s+UTC\s+(|NEAR\s+)(\d*\.\d*)\s+(NORTH|SOUTH)\s+(\d*\.\d*)\s+(EAST|WEST)/i) {
            $forecast_active=1;  # this is a forecast
            $nforecast = $#for_day+1;
            $for_name[$nforecast] = "\U$name";
            $for_day[$nforecast] = substr($1,0,2);
            $for_time[$nforecast] = substr($1,2,2);
            $for_minute[$nforecast] = substr($1,4,2);
            ($for_year[$nforecast],$for_nmonth[$nforecast]) = &GETNMONTH($year,$nmonth,$day,$for_day[$nforecast]);
            $for_year2[$nforecast] = substr($for_year[$nforecast],2,2);
            $for_lat[$nforecast] = $3;
            $for_lon[$nforecast] = $5;
            $for_latns[$nforecast] = substr($4,0,1);
            $for_lonew[$nforecast] = substr($6,0,1);
            if (uc($for_lonew[$nforecast]) eq "W") { $for_lon[$nforecast] = 360 - $for_lon[$nforecast]; }

            #  Update the storm direction & speed.
            if ($nforecast >= 1) {
                if ($nforecast == 1) {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($longitude), deg2rad(90 - $latitude));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($year,$nmonth,$day,$time,0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$longitude,$latitude,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$year,$nmonth,$day,$time);
                } else {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($for_lon[$nforecast-1]), deg2rad(90 - $for_lat[$nforecast-1]));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1],0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$for_lon[$nforecast-1],$for_lat[$nforecast-1],$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1]);
                }
                #  Ending lon, lat (from the current record), converted to radians:
                #@end = (deg2rad($for_lon[$nforecast]), deg2rad(90 - $for_lat[$nforecast]));
                #  Great circle distance, in km:
                #$km = great_circle_distance(@start, @end, 6378.137);
                #  Bearing from start to end point (radians):
                #$dir = great_circle_bearing(@start, @end);
                #  Convert direction to degrees:
                #$tc_dir[$nforecast] = rad2deg($dir);
                #if ($Dd == 0) {
                #  $dtime = $Dh;
                #} else {
                #  $dtime = $Dh + ($Dd * 24);
                #}
                #  TC translation speed is distance / time:
                #$tc_speed[$nforecast] = ($km/$nm2km)/$dtime;

                if ($DEBUG) {warn " ^^^^^^vvv FI = $fi; TC_Dir3 = $tc_dir[$nforecast]; TC_Speed = $tc_speed[$nforecast]; Delta_Time = $dtime hr, NForecast = $nforecast.\n";}
                @start = (); @end = (); $km = $INIT; $dir = $INIT; $dtime = $INIT;
            }

            if ($DEBUG) {warn "-F/OA $nforecast: $for_lat[$nforecast] $for_latns[$nforecast] $for_lon[$nforecast] $for_lonew[$nforecast]\n";}

            #  Find radius of max. winds.
            $_ = <INPUT>;
            @line = split ' ', $_;

            #  Find the max. wind speed.
            if ($_ =~ /(MAX\s+WIND)/) {
                $for_winds[$fi] = $line[2];
                $for_rad[$fi] = get_radius($for_winds[$fi]) * $nm2km;
                #  Use both RMW algorithms to check whether either is better than the above.
                $rmw0 = 66.785 - (0.09102 * $for_winds[$fi]) + (1.0619 * (abs($for_lat[$nforecast]) - 25));  # Knaff & Zehr (2007) in km
                $rmw1 = (35.37 - (0.111 * $for_winds[$fi]) + (0.570 * (abs($for_lat[$nforecast]) - 25)))  * $nm2km;  # Gross et al. (2004) in nmi converted to km
                $rmw2 = 51.6 * exp((-0.0223 * $for_winds[$fi] * $kn2mps) + (0.0281 * abs($for_lat[$nforecast])));     # Willoughby & Rahn (2004) in km
                #warn "Forecast $fi: Radius = $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n";
                if ($DEBUG) { warn " Forecast $fi Wind, Radius = $for_winds[$fi], $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n"; }
            }  # end if($_
            $fi++;
        }
        if (/LATITUDE +(\D*)(\d*\.\d*) +(NORTH|SOUTH)/i && !$forecast_active) {    # Position
            $latitude = $2;
            $latns = substr($3,0,1);
            if ($DEBUG) {warn "*** $latitude $latns\n";}
            #last switch;
        }
        if (/LONGITUDE +(\D*)(\d*\.\d*) +(WEST|EAST)/i && !$forecast_active) {    # Position
            $longitude = $2;
            $lonew = substr($3,0,1);
            if (uc($lonew) eq "W") { $longitude = 360 - $longitude; }
            if ($DEBUG) {warn "*** $longitude $lonew\n";}
            #last switch;
        }
        if ($doubleline =~ /WIND(\D*)(\d*)(\D*)(KT|KTS|MPH|KNOT)/i && !($doubleline =~ /REPORT/i) && ($winds eq "???" || $forecast_active)) {
            if ($forecast_active) {  # Assumes that forecast position came first
                $for_winds[$nforecast] = $2;
                if (!&ISNUMBER($for_gusts[$nforecast])) {
                    $for_gusts[$nforecast]="???";
                }
                if ("\U$4" eq "MPH") {$for_winds[$nforecast] = int($for_winds[$nforecast]*0.869+0.5);}
                if ($DEBUG) {warn "Forecast Wind: ".$for_winds[$nforecast]."\n";}
            } else {
                $winds = $2;
                if ("\U$4" eq "MPH") {$winds = int($winds*0.869+0.5);}
                if ($DEBUG) {warn "--Winds: $winds  Gusts: $gusts\n";}
            }
        }
        if (/(MAX\.|MAX|MAXIMUM)\s+(SUSTAINED\s+|)(WIND|WINDS)\s*\.\.\.\s*(\d+)\s*(KT|KTS|MPH|KNOT|KM\/HR)/i && (!$forecast_active)) {
            $winds = $4;
            if ("\U$5" eq "MPH") {$winds = int($winds*0.869+0.5);}
            if ("\U$5" eq "KM/HR") {$winds = int($winds*0.55+0.5);}
            if ($DEBUG) {warn "Winds: $winds  Gusts: $gusts\n";}
        }
        if ($doubleline =~ /GUST(\D*)(\d*)(\D*)(KT|KTS|MPH|KNOT)/i && !($doubleline =~ /REPORT/i) && ($gusts eq "???" || $forecast_active)) {
            if ($forecast_active) {  # Assumes that forecast position came first
                $for_gusts[$nforecast] = $2;
                if (!&ISNUMBER($for_winds[$nforecast])) {
                    $for_winds[$nforecast]="???";
                }
                if ("\U$4" eq "MPH") {
                    $for_gusts[$nforecast] = int($for_gusts[$nforecast]*0.869+0.5);
                }
                if ($DEBUG) {warn "Forecast Gusts: ".$for_gusts[$nforecast]."\n";}
            } else {
                $gusts = $2;
                if ("\U$4" eq "MPH") {$gusts = int($gusts*0.869+0.5);}
                if ($DEBUG) {warn "---Winds: $winds  Gusts: $gusts\n";}
            }
        }
        if (m#[FORECAST|OUTLOOK] +VALID +(\D*)(\d*)/(\d*)Z +(\D*)(\d*\.\d*)( +N| +S|N|S)(\D*)(\d*\.\d*)( +E| +W|E|W)#i) {  # Forecast
            $forecast_active=1;  # this is a forecast
            $nforecast = $#for_day+1;
            $for_name[$nforecast] = "\U$name";
            $for_day[$nforecast] = $2;
            $timelen = length($3);
            if ($timelen == 4) {
                $for_time[$nforecast] = substr($3,0,2);
                $for_minute[$nforecast] = substr($3,2,2)
            } elsif ($timelen == 3) {
                $for_time[$nforecast] = substr($3,0,1);
                $for_minute[$nforecast] = substr($3,1,2)
            } elsif ($timelen == 2) {
                $for_time[$nforecast] = substr($3,0,2);
                $for_minute[$nforecast] = "00";
            } elsif ($timelen == 1) {
                $for_time[$nforecast] = $3;
                $for_minute[$nforecast] = "00";
            }
            ($for_year[$nforecast],$for_nmonth[$nforecast]) = &GETNMONTH($year,$nmonth,$day,$for_day[$nforecast]);
            $for_year2[$nforecast] = substr($for_year[$nforecast],2,2);
            $for_lat[$nforecast] = $5;
            $for_lon[$nforecast] = $8;
            $for_latns[$nforecast] = $6;
            $for_lonew[$nforecast] = $9;
            if (uc($for_lonew[$nforecast]) eq "W") { $for_lon[$nforecast] = 360 - $for_lon[$nforecast]; }

            #  Update the storm direction & speed.
            if ($nforecast >= 1) {
                if ($nforecast == 1) {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($longitude), deg2rad(90 - $latitude));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($year,$nmonth,$day,$time,0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$longitude,$latitude,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$year,$nmonth,$day,$time);
                } else {
                    #  Starting lon, lat (from previous record) converted to radians:
                    #@start = (deg2rad($for_lon[$nforecast-1]), deg2rad(90 - $for_lat[$nforecast-1]));
                    #  Calculate time difference between the current and previous records.
                    #($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1],0,0,$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],0,0);
                    ($tc_speed[$nforecast],$tc_dir[$nforecast]) = storm_velocity($for_lon[$nforecast],$for_lat[$nforecast],$for_lon[$nforecast-1],$for_lat[$nforecast-1],$for_year[$nforecast],$for_nmonth[$nforecast],$for_day[$nforecast],$for_time[$nforecast],$for_year[$nforecast-1],$for_nmonth[$nforecast-1],$for_day[$nforecast-1],$for_time[$nforecast-1]);
                }
                #  Ending lon, lat (from the current record), converted to radians:
                #@end = (deg2rad($for_lon[$nforecast]), deg2rad(90 - $for_lat[$nforecast]));
                #  Great circle distance, in km:
                #$km = great_circle_distance(@start, @end, 6378.137);
                #  Bearing from start to end point (radians):
                #$dir = great_circle_bearing(@start, @end);
                #  Convert direction to degrees:
                #$tc_dir[$nforecast] = rad2deg($dir);
                #if ($Dd == 0) {
                #  $dtime = $Dh;
                #} else {
                #  $dtime = $Dh + ($Dd * 24);
                #}
                #  TC translation speed is distance / time:
                #$tc_speed[$nforecast] = ($km/$nm2km)/$dtime;

                if ($DEBUG1) {warn " ^^^^^^vvv FI = $fi; TC_Dir4 = $tc_dir[$nforecast]; TC_Speed = $tc_speed[$nforecast]; Delta_Time = $dtime hr, NForecast = $nforecast.\n";}
                @start = (); @end = (); $km = $INIT; $dir = $INIT; $dtime = $INIT;
            }

            #Original:         if ($DEBUG) {print "F/O $5 $6 $8 $9\n";}
            if ($DEBUG) {warn "F/O $nforecast: $for_lat[$nforecast] $for_latns[$nforecast] $for_lon[$nforecast] $for_lonew[$nforecast]\n";}

            #  Find radius of max. winds.
            $_ = <INPUT>;
            @line = split ' ', $_;

            #  Find the max. wind speed.
            if ($_ =~ /(MAX\s+WIND)/) {
                $for_winds[$fi] = $line[2];
                $for_rad[$fi] = get_radius($for_winds[$fi]) * $nm2km;
                #  Use both RMW algorithms to check whether either is better than the above.
                $rmw0 = 66.785 - (0.09102 * $for_winds[$fi]) + (1.0619 * (abs($for_lat[$nforecast]) - 25));  # Knaff & Zehr (2007) in km
                $rmw1 = (35.37 - (0.111 * $for_winds[$fi]) + (0.570 * (abs($for_lat[$nforecast]) - 25)))  * $nm2km;  # Gross et al. (2004) in nmi converted to km
                $rmw2 = 51.6 * exp((-0.0223 * $for_winds[$fi] * $kn2mps) + (0.0281 * abs($for_lat[$nforecast])));     # Willoughby & Rahn (2004) in km
                #warn "Forecast $fi: Radius = $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n";
                if ($DEBUG) { warn " Forecast $fi Wind, Radius = $for_winds[$fi], $for_rad[$fi]; RMW0 = $rmw0; RMW1 = $rmw1; RMW2 = $rmw2.\n"; }
            }  # end if($_
            $fi++;
        }

        #Original:      $nothing = 1;
        if ($for_year[$nforecast] ne '???') {
            if ($DEBUG) { warn "NForecast = $nforecast: Y,m,d,time,min = $for_year[$nforecast], $for_nmonth[$nforecast], $for_day[$nforecast], $for_time[$nforecast],$for_minute[$nforecast].\n"; }
            $datenum[$nforecast] = Mktime($for_year[$nforecast], $for_nmonth[$nforecast], $for_day[$nforecast], $for_time[$nforecast], $for_minute[$nforecast] + 0, 0.0);

        }

        #  Check for missing (i.e., zero) pressure & replace with the previous value or ambient pressure.
        if ($for_cp[$nforecast] <= 0) {
            #  Invalid value, so check for a previous valid value.
            if ($nforecast == 1) {
                #  First forecast, so check the analysis pressure.
                if ($cp > 0) {
                    #  Valid, so use it.
                    $for_cp[$nforecast] = $cp;
                } else {
                    #  Invalid, so use ambient.
                    $for_cp[$nforecast] = $pn;
                }
            } else {
                #  Past the 1st forecast, so check the previous value.
                if ($for_cp[$nforecast-1] > 0) {
                    #  Valid, so use it.
                    $for_cp[$nforecast] = $for_cp[$nforecast-1];
                } else {
                    #  Invalid, so use ambient.
                    $for_cp[$nforecast] = $pn;
                }
            }
        }

    }  # end of MAIN while loop...

    #  Calculate Holland B parameter based on current observations.  Warn the user
    #  if B is invalid, and reset the value to either the min or max valid value.
    #  Calculation is performed only if not appending to the output flag.
    if ($appendflag != 1 && $winds > 0 && $cp != $pn) { $B = ($rho  * $e  * (($winds * $kn2mps)**2)) / (abs($pn - $cp) * $mb2npm2); }
    if ($B > $BMax) {
        warn "  NOTICE: Holland B parameter exceeds $BMax; resetting to $BMax.\n";
        $B = $BMax;
    } elsif ($B < $BMin) {
        warn "  NOTICE: Holland B parameter less than $BMin; resetting to $BMin.\n";
        $B = $BMin;
    }

    if ($idump == 0) {
        if ($DEBUG) { warn "CP: $cp Lat: $latitude Lon: $longitude LonHem: $lonew RMW: $radius.\n"; }
        &DUMP();
    }

    if ($DEBUG) { warn "Datenum = @datenum\n"; }

    return;
}

#********************************************************************************
sub get_radius {
    #
    #  This subroutine obtains an averaged radius of maximum winds (RMW) based on present
    #  or forecast wind speeds and radii thereof for 4 quadrants.  The radii of the
    #  4 quadrants are averaged for each wind speed.  The present (or forecast) wind
    #  speeds and averaged radii are then used to extrapolate to the maximum sustained
    #  wind speed (or max. forecast wind speed) to obtain the RMW.  Extrapolation is
    #  performed using an external linear interpolation routine (cf. the 'require'
    #  statement above).  If only 1 forecast record with radii is available, the rmw
    #  is estimated using the ratio of the input maximum wind speed and the forecast
    #  wind speed corresponding to the radii.  Units are those of the advisory (currently,
    #  traditional nautical units -- i.e., speed in knots and radius in nautical miles).
    #  If Delft3D output is being generated, then the R34/50/64 radii values are interpolated
    #  to R35/50/65/100 values using an external subroutine (interp_rads).
    #
    #  Input(s):      $maxspd     -- Maximum wind speed
    #
    #  Output(s):     $rad        -- Radius of maximum wind speed
    #                 @r100_65_50_35-Radius of 100, 65, 50, 35 kt. winds, for Delft3D output
    #
    #  Called by:     process_data
    #
    #  Calls:  interp (interpl.pl), interp_rads (wes_routines.pl)
    #
    #  References:  Holland, G.J., J.I. Belanger, and A. Fritz (2010). A revised model for radial
    #               profiles of hurricane winds, Mon. Wea. Rev, 128, 4393-4401.
    #               Sampson, C.R., P.A. Wittmann, and H.L. Tolman (2010). Consistent tropical
    #               cyclone wind and wave forecasts for the U.S. Navy, Wea. Forecasting, 25,
    #               1293-1306.
    #
    #  Revision History:
    #  10 Apr 2006  Initial coding.  R.S. Linzell, PSI/Neptune Sciences Div.
    #  18 Apr 2006  Upgraded to ignore zero radius values when calculating the average
    #               radius, and to estimate the radius where only 1 forecast record has
    #               radii instead of ignoring such records.  (RSL)
    #  27 Jun 2006  Added code to check for invalid RMW values.  (RSL)
    #  02 Nov 2011  Added code to perform interpolation of R34/50/64 to R35/50/65/100
    #               in support of Delft3D output, including the call to sub. interp_rads,
    #               and added several new variables for this new capability; updated the
    #               help content.  (RSL)
    #  03 Nov 2011  Performed some bug fix on the Delft3D code.  (RSL)
    #  17 Feb 2012  Fixed a bug in the application of the wind speed ratio (changed the
    #               erroneous multiplication to division).  (RSL)
    #  10 May 2012  Added a nested IF block to always initialize the R100 arrays if
    #               Delft3D output is being produced; added parameters and code to
    #               store the wind radii in a 2D array (@dist2D) and to load those values
    #               in @100_65_50_35_3D; added some of the parameters in the my() commands
    #               at the beginning of this sub.  (RSL)
    #  05 Sep 2012  Changed the IF block after the main FOR loop, so that the output $rad is
    #               _always_ calculated based on the ratio of max. wind speed to max. of
    #               R34/R50/R64; changed debug 'print' statements to 'warn' statements;
    #               changed the IF block in which the value of $rad is compared to the
    #               value of $dmin so that if it is greater than or equal to $dmin, it
    #               is set to $dmin - 5 (nmi), per Sampson et al. (2010); modified the
    #               code for the RMW value ($rad) to use the $ratio squared (see comments).
    #               (RSL)
    #  13 Sep 2012  Commented-out the call to sub. interp_rads and the last FOR loop
    #               inside the last IF block for Delft3D output, so that interpolation
    #               of the R34/... wind radii is not performed.  (RSL)
    #  14 Sep 2012  Modified the IF statements within the FOR $k block to include any
    #               wind radii of zero, and added a WARN statement to notify the user
    #               that zero radii are included in the output.  (RSL)
    #  18 Sep 2012  Enclosed the new zero radii warning in an IF block.  (RSL)
    #  01 Feb 2013  Changed the 1st IF block to check abs($wes_output), to support
    #               best track output format; added an IF block to the last FOR loop
    #               to define $m3d according to the output type (WES or best track).
    #               (RSL)
    #  22 Feb 2013  Added an IF block to check for an empty @spd array, & to set the
    #               $dmin, $rad & $ratio values to zero if it is empty.  (RSL)
    #
    #********************************************************************************

    #  Input max. wind speed.
    my ($maxspd) = @_;

    #  Print statement for inclusion of zero wind radii:
    if ($DEBUG1) { warn "\n\n*******  NOTE: THE OUTPUT FILE INCLUDES WIND RADII OF ZERO!  *******\n\n"; }

    #  Local definitions and initialization:
    my $RMIN = 6/$nm2km;  # Minimum RMW in NM (6 km) per Deltares WES documentation (WES_Sci_tech_v8.doc, p. 7).
    my @line = ();
    my @avg_dist = ();
    my @spd = ();
    my ($ir100, @r100, @sorted_dist, @dist2D);
    my ($i, $j, $jj, $k, $m, $m3d, $num, $rad, $ratio, $sum, $dmin) = (0) x 11;
    my ($dist1, $dist2, $dist3, $dist4) = ($MISS) x 4;
    my $nquad = 4;

    if ($DEBUG1) { warn " -- -- in sub. get_radius; max spd = $maxspd, FI = $fi, NForecast = $nforecast...\n"; }

    #  Perform the initialization only if WES.EXE output (Delft3D) is needed.
    #  Initialize the R100 array with missing value codes.
    for $ir100 ( 0..$#spd_needed ) {
        $r100_65_50_35[$nforecast][$ir100] = $MISS;
        for ($i3d = 0; $i3d < $nquad; $i3d++) {
            $dist2D[$ir100][$i3d] = $MISS;
            $r100_65_50_35_3D[$nforecast][$ir100][$i3d] = $MISS;
        }
    }

    #  Loop over input records to find a list of wind speeds and radii of those speeds.
    for $m ( 0..10 ) {
        $_ = <INPUT>;
        #print "----Line $m: $_ -- Length = ",length($_),".\n";

        #  Check for initial blank line (usually has 1 space) if forecast wind speed not present.
        if ($m == 0 && length($_) <= 2) {
            #  Handle the case where Delft3D output is being generated; fill this
            #  row of the R100 array with missing values.
            return ($rad);
        }

        @line = split ' ', $_;
        if (length($_) > 2 && $line[0] =~ /\d+/ && $line[1] =~ /^KT/) {
            # get wind, radii;
            #warn "J = $j; RADII: @line\n";
            $spd[$j] = $line[0];
            for $k ( 0..$#line ) {
                #warn "  K + $k; Word = $line[$k]\n";

                if (uc($line[$k]) =~ /(NE$|NORTHEAST)/) {
                    ($dist1 = $line[$k]) =~ s/\D+//g;
                    #if ($dist1 > 0) { $dist2D[$j][0] = $dist1; }
                    if ($dist1 >= 0) { $dist2D[$j][0] = $dist1; }
                }
                if (uc($line[$k]) =~ /(SE|SOUTHEAST)/) {
                    $dist2 = substr($line[$k], 0, -2);
                    #if ($dist2 > 0) { $dist2D[$j][1] = $dist2; }
                    if ($dist2 >= 0) { $dist2D[$j][1] = $dist2; }
                }
                if (uc($line[$k]) =~ /(SW|SOUTHWEST)/) {
                    $dist3 = substr($line[$k], 0, -2);
                    #if ($dist3 > 0) { $dist2D[$j][2] = $dist3; }
                    if ($dist3 >=0) { $dist2D[$j][2] = $dist3; }
                }
                if (uc($line[$k]) =~ /(NW|NORTHWEST)/) {
                    $dist4 = substr($line[$k], 0, -3);
                    #if ($dist4 > 0) { $dist2D[$j][3] = $dist4; }
                    if ($dist4 >= 0) { $dist2D[$j][3] = $dist4; }
                }  # end if($line....
            }  # end for my $k...

            #  Get the sum & number of valid (i.e., nonzero) distances:
            $sum = 0;
            $num = 0;
            #if ($dist1 > 0) {
            if ($dist1 >= 0) {
                $sum += $dist1;
                $num++;
            }
            #if ($dist2 > 0) {
            if ($dist2 >= 0) {
                $sum += $dist2;
                $num++;
            }
            #if ($dist3 > 0) {
            if ($dist3 >= 0) {
                $sum += $dist3;
                $num++;
            }
            #if ($dist4 > 0) {
            if ($dist4 >= 0) {
                $sum += $dist4;
                $num++;
            }
            #Orig:         $avg_dist[$j] = ($dist1 + $dist2 + $dist3 + $dist4)/4;
            $avg_dist[$j] = $sum/$num;

            if ($DEBUG1) { warn " -- Spd $j: $spd[$j]; AvgDist = $avg_dist[$j] -- 4 dists: $dist1, $dist2, $dist3, $dist4.\n"; }

            $j++;
        } else {
            last;
        } # end if ($line...
    }  # end for my $m...

    #  Now we have a list of wind speeds and average radii of those speeds.  Extrapolate
    #  (using a linear interpolation routine) to determine the radius of the max. wind speed.
    #  Currently (as of Apr 2006), use only the 1st, last wind speeds & radii to extrapolate.
    #RSL      if ($#spd > 0) {
    #RSL        @spd = reverse @spd;
    #RSL        @avg_dist = reverse @avg_dist;
    #RSL      print "  about to call interp using s1 = $spd[$#spd-1], s2 = $spd[$#spd], sInt = $maxspd, d1 = $avg_dist[$#avg_dist-1], d2 = $avg_dist[$#avg_dist].\n";
    #RSL        $rad = interp($spd[$#spd-1], $spd[$#spd], $maxspd, $avg_dist[$#avg_dist-1], $avg_dist[$#avg_dist]);
    #RSL        print "  just returned from interp; RAD = $rad.\n";

    #$rad = abs(interp($spd[0], $spd[$#spd], $maxspd, $avg_dist[0], $avg_dist[$#avg_dist]));

    #  Find the minimum average distance of this data set by simply sorting the values
    #  and selecting the 1st value.
    #RSL        @sorted_dist = sort {$a <=> $b} @avg_dist;
    #RSL        $dmin = $sorted_dist[0];
    #RSL      } elsif ($#spd == 0) {

    #  If only 1 valid forecast record (including radii & wind speeds), use the ratio of the
    #  max. forecast wind speed to the wind speed w/radii to estimate the radius of max. winds.
    #  Note that the ratio is squared; cf. Holland et al. (2010), p. 4394.
    if (scalar(@spd) > 0) {
        $ratio = $maxspd/$spd[0];
        $rad = $avg_dist[0] / ($ratio*$ratio);
        $dmin = $avg_dist[0];
    } else {
        #return(0);
        $ratio = 0;
        $rad = 0;
        $dmin = 0;
    }
    #RSL      }

    #  Gross error check
    if ($rad < $RMIN) {
        if ($DEBUG1) { warn "Calculated RMW ($rad) is less than the minimum ($RMIN); setting to the minimum.\n"; }
        $rad = $RMIN;
    } elsif ($rad >= $dmin) {
        #  If calculated RMAX is greater than the (average) minimum of the R34/R50/R64 radii, then it is set
        #  to the minimum - 5 nmi (Sampson et al., 2010, p. 1296).
        if ($DEBUG1) { warn "Calculated RMW ($rad) is greater than the (average) minimum wind radius for this data set ($dmin); setting to the minimum - 5 nmi.\n"; }
        $rad = $dmin - 5;
    } elsif ($rad > $RMW_MAX) {
        if ($DEBUG1) { warn "  WARNING: RMW > Maximum \($RMW_MAX km\); resetting to Maximum.\n"; }
        $rad = $RMW_MAX;
    }

    if ($DEBUG1) { warn "  -- Final radius: $rad.\n"; }
    #  Interpolate from NHC wind radii to WES.EXE radii if Delft3D output is being
    #  generated.

    #  Loop over the number of wind radii to store the radii values in the output array.
    for $m (0..$#dist2D) {
        #  Calculate the wind value (i.e., R100 is 0, R64/65 is 1, ...) index of the 3D array
        #  since the values are written in that order vice the order in the input file.
        $m3d = $#dist2D - $m;
        #print " * * * * * * M = $m * * * M3D = $m3d* * * * * *Dist2D = $dist2D[$m][0] $dist2D[$m][1] $dist2D[$m][2] $dist2D[$m][3]\n";
        $r100_65_50_35_3D[$nforecast][$m3d][0] = $dist2D[$m][0];
        $r100_65_50_35_3D[$nforecast][$m3d][1] = $dist2D[$m][1];
        $r100_65_50_35_3D[$nforecast][$m3d][2] = $dist2D[$m][2];
        $r100_65_50_35_3D[$nforecast][$m3d][3] = $dist2D[$m][3];
    }

    return($rad);
}

#********************************************************************************
sub ISNUMBER {
    #********************************************************************************

    # Is the argument a number, e.g. 1.0, .1, 1, 1e4, 1.0e4, etc.?
    # Returns true or false.

    $_[0] =~ /^\s*(\+|-|)(\d+|\d+\.\d*|\d*\.\d+)(e\d+|)\s*$/i;
}

#********************************************************************************
sub RESET {
    #********************************************************************************

    # Reset global variables

    if ($DEBUG) {print "RESET\n";}

    # The following hold data for the actual observation
    $time="???"; $minute="00"; $dow="???"; $month="???"; $nmonth="???"; $day="???"; $year="???";
    $latitude="???"; $longitude="???"; $latns="???"; $lonew="???";
    $course="???"; $speed="???";
    $pressure="???";
    $winds="???"; $gusts="???";
    $name="???"; $type="???";
    $wmoheader="???";

    @for_name="???";        # The for_ arrays hold the data for the forecasts
    @for_time="???";
    @for_minute="00";
    @for_nmonth="???";
    @for_day="???";
    @for_year="???";
    @for_year2="???";
    @for_lat="???";
    @for_lon="???";
    @for_latns="???";
    @for_lonew="???";
    @for_winds="???";
    @for_gusts="???";
    $nforecast=0;
    $forecast_active = 0;  # Flag to check whether data is a FORECAST

    $abcmessage = 0;  # Flag to check if message is of the A., B. etc. type
    $exercise=0;      # Flag to check if message is an exercise
}

#********************************************************************************
sub DUMP {
    #********************************************************************************

    # Check that name is no in the bad list
    $name =~ s/\(//;
    $name =~ s/\)//;
    $name = "\U$name";
    $metaname = $name;
    $metaname =~ s/(\W)/\\\1/g;     # make any metacharacters literal
    $nbadname = grep(/^$metaname$/i,@badnames);
    my ($num, $sum, $ratio) = (0) x 3;
    $MeanRMW = (20/2) * $nm2km;	# Mean eye diameter is ~20 nmi (http://www.srh.weather.gov/hgx/tropical/geninfo.htm)

    # Dump the actual data to STDOUT
    if (&ISNUMBER($year) && &ISNUMBER($nmonth) && &ISNUMBER($day) &&
        &ISNUMBER($time) &&
        &ISNUMBER($latitude) && &ISNUMBER($longitude) &&
        $latns ne "???" && $latns=~/\w+/ &&
        $lonew ne "???" && $lonew=~/\w+/ &&
        $name ne "???" && $name=~/\w+/ && $nbadname<=0 &&
        $name ne "EXERCISE" &&
        $wmoheader ne "???" &&
        $type ne "???" && $type=~/\w+/ && !$exercise) {
        $doy = $dayyear[$nmonth-1]+$day;
        #Original:      if ($nmonth > 2 && ($year % 4) == 0) {$doy++;}  # leap year
        if ($nmonth > 2 && isleap($year)) {$doy++;}  # leap year
        if (!&ISNUMBER($minute)) {$minute="00";}
        if ($minute < 0) {$minute = -$minute;}
        if ($minute > 60) {$minute = "00";}
        if (isleap($year)) {
            $ddate = $year + ($doy + $time/24.0 + $minute/1440.0)/366.0;
        } else {
            $ddate = $year + ($doy + $time/24.0 + $minute/1440.0)/365.0;
        }
        if (&ISDEPRESSION($name)) {$name = &CHECKDEPNAME($name);}
        if (!&ISNUMBER($winds)) {$winds = "???";}
        if (!&ISNUMBER($gusts)) {$gusts = "???";}
        if (&ISNUMBER($winds) && &ISNUMBER($gusts) && $gusts < $winds) {$gusts = "???";}
        if (!&ISNUMBER($speed)) {$speed = "???";}
        if (!&ISNUMBER($course)) {$course = "???";}
        $outstr = "$year $nmonth $day $time.$minute UT $latitude $latns $longitude $lonew ${name}-$year2 $course T $speed kt $pressure mb $winds $gusts kt $type ACT $ddate $wmoheader $now";
        $outstr =~ s/  / /g;
        # print only if the data is older than now plus twleve hours
        if ($DEBUG && $now+12.1*0.000113843 > $ddate) {
            warn "$outstr\n";
        }

        # Now dump the forecast data to STDOUT.  Note: No forecasts are output if
        # there was no actual observation.  This is a safety feature to ensure
        # that garbled messages are not output.

        if ($#for_year < $#for_winds) {  # How many forecasts?
            $nforecast=$#for_year;
        } else {
            $nforecast=$#for_winds;
        }
        for ($i=1; $i<=$nforecast; $i++) {
            if (&ISNUMBER($for_year[$i]) && &ISNUMBER($for_nmonth[$i]) &&
                &ISNUMBER($for_day[$i]) && &ISNUMBER($for_time[$i]) &&
                &ISNUMBER($for_lat[$i]) && &ISNUMBER($for_lon[$i]) &&
                $for_latns[$i] ne "???" && $for_latns[$i]=~/\w+/ &&
                $for_lonew[$i] ne "???" && $for_lonew[$i]=~/\w+/ &&
                &ISNUMBER($for_winds[$i]) &&
                $for_name[$i] ne "???" && $for_name[$i]=~/\w+/ && !$exercise) {
                if ($DEBUG) {warn "dumping forecast $i\n";}
                $for_name[$i] =~ s/\(//;
                $for_name[$i] =~ s/\)//;
                if (&ISDEPRESSION($for_name[$i])) {$for_name[$i] = &CHECKDEPNAME($for_name[$i]);}
                $foryear = $for_year[$i];
                $foryear2 = $for_year2[$i];
                $fornmonth = $for_nmonth[$i];
                $forday = $for_day[$i];
                $fortime = $for_time[$i];
                $forminute = $for_minute[$i];
                if (&ISNUMBER($forminute)) { $forminute = $for_minute[$i];}
                else { $forminute = "00";}
                if ($forminute < 0) {$forminute = -$forminute;}
                if ($forminute > 60) {$forminute = "00";}
                $forlat = $for_lat[$i];
                $forlon = $for_lon[$i];
                $forlatns = $for_latns[$i];
                $forlonew = $for_lonew[$i];
                $forwinds = $for_winds[$i];
                $forgusts = $for_gusts[$i];
                $forname = $for_name[$i];
                $fordoy = $dayyear[$fornmonth-1]+$forday;
                #Original:            if ($fornmonth > 2 && ($foryear % 4) == 0) {$fordoy++;}  # leap year
                if ($fornmonth > 2 && isleap($foryear)) {$fordoy++;}  # leap year
                if (isleap($foryear)) {
                    $forddate[$i] = $foryear + ($fordoy + $fortime/24.0 + $forminute/1440.0)/366.0;
                } else {
                    $forddate[$i] = $foryear + ($fordoy + $fortime/24.0 + $forminute/1440.0)/365.0;
                }
                if (!&ISNUMBER($forwinds)) {$forwinds = "???";}
                if (!&ISNUMBER($forgusts)) {$forgusts = "???";}
                if (&ISNUMBER($forwinds) && &ISNUMBER($forgusts) && $forgusts < $forwinds) {$forgusts = "???";}
                $outstr = "$foryear $fornmonth $forday $fortime.$forminute UT $forlat $forlatns $forlon $forlonew ${forname}-$foryear2 ??? T ??? kt ??? mb $forwinds $forgusts kt ??? FOR $forddate[$i] $wmoheader $ddate";
                $outstr =~ s/  / /g;
                # Print if the forecast is no more than 7 days + 12 h = 180 hours after now
                if ($DEBUG && $now+180.1*0.000113843 > $forddate[$i]) {
                    warn "$outstr\n";
                }
            } else {
                if ($DEBUG) {warn "Can't dump forecast $i\n";}
            }
        }

        #  Calculate central pressures (CPs) at forecast periods based on wind speeds & decimal dates.
        if ($nforecast >= 2) {
            if ($DEBUG) {warn "Estimating forecast CPs...\n"; }

            for ($i=1; $i<=$nforecast; $i++) {

                #  Ratio of present wind speed to this forecast wind speed:
                $ratio = $winds/$for_winds[$i];

                #  If the present RMW wasn't found, estimate from the 1st forecast value.
                if ($radius == 0 && $i == 1) { $radius = $for_rad[$i] * $ratio; }

                #  This forecast central pressure is simply the present CP * the wind speed ratio.
                #$for_cp[$i] = $cp * $ratio;

                #  Upgrade of 28 Jul 2011: central pressure is calculated per Dvorak(1984) based
                #  on max. steady wind (MSW).
                $for_cp[$i] = wind_pressure_relationship($for_winds[$i],$BINDEX);

                #  Perform gross error check on forecast CP.
                if ($for_cp[$i] < $CP_MIN) {
                    warn "  WARNING: Forecast Central Pressure \($for_cp[$i]\) < Minimum \($CP_MIN mb\); resetting to Minimum.\n";
                    $for_cp[$i] = $CP_MIN;
                } elsif ($for_cp[$i] > $pn) {
                    warn "  WARNING: Forecast Central Pressure \($for_cp[$i]\) > Ambient \($pn mb\); resetting to Ambient.\n";
                    $for_cp[$i] = $pn;
                }

                ##print "DDate: $ddate CP: $cp  WInds: $winds For_Wind:  $for_winds[$i]  Ratio: $ratio:  For_CP: $for_cp[$i].\n";
                ##          if ($i < $nforecast) {

                ##print "I: $i  Next ForDdate: $forddate[($i+1)]  Next For_Winds: $for_winds[($i+1)].\n";

                #Original:            $forcp[$i] = interp($ddate, $forddate[($i+1)], $cp, $winds, $for_winds[($i+1)]);
                ##            $for_cp[$i] = interp($forddate[($i+1)], $ddate, $cp, $for_winds[($i+1)], $winds);
                ##          } else {

                ##print "I: $i  Last ForDdate: $forddate[($i-1)]  Last For_Winds: $for_winds[($i-1)].\n";

                ##            $for_cp[$i] = interp($ddate, $forddate[($i-1)], $cp, $winds, $for_winds[($i-1)]);
                ##          }

                #  Sum the forecast RMW to get an average.
                #  Start w/the analysis value.
                if ($i == 1) {
                    $sum = $radius;
                    $num = 1;
                }
                if ($for_rad[$i] > 0) {
                    $sum += $for_rad[$i];
                    $num++;
                }
            } # end for($i...
            #  Define mean of forecast RMWs.
            if ($num > 0) {
                $MeanRMW = $sum/$num;
                if ($DEBUG) { warn ".............. MeanRMW = $MeanRMW\n"; }
            }
            ##      } else {
            ##        $for_cp[0] = 0;
        } # end if($nforecast....

    } else {
        if ($DEBUG) {
            warn "dump failed:\n";
            warn "year=$year nmonth=$nmonth day=$day time=$time\n";
            warn "lat=$latitude lon=$longitude latns=$latns lonew=$lonew\n";
            warn "name=$name wmo=$wmoheader type=$type exercise=$exercise\n";
        }
    }
}

#********************************************************************************
sub GETNMONTH {
    #********************************************************************************
    local($nowyear,$nowmonth,$nowday,$day) = @_;

    # Set the month number given the day of month and the current year,
    # current month, and current day.  Returns adjusted year and month.

    # $nowyear is like 1994
    # $nowmonth is like 8 (for August)
    # $nowday is like 31
    # $day is like 31

    # Use same month as now if within +7 or -20 days from now.  Not symmetric
    # because there is not likely to be a time far in the future (the
    # longest forecast should be 3 days), but there is likely to be a time
    # in the past (old observation).

    local($year,$month) = ($nowyear,$nowmonth);
    if (($nowday-$day) > 20) {  # Is a month border crossed forwards?
        $month++;
    }
    if ($month > 12) {
        $year++;
        $month=1;
    }
    if (($day-$nowday) > 7) {  # Is a month border crossed backwards?
        $month--;
    }
    if ($month < 1) {
        $year--;
        $month=12;
    }
    if ($DEBUG) {warn "GETNMONTH: ",$nowyear," ",$nowmonth," ",$nowday," ",$year," ",$month," ",$day,"\n";}

    return($year,$month);   # return values

}

#********************************************************************************
sub ISDEPRESSION {
    #********************************************************************************

    # Check if a name looks like a tropical depression name

    local($newname) = @_;
    local($junk) = 0;
    local($HyphenNumNew) = 0;

    $HyphenNumNew = split(/\-/,$newname);   # How many hyphens in the new name?
    if ($HyphenNumNew > 2 ||
        $newname =~ /\d+/ ||
        $newname =~ /[^a-z0-9 -]/i ||
        $newname =~ /^\d+[A-Z|]/i ||
        $newname =~ /^one$/i ||
        $newname =~ /^one\-/i ||
        $newname =~ /^two$/i ||
        $newname =~ /^two\-/i ||
        $newname =~ /^three$/i ||
        $newname =~ /^three\-/i ||
        $newname =~ /^four$/i ||
        $newname =~ /^four\-/i ||
        $newname =~ /^five$/i ||
        $newname =~ /^five\-/i ||
        $newname =~ /^six$/i ||
        $newname =~ /^six\-/i ||
        $newname =~ /^seven$/i ||
        $newname =~ /^seven\-/i ||
        $newname =~ /^eight$/i ||
        $newname =~ /^eight\-/i ||
        $newname =~ /^nine$/i ||
        $newname =~ /^nine\-/i ||
        $newname =~ /^ten$/i ||
        $newname =~ /^ten\-/i ||
        $newname =~ /^eleven$/i ||
        $newname =~ /^eleven\-/i ||
        $newname =~ /^twelve$/i ||
        $newname =~ /^twelve\-/i ||
        $newname =~ /^thirteen$/i ||
        $newname =~ /^thirteen\-/i ||
        $newname =~ /^fourteen$/i ||
        $newname =~ /^fourteen\-/i ||
        $newname =~ /^fifteen$/i ||
        $newname =~ /^fifteen\-/i ||
        $newname =~ /^sixteen$/i ||
        $newname =~ /^sixteen\-/i ||
        $newname =~ /^seventeen$/i ||
        $newname =~ /^seventeen\-/i ||
        $newname =~ /^eighteen$/i ||
        $newname =~ /^eighteen\-/i ||
        $newname =~ /^nineteen$/i ||
        $newname =~ /^nineteen\-/i ||
        $newname =~ /^twenty\-/i ||      # is there a dash after twenty etc.?
        $newname =~ /^twenty$/i ||
        $newname =~ /^thirty\-/i ||
        $newname =~ /^thirty$/i ||
        $newname =~ /^forty\-/i ||
        $newname =~ /^forty$/i ||
        $newname =~ /^forty\-/i ||
        $newname =~ /^forty$/i ||
        $newname =~ /^fifty\-/i ||
        $newname =~ /^fifty$/i ||
        $newname =~ /^sixty\-/i ||
        $newname =~ /^sixty$/i ||
        $newname =~ /^seventy\-/i ||
        $newname =~ /^seventy$/i ||
        $newname =~ /^eighty\-/i ||
        $newname =~ /^eighty$/i ||
        $newname =~ /^ninety\-/i ||
        $newname =~ /^ninety$/i) {$junk=1;} else {$junk=0;}
    $junk;
}

#********************************************************************************
sub CHECKDEPNAME {
    #********************************************************************************
    #Original:         local($name) = @_[0];
    local($name) = $_[0];
    $name =~ s/^ONE\-/01/i;
    $name =~ s/^ONE/01/i;
    $name =~ s/^TWO\-/02/i;
    $name =~ s/^TWO/02/i;
    $name =~ s/^THREE\-/03/i;
    $name =~ s/^THREE/03/i;
    $name =~ s/^FIVE\-/05/i;
    $name =~ s/^FIVE/05/i;
    $name =~ s/^TEN\-/10/i;
    $name =~ s/^TEN/10/i;
    $name =~ s/^ELEVEN\-/11/i;
    $name =~ s/^ELEVEN/11/i;
    $name =~ s/^TWELVE\-/12/i;
    $name =~ s/^TWELVE/12/i;
    $name =~ s/^THIRTEEN\-/13/i;
    $name =~ s/^THIRTEEN/13/i;
    $name =~ s/^FOURTEEN\-/14/i;
    $name =~ s/^FOURTEEN/14/i;
    $name =~ s/^FOUR\-/04/i;
    $name =~ s/^FOUR/04/i;
    $name =~ s/^FIFTEEN\-/15/i;
    $name =~ s/^FIFTEEN/15/i;
    $name =~ s/^SIXTEEN\-/16/i;
    $name =~ s/^SIXTEEN/16/i;
    $name =~ s/^SIX\-/06/i;
    $name =~ s/^SIX/06/i;
    $name =~ s/^SEVENTEEN\-/17/i;
    $name =~ s/^SEVENTEEN/17/i;
    $name =~ s/^SEVEN\-/07/i;
    $name =~ s/^SEVEN/07/i;
    $name =~ s/^EIGHTEEN\-/18/i;
    $name =~ s/^EIGHTEEN/18/i;
    $name =~ s/^EIGHT\-/08/i;
    $name =~ s/^EIGHT/08/i;
    $name =~ s/^NINETEEN\-/19/i;
    $name =~ s/^NINETEEN/19/i;
    $name =~ s/^NINE\-/09/i;
    $name =~ s/^NINE/09/i;
    $name =~ s/^TWENTY\-ONE\-/21/i;
    $name =~ s/^TWENTY\-ONE/21/i;
    $name =~ s/^TWENTY\-TWO\-/22/i;
    $name =~ s/^TWENTY\-TWO/22/i;
    $name =~ s/^TWENTY\-THREE\-/23/i;
    $name =~ s/^TWENTY\-THREE/23/i;
    $name =~ s/^TWENTY\-FOUR\-/24/i;
    $name =~ s/^TWENTY\-FOUR/24/i;
    $name =~ s/^TWENTY\-FIVE\-/25/i;
    $name =~ s/^TWENTY\-FIVE/25/i;
    $name =~ s/^TWENTY\-SIX\-/26/i;
    $name =~ s/^TWENTY\-SIX/26/i;
    $name =~ s/^TWENTY\-SEVEN\-/27/i;
    $name =~ s/^TWENTY\-SEVEN/27/i;
    $name =~ s/^TWENTY\-EIGHT\-/28/i;
    $name =~ s/^TWENTY\-EIGHT/28/i;
    $name =~ s/^TWENTY\-NINE\-/29/i;
    $name =~ s/^TWENTY\-NINE/29/i;
    $name =~ s/^TWENTY\-/20/i;
    $name =~ s/^TWENTY/20/i;
    $name =~ s/^THIRTY\-ONE\-/31/i;
    $name =~ s/^THIRTY\-ONE/31/i;
    $name =~ s/^THIRTY\-TWO\-/32/i;
    $name =~ s/^THIRTY\-TWO/32/i;
    $name =~ s/^THIRTY\-THREE\-/33/i;
    $name =~ s/^THIRTY\-THREE/33/i;
    $name =~ s/^THIRTY\-FOUR\-/34/i;
    $name =~ s/^THIRTY\-FOUR/34/i;
    $name =~ s/^THIRTY\-FIVE\-/35/i;
    $name =~ s/^THIRTY\-FIVE/35/i;
    $name =~ s/^THIRTY\-SIX\-/36/i;
    $name =~ s/^THIRTY\-SIX/36/i;
    $name =~ s/^THIRTY\-SEVEN\-/37/i;
    $name =~ s/^THIRTY\-SEVEN/37/i;
    $name =~ s/^THIRTY\-EIGHT\-/38/i;
    $name =~ s/^THIRTY\-EIGHT/38/i;
    $name =~ s/^THIRTY\-NINE\-/39/i;
    $name =~ s/^THIRTY\-NINE/39/i;
    $name =~ s/^THIRTY\-/30/i;
    $name =~ s/^THIRTY/30/i;

    $name;
}

#********************************************************************************
sub isleap {
    #
    #  This subroutine determines whether an input year is a leap year.  See the URL
    #  http://aa.usno.navy.mil/faq/docs/leap_years.html for an explanation.  Note that
    #  the input year should be 4 digits in length (e.g., 2006 instead of 06).
    #
    #  Input(s):      $testyr     -- Year to be checked
    #
    #  Output(s):     $ans        -- Indicator of leap year (=1) or not (=0)
    #
    #  Called by:     DUMP
    #
    #********************************************************************************
    my ($testyr) = @_;
    my $ans = 0;

    if (($testyr % 4) == 0) {
        $ans = 1;
        if (substr($testyr, -2) eq "00" && ($testyr % 400) != 0) {
            $ans = 0;
        }
    } else {
        $ans = 0;
    }
    return($ans);
}

#*******************************************************************************
sub wind_pressure_relationship {
    #
    #  This subroutine calculates a minimum sea level pressure (MSLP) given a maximum
    #  steady wind (MSW) speed using the Dvorak WestPac wind-pressure relationship.
    #
    #  References:
    #  Dvorak, V.F., 1975: Tropical cyclone intensity analysis using satellite data.
    #    NOAA Tech. Rep. NESDIS 11, 45 pp. [presently not available on-line]
    #  Velden, C., B. Harper, F. Wells, J.L. Beven II, R. Zehr, T. Olander, M. Mayfield,
    #    C. Guard, M. Lander, R. Edson, L. Avila, A. Burton, M. Turk, A. Kikuchi, A.
    #    Christian, P. Caroff, and P. McCrone, 2006: The Dvorak tropical cyclone intensity
    #    estimation technique: A satellite-based method that has entured for over 30 years,
    #    Bull. Amer. Meteor. Soc., 87, 1195-1210.
    #
    #  Syntax: $mslp = wind_pressure_relationship($msw,$bindex);
    #  where:  $mslp is the returned MSLP value (in hPa),
    #          $msw is the user-supplied MSW value (in kt), and
    #          $bindex is the basin index (1 = Atlantic, 2 = Pacific)
    #
    #  Calls:  interp (subroutine)
    #
    #  Called by: main program
    #
    #  Revision History:
    #  18 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  22 Jul 2011  Added an ELSIF block to extrapolate to wind speeds below the
    #               minimum wind speed of the table.  (RSL)
    #  28 Jul 2011  Copied from parse_jtwc_warning.pl; added middle column of @WPR
    #               for the Atlantic basin; added input of $bindex to indicate which
    #               basin column to use.  (RSL)
    #  05 Sep 2012  Added 'last' statements in the IF block where appropriate.  (RSL)
    #
    #*******************************************************************************

    my ($msw,$bindex) = @_;
    my $mslp = 0;

    #  Define the wind-pressure relationship for the Atlantic & WestPac basins,
    #  based on Dvorak (1984) as summarized by Velden et al. (2006).
    #  The 1st column is maximum steady wind speed (MSW) in knots;
    #  the 2nd column is minimum sea level pressure (MSLP) in hPa (= mb) for the
    #  Atlantic basin, and the 3rd column is MSLP for the WestPac basin.
    my @WPR = (
        [30, 1009, 1000],
        [35, 1005,  997],
        [45, 1000,  991],
        [55,  994,  984],
        [65,  987,  976],
        [77,  979,  966],
        [90,  970,  954],
        [102,  960,  941],
        [115,  948,  927],
        [127,  935,  914],
        [140,  921,  898],
        [155,  906,  879],
        [170,  890,  858],
    );

    #  Loop over wind speeds....
    for my $i (0..$#WPR) {
        #  Check whether this MSW value equals the user-supplied value.
        if ($WPR[$i][0] == $msw) {
            #  The values are equal.  Simply return the corresponding MSLP.
            $mslp = $WPR[$i][$bindex];
            last;
        } elsif ($i > 0 && ($WPR[$i][0] > $msw && $WPR[$i-1][0] < $msw) ||
            ($i == $#WPR && $WPR[$i][0] < $msw)) {
            #  Interpolate between wind speeds, or extrapolate above the max. wind speed.
            $mslp = interp($WPR[$i][0],$WPR[$i-1][0],$msw,$WPR[$i][$bindex],$WPR[$i-1][$bindex]);
            last;
        } elsif ($i == 0 && $WPR[$i][0] > $msw) {
            #  Extrapolate to low wind speed below the minimum.
            $mslp = interp($WPR[$i+1][0],$WPR[$i][0],$msw,$WPR[$i+1][$bindex],$WPR[$i][$bindex]);
            last;
        } else {
            #  Some sort of fail-safe may go here.
        }
    }  # end for my $i....

    return $mslp;
}

#*******************************************************************************
sub storm_velocity {
    #
    #  This subroutine calculates the velocity of the storm itself based on the current
    #  and previous positions, dates & times.
    #
    #  Syntax: ($speed, $dir) = storm_velocity($this_lon,$this_lat,$last_lon,$last_lat,
    #          $yy,$mm,$dd,$hh,$yy_last,$mm_last,$dd_last,$hh_last)
    #  where:  $speed is the returned storm translational speed,
    #          $dir is the returned storm direction,
    #          $this_lon is the current longitude,
    #          $this_lat is the current latitude,
    #          $last_lon is the previous longitude,
    #          $last_lat is the previous latitude,
    #          $yy is the current year,
    #          $mm is the current month,
    #          $dd is the current day,
    #          $hh is the current hour,
    #          $yy_last is the previous year,
    #          $mm_last is the previous month,
    #          $dd_last is the previous day, and
    #          $hh_last is the previous hour.
    #
    #  Calls:  Date::Calc::Delta_DHMS, Math::Trig::(great_circle_distance deg2rad
    #          rad2deg great_circle_bearing)
    #
    #  Called by: sub. process_data
    #
    #  Revision History:
    #  04 Sep 2013  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  05 Sep 2013  Wrapped up initial coding; added help content; renamed from
    #               calc_speed_dir to storm_velocity.  (RSL)
    #
    #*******************************************************************************

    my ($this_lon,$this_lat,$last_lon,$last_lat,$yy,$mm,$dd,$hh,$yy_last,$mm_last,$dd_last,$hh_last) = @_;
    my(@start,@end,$dir,$km,$dtime,$Dd,$Dh,$Dm,$Ds,$speed);

    #  Import the great circle & other trig. functions, date calculation function.
    use Date::Calc qw(Delta_DHMS);
    use Math::Trig qw(great_circle_distance deg2rad rad2deg great_circle_bearing);

    #  Starting lon, lat (from previous record) converted to radians:
    @start = (deg2rad($last_lon), deg2rad(90 - $last_lat,));
    #  Calculate time difference between the current and previous records.
    ($Dd,$Dh,$Dm,$Ds) = Delta_DHMS($yy_last,$mm_last,$dd_last,$hh_last,0,0,$yy,$mm,$dd,$hh,0,0);
    #  Ending lon, lat (from the current record), converted to radians:
    @end = (deg2rad($this_lon), deg2rad(90 - $this_lat));
    #  Great circle distance, in km:
    $km = great_circle_distance(@start, @end, 6378.137);
    #  Bearing from start to end point (radians), converted to degrees:
    $dir = rad2deg(great_circle_bearing(@start, @end));
    #  Calculate delta time.
    if ($Dd == 0) {
        $dtime = $Dh;
    } else {
        $dtime = $Dh + ($Dd * 24);
    }
    #  TC translation speed is distance / time:
    $speed = ($km/$nm2km)/$dtime;

    return($speed, $dir);
}
