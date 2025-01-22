#!/usr/bin/perl -w
#*******************************************************************************
#
#  parse_jtwc_warning.pl
#
#  This Perl script reads a tropical cyclone (TC) warning file from the Joint Typhoon
#  Warning Center (JTWC) and extracts the pertinent information for use in an
#  external meteorological model.  Output is in one of two formats: if only 2
#  input parameters are specified, then output is a single file in "cyclone.dat"
#  format for use with the PCTides program, cyclone.exe.  If six input parameters
#  are specified, then two output files are created, one in wes.inp format, and one
#  in trackfile.trk format for use in the Delft3D program, WES.EXE.  If seven input
#  parameters are specified, then a JTWC best track must be the 7th parameter, and
#  the data from it are prepended to the forecast values from the TC warning file.
#
#  NB: The external Perl module, Date::Calc, is required.  This module is NOT generally
#      included in a standard Perl installation, so it may require separate installation.
#      Such installation may require administrator privileges.  Another external Perl
#      module, Math::Trig, also is required.  This module IS generally included in a
#      standard Perl installation.  The presence of these modules can be checked by
#      entering this command from a command prompt: perl -cW parse_jtwc_warning.pl
#      and noting the output.
#
#  Syntax: parse_jtwc_warning.pl infile outfile [best_track_file] --OR--
#          parse_jtwc_warning.pl infile outfile [trkfile num_rad_bins num_dir_bins tc_radius [best_track_file]]
#  where:  infile is the input ASCII text file with the warning text,
#          outfile is the output ASCII text file with the needed TC parameters,
#          trkfile is the optional output ASCII text track file for WES.EXE,
#          num_rad_bins is the number of Delft3D Spiderweb file radial bins (required if trkfile is specified),
#          num_dir_bins is the number of Delft3D Spiderweb file directional bins (required if trkfile is specified),
#          tc_radius is the tropical cyclone radius (km) (required if trkfile is specified), and
#          best_track_file is an optional JTWC best track file corresponding to the infile.
#
#  NOTE: The NRL extended best track output option is enabled using the 2nd syntax above, in which
#        the 'trkfile' parameter has a name that does not have a file extension of ".trk" -- the
#        NRL convention is to use a file extension of ".ebt" for the extended best track files.
#        See the best_track_routines.pl source code for further details.
#
#  Revision History:
#  18 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
#  28 Oct 2011  Began adding support for the Delft3D program, WES.EXE (i.e., added
#               creation of the *.inp and *.trk files for that program).  (RSL)
#  31 Oct 2011  Added arrays for TC direction & speed, which are needed for the WES.EXE
#               track file; upgraded sub. read_jtwc_file to determine the direction
#               & speed; added sub. interp_rads to perform interpolation of wind
#               radii from R64/R50/R34 to R100/R65/R50/R35, and added a call to that
#               sub. to sub. read_jtwc_file; upgraded sub. print_track_file to write
#               the actual values.  (RSL)
#  01 Nov 2011  Fixed a bug in date determination in sub. read_jtwc_file; modified sub.
#               print_track_file to build the output format string on the fly based
#               on the values of @r100_65_50_35, and removed the printing of TC
#               direction & speed, as WES.EXE does not expect these values for
#               Method 3; updated the help content; added further validity checks
#               of the R100 values to subs. read_jtwc_file & interp_rads;
#               removed sub. interp and added the 'require' statement near the top of
#               the code; removed the WES.EXE-specific subs., print_wes_inp_file,
#               print_track_file and interp_rads, put them in wes_routines.pl, and
#               added the corresponding 'require' statement; modified the call to
#               sub. print_track_file to include the needed array references vice
#               using the global variables.  (RSL)
#  16 Dec 2011  Added support to extract the storm name in $tc_name, and print it
#               in the .trk file if that option is used.  (RSL)
#  06 Feb 2012  Changed the call to die() to a call to warn() in sub. read__jtwc_file
#               when an TC formation alert is being processed, as such an alert may
#               still have TC product information.  (RSL)
#  09 Feb 2012  Added sub. check_longitude to check the signs of the longitude values
#               and to ensure that all have the same sign; added a call to sub. check_
#               longitude to the main program immediately after the call to sub.
#               read_jtwc_file.  (RSL)
#  24 Feb 2012  Added sub. read_best_track (copied from tst_parse_best.pl); began
#               adding support for reading JTWC best track files and combining
#               those parameters with the forecast parameters from a "standard" JTWC
#               TC warning message; updated the help content and usage message.  (RSL)
#  02 Mar 2012  Upgraded sub. wind_pressure_relationship: modified the @WPR array
#               to use the values provided by JTWC personnel (see the comments in
#               that subroutine for details).  (RSL)
#  12 Mar 2012  Added code in sub. write_tc_parameters to correct negative longitudes
#               for PCTides output.  (RSL)
#  03 May 2012  Upgraded sub. read_best_track to read, store the R34/R50/R64
#               parameters and correctly handle the duplicate records present
#               when R50/R64 parameters are present; upgraded the main program
#               to use the new call syntax for sub. read_best_track & to print a
#               table of values after returning from sub. read_best_track, and to
#               build a hash with the date/time groups for the different data sets
#               to provide a means of sorting the data records chronologically;
#               upgraded sub. write_tc_parameters to add support for best track
#               file output.  (RSL)
#  04 May 2012  Added an IF block to sub. write_tc_parameters to print the best
#               track environmental pressure if it is available; changed the type
#               descriptor from 'WARN' to 'W' and 'BEST' to 'B'; added new sub.
#               combine_data to combine best track and warning file data into the
#               original warning data arrays in chronological order; added a call
#               to sub. combine_data to the main program in the IF block in which
#               sub. read_best_track is called; removed the IF block in sub.
#               write_tc_parameters in the main FOR loop, since all data are now
#               combined in the original data arrays; moved the 'use Math::Trig'
#               pragma from sub. read_jtwc_file to the main program; move the code
#               to calculate storm translational speed & direction from sub.
#               read_jtwc_file to newly added sub. storm_dir_speed; added a call to
#               sub. storm_dir_speed to subs. read_jtwc_file, read_best_track, and
#               added best storm direction & speed parameters to subs. read_best_track,
#               combine_data and the main program; moved the 'use Date::Calc qw(Add_Delta_DHMS)'
#               pragma from sub. read_jtwc_file to the main program; added calls to
#               Date::Calc::Mktime to subs. read_jtwc_file, read_best_track; added
#               the date number arrays to sub. combine_data.  (RSL)
#  07 May 2012  Added a line to sub. read_best_track to skip best track records
#               with too few parameters (typical for the early stages of a storm).  (RSL)
#  08 May 2012  Upgraded sub. read_best_track to use the Delft3D missing value code
#               (1.0e30) instead of zero for missing data or initialization; added
#               the array @r100_65_50_35_3D to store all R34/R50/R64 values for
#               Delft3D-format output; upgraded sub. read_jtwc_file to load this
#               new array; upgraded sub. combine_data to include support for this
#               new array; added $envpress to the input parameter list for sub.
#               print_track_file; added an ELSE block for checking the value of
#               $want_best to load the @data_type array if best track handling
#               is not chosen and cyclone.dat output is chosen.  (RSL)
#  14 May 2012  Fixed an index bug in sub. read_jtwc_file pertaining to loading
#               the array @r100_65_50_35_3D.  (RSL)
#  18 May 2012  Modified the call to sub. get_year_month to obtain the storm ID
#               string in sub. read_jtwc_file; upgraded sub. get_year_month to create
#               & return a storm ID string; added creation of a TC warning ID string
#               from the input file name in sub. read_jtwc_file.  (RSL)
#  22 May 2012  Added a check to sub. storm_dir_speed for identical latitude,
#               longitude values.  (RSL)
#  25 May 2012  In sub. read_jtwc_file: added def. of $RMIN (minimum RMW value),
#               and added 2 IF blocks to check the calculated value of RMW and
#               reset any that are less than the min. to the min.; added def. of
#               $CP_MIN (minimum central pressure/MSLP), and added an IF block to
#               check the value of MSLP and reset any that are < the min. to the
#               min.  (RSL)
#  20 Aug 2012  Upgraded sub. write_tc_parameters: added writing of the storm name
#               to a file named 'cyclone_name.txt' for use by an external plotting
#               script.  (RSL)
#  13 Sep 2012  Upgraded sub. read_jtwc_file: commented-out the call to sub. interp_rads
#               and the last FOR loop inside the last IF block for Delft3D output, so
#               that interpolation of the R34/... wind radii is not performed.  (RSL)
#  19 Sep 2012  Upgraded sub. read_jtwc_file: changed the IF block within the main FOR
#               loop, so that the output @rmw is _always_ calculated based on the
#               ratio of max. wind speed to max. wind speed of R34/R50/R64; changed
#               debug 'print' statements to 'warn' statements & enclosed them in IF
#               blocks to print only if $DEBUG = 1; changed the IF block in which
#               the value of $rmw[$nrec] is compared to the value of $avg_dist[0]
#               so that if it is greater than or equal to that value, it is set to
#               $avg_dist[0] - 5 (nmi), per Sampson et al. (2010); modified the
#               code for the RMW value (@rmw) to use the wind speed ratio squared
#               (see comments below in sub. read_jtwc_file).  (RSL)
#  27 Sep 2012  Fixed an issue in sub. read_jtwc_file pertaining to when the $nrec
#               counter/index was incremented (perhaps there was a slight mod.
#               to the JTWC warning messages).  (RSL)
#  28 Sep 2012  Fixed an issue in sub. read_jtwc_file pertaining to when the $nrec
#               counter/index was incremented (modified the surrounding IF block to
#               check a different string); added an IF block to sub. read_jtwc_file
#               to skip the processing in the main FOR loop if a certain data
#               record is read, and added an IF block to decrement $nrec if it is
#               incorrectly incremented.  (RSL)
#  15 Feb 2013  Added new parameters, several IF blocks & several loops, as well as
#               a call to sub. write_ext_best_track (in best_track_routines.pl), to
#               add support for output in NRL extended best track format; updated
#               the help content.  (RSL)
#  22 Feb 2013  Added new parameters (@pres_outer_iso, @rad_outer_iso, @eye_diam,
#               @pres_outer_iso_w, @rad_outer_iso_w, @eye_diam_w) to the main
#               program for extended best track output; added an IF block to the
#               main program to call sub. read_std_best_track if extended best
#               track output format is requested; moved init. of @pres_outer_iso_w,
#               @rad_outer_iso_w, @eye_diam_w to just after the call to sub.
#               read_jtwc_file; upgraded sub. combine_data to add support for the
#               newly added arrays (see above for todays's mods).  (RSL)
#  08 Mar 2013  In the main program, moved the FOR loop that loads the %time_type
#               hash with the best track times to before the FOR loop that loads
#               the same hash with the TC warning file times; upgraded sub.
#               combine_data to load output arrays using best track records that
#               precede the TC warning 1st record (i.e., analysis), vice all best
#               track records.  (RSL)
#  12 Mar 2013  Added IF blocks within the best track code of the main program to
#               check the best track file type according to the file extension,
#               to call new sub. read_ext_best_track_file (in best_track_routines.pl),
#               and to provide support for NRL extended best track files.  (RSL)
#  15 May 2013  Added $READ_ALL & set it to 0 to indicate that only complete input best
#               track records are to be read & processed, & added it to the input
#               parameter list of the call to sub. read_std_best_track, with a check
#               for extended best track output (then the flag is set to 1); fixed
#               a typo in the IF block about the call to sub. print_wes_inp_file;
#               fixed an issue in sub. read_jtwc_file to correctly handle corrected
#               or relocated track files.  (RSL)
#
#*******************************************************************************

#  Enable strict syntax checking, and load the needed external routines.
use strict;
use File::Basename;
my $BaseDir = dirname($0);                  # Get name of directory in which this program resides
require "$BaseDir/cyc_routines.pl";  # External best track I/O routines
require "$BaseDir/interpl.pl";              # External linear interpolation routine
use Math::Trig qw(great_circle_distance deg2rad rad2deg great_circle_bearing);  # GC Distance, bearing

#  Use the following module to calculate date number (sec. since 1970-01-01 00:00:00).
use Date::Calc qw(Add_Delta_DHMS Mktime);

#  Declare the input, output file names & other parameters.
my ($jtwc_file, $trkfile, $best_track_file, $num_rad_bins, $num_dir_bins, $tc_radius);

#  Define a flag for whether Delft3D WES.EXE files will be generated.
my $want_best = 0;   # Default: 0 means no best track file will be used.

#  Declare (global) output parameters: date (YYYYMMDD), time (HHMM),
#  longitude (decimal deg.), latitude (dec. deg.), minimum sea level pressure
#  (MSLP) (hPa), radius of maximum winds (km), and maximum sustained winds (m/s),
#  analysis or forecast (flag value; 1 for analysis or 0 for forecast), TC direction,
#  TC speed (kt) and radii for R100/R65/R50/R35.
our ($tc_name, $tc_id, @date, @time, @lon, @lat, @mslp, @rmw, @msw_kt, @an_fcst, @tc_dir, @tc_speed, @r100_65_50_35, @datenum, $warn_id);
our $DEBUG = 0;  #  Debug print flag; set to zero to disable print/warn statements.

#  Similar parameters from the best track file:
our (@latt,@long,@RMW,@MSLP,@datestr,@max_wind_kt,@storm_name,@r100_65_50_35_3D,%time_type,@yr,@mon);
our (@R34,@R34_Rad,@R50,@R50_Rad,@R64,@R64_Rad,@data_type,@bt_dir,@bt_speed,@bt_datenum,@day,@hour);
our (@pres_outer_iso,@rad_outer_iso,@eye_diam,@pres_outer_iso_w,@rad_outer_iso_w,@eye_diam_w);
our $envpress = 1013.25;    # Environmental pressure (i.e., ambient press outside storm)
my $MISS = 1.0e30;          # Missing value code, per Deltares WES documentation (WES_UM_13.doc, p. 11).
my $EBT_MISS = -99;         # Missing value code for optional extended best track output
my $READ_ALL = 0;           # Flag to read only complete input records
#my @rad_outer_iso = ();     # Radius of the outermost closed isobar (ext. best track option)
#my @pres_outer_iso = ();    # Pressure of the outermost closed isobar (ext. best track option)
#my @eye_diam = ();         # Array of eye diameter values (ext. best track option)
my (@yyyy, @mm, @dd, @hr);  # Arrays of year, month, day, hour values (ext. best track option)
my ($sum34, $sum50, $sum64, @dist, @rad3d, $i);
my $FMT = "%d %.1f %.1f %.1f %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %.1f  %s\n";
my $usage = " Usage: $0 InFile num_rad_bins num_dir_bins tc_radius [best_track_file]\n for Extended Best Track or Delft3D output.\n";


#  Parse the command-line arguments.
if ( ($#ARGV+1) != 4 && ($#ARGV+1) != 5) {
    die $usage;

} else {

    #  Delft3D output:
    $jtwc_file = $ARGV[0];
    #  The optional inputs were entered, so parse them out.
    $num_rad_bins = $ARGV[1];  # Number of radial bins
    $num_dir_bins = $ARGV[2];  # Number of directional bins
    $tc_radius = $ARGV[3];     # TC radius (km)
    if (($#ARGV+1) == 5) {
        $best_track_file = $ARGV[4];  # Input best track file name
        $want_best = 1;
    }
}

#  Read the JTWC TC warning data file.
read_jtwc_file($jtwc_file);

#  Remove parentheses from TC name string.
$tc_name =~ s/[\(\)]//g;
#  Get the storm "human" name (e.g., "ALICE")
my @s_name = split ' ', $tc_name;
my $sname = uc($s_name[$#s_name]); #Change name to upper case

# Create file name to store cyclone info
my $DIR = dirname(${jtwc_file}); # Get the directory where the forecast file resides
# Get the 2 digit basin, 2 digit storm # and 2 digit yr from the JTWC warning filename
my $fname = fileparse($jtwc_file);
my $sid = substr($fname,0,6);
# Get the forecast date/time from the filename
my @tmp = split /[_.]/, $fname;
$trkfile="${DIR}/${sid}_${sname}_$tmp[1].cyc\n";
# New naming convention
#my $tmp1 = substr($fname,4,10);
#$trkfile="${DIR}/${sid}_${sname}_${tmp1}.cyc\n";

#  Load output arrays that have no values (or a default value) from the input file.
@rad_outer_iso_w = ($EBT_MISS) x scalar(@lat);
@eye_diam_w = ($EBT_MISS) x scalar(@lat);
@pres_outer_iso_w = ($envpress) x scalar(@lat);

#  DEBUG print of 3D radii;;;
#print " %%%%%%%%%%%%%%%%%% Before  read_best_track %%%%%%%%%%%%%%%%%%%\n";
#for my $i (0..$#r100_65_50_35_3D){
#  print "\nI = $i  ";
#  for my $j (0..$#{$r100_65_50_35_3D[$i]} ) {
#    print "J = $j 100_65_50_35_3D = ";
#    for my $k (0..3) {
#      print " $r100_65_50_35_3D[$i][$j][$k]";
#    }
#    print "\n";
#  }
#}
#print " %%%%%%%%%%%%%%%%%% Before  read_best_track %%%%%%%%%%%%%%%%%%%\n";


#  Check whether a best track file was supplied.
if ($want_best == 1) {
    #warn "BEST Track support...File is $best_track_file\n";

    #  Check whether output is extended best track format; if so, then read all input best track records.
    $READ_ALL = 1;
    ($envpress,*latt,*long,*RMW,*MSLP,*datestr,*max_wind_kt,*storm_name,*R34,*R34_Rad,*R50,*R50_Rad,
        *R64,*R64_Rad,*bt_datenum,*pres_outer_iso,*rad_outer_iso,*eye_diam) = read_std_best_track($READ_ALL, $best_track_file);

    #  Loop over best track data records....
    for $i (0..$#R34_Rad) {
        #  Load the hash with the current date/time group as the key and a string description & index as the values.
        $time_type{$datestr[$i]} = ['B', $i];

        #  Debug printing:
        # printf $FMT, $datestr[$i], $latt[$i], $long[$i], $max_wind_kt[$i], $MSLP[$i], $RMW[$i], $R34[$i], $R34_Rad[$i][0], $R34_Rad[$i][1],  $R34_Rad[$i][2], $R34_Rad[$i][3], $R50[$i], $R50_Rad[$i][0], $R50_Rad[$i][1], $R50_Rad[$i][2], $R50_Rad[$i][3], $R64[$i], $R64_Rad[$i][0], $R64_Rad[$i][1], $R64_Rad[$i][2], $R64_Rad[$i][3], $storm_name[$i];
        # printf join("\n". @storm_name);
    }

    #  Load the hash with the dates & description + index of the "normal" warning file.
    #  This is loaded 2nd so that any duplicate date/time values from the best
    #  track file will be overwritten by those from the warning file (since they
    #  are forecasts, they should be more recent).
    #print "] ] ] ] ] ] Loading the time_type hash with warning data....\n";
    for $i (0..$#date) {
        my $dtg = sprintf('%d%02d', $date[$i], substr($time[$i],0,2));
        $time_type{$dtg} = ['W', $i];
        #print "I = $i  DTG = $dtg\n";
    }

    #  Print the data in a tabular format.

    #print "Back in the MAIN program; BG or Env. Pressure = $envpress; ", scalar @latt, " Best Track data records:\n";

    #print "Date/time  Lat  Lon  VMax  MSLP  RMW  R34  R34_1  R34_2  R34_3  R34_4  R50  R50_1  R50_2  R50_3  R50_4  R64  R64_1  R64_2  R64_3  R64_4   Name\n";

    #  Print the hash data for debugging.
    #for my $key (sort keys %time_type) {
    #  print "Date: $key  Desc: $time_type{$key}[0]  Index:$time_type{$key}[1] \n";
    #}

    #  Combine the data sets in chronological order in the original warning file data arrays.
    combine_data();

}

#  Check the signs of the longitudes, and correct any discrepancies.
check_longitude();

#  DEBUG print of 3D radii;;;
#for my $i (0..$#r100_65_50_35_3D){
#  print "\nI = $i  ";
#  for my $j (0..$#{$r100_65_50_35_3D[$i]} ) {
#    print "J = $j 100_65_50_35_3D = ";
#    for my $k (0..3) {
#      print " $r100_65_50_35_3D[$i][$j][$k]";
#    }
#    print "\n";
#  }
#}


#  Write the pertinent TC values.

#print "Dir = @tc_dir\nSpd = @tc_speed\n";
#  Check the output option.

#  Create an array of the storm name.
@storm_name = ($s_name[$#s_name - 1]) x scalar(@lat);

#  Loop over data records....
for my $i (0..$#date) {
    #  Break the date/time group into year, month, day, hour.
    $yyyy[$i] = substr($date[$i], 0, 4) + 0;  # Year
    $mm[$i] = substr($date[$i], 4, 2) + 0;    # Month
    $dd[$i] = substr($date[$i], 6, 2) + 0;    # Day
    $hr[$i] = substr($time[$i], 0, 2) + 0;    # Hour

    #  Loop over wind radii to convert missing values to zero.
    for my $ii (0..$#{$r100_65_50_35_3D[$i]}) {
        #  Loop over directions....
        for my $jj (0..3) {
            #  Check whether this radius value is a missing value.
            if ($r100_65_50_35_3D[$i][$ii][$jj] == $MISS) {
                #  It is, so convert to zero.
                $r100_65_50_35_3D[$i][$ii][$jj] = 0;
            }
            #  Load the output 3D array with the current wind radii values.
            $rad3d[$i][$ii][$jj] = $r100_65_50_35_3D[$i][$ii][$jj];
        }
    }

    #  Determine which of the R34/50/64 wind radii are nonzero by summing each.
    $sum34 = $r100_65_50_35_3D[$i][3][0] + $r100_65_50_35_3D[$i][3][1] + $r100_65_50_35_3D[$i][3][2] + $r100_65_50_35_3D[$i][3][3];
    $sum50 = $r100_65_50_35_3D[$i][2][0] + $r100_65_50_35_3D[$i][2][1] + $r100_65_50_35_3D[$i][2][2] + $r100_65_50_35_3D[$i][2][3];
    $sum64 = $r100_65_50_35_3D[$i][1][0] + $r100_65_50_35_3D[$i][1][1] + $r100_65_50_35_3D[$i][1][2] + $r100_65_50_35_3D[$i][1][3];

    #  Define the indices of the R34/50/64 wind radii, depending on which values are nonzero.
    if ($sum64 == 0 && $sum50 == 0 && $sum34 >= 0) {
        @dist = (1);
        #  #  R34 values may be present, R50 & R64 are not.
    } elsif ($sum64 == 0 && $sum50 > 0 && $sum34 > 0) {
        @dist = (3,2);
        #  #  R34 & R50 values are present, R64 are not.
    } elsif ($sum64 > 0 && $sum50 > 0 && $sum34 > 0) {
        @dist = (3,2,1);
        #  All 3 radii are present.
    }
    #  If more than 1 set of wind radii exists, swap the 1st & last ones for the output routine.
    if (scalar(@dist) > 1) {
        #  Loop over directions....
        for my $j (0..3) {
            #  Set the last of the list to the 1st.
            $rad3d[$i][$dist[$#dist]][$j] = $r100_65_50_35_3D[$i][$dist[0]][$j];
            #  Set the 1st of the list to the last.
            $rad3d[$i][$dist[0]][$j] = $r100_65_50_35_3D[$i][$dist[$#dist]][$j];
        }
    }
    #  Clear out local vars.
    @dist = {}; $sum34 = 0; $sum50 = 0; $sum64 = 0;
    #warn "\n";
}



#  Output to extended best track format.
write_cyc($trkfile,$tc_name,\@eye_diam_w,\@pres_outer_iso_w,\@rad_outer_iso_w,\@yyyy,\@mm,\@dd,\@hr,
    \@lat,\@lon,\@msw_kt,\@rmw,\@rad3d,\@mslp,$tc_radius,$num_rad_bins,$num_dir_bins);

print " Cyclone parameters written to: ${trkfile}";

exit;

#*******************************************************************************
sub read_jtwc_file {
    #
    #  This subroutine reads a JTWC warning bulletin file and extracts tropical
    #  cyclone (TC) parameters needed by an external meteorological model.
    #
    #  Syntax: read_jtwc_file($infile);
    #  where:  $infile is the name of the input warning filewp9411web_201107242300.txt
    #
    #  Calls:
    #
    #  Called by: main program
    #
    #  References:  Holland, G.J., J.I. Belanger, and A. Fritz (2010). A revised model for radial
    #               profiles of hurricane winds, Mon. Wea. Rev, 128, 4393-4401.
    #               Sampson, C.R., P.A. Wittmann, and H.L. Tolman (2010). Consistent tropical
    #               cyclone wind and wave forecasts for the U.S. Navy, Wea. Forecasting, 25,
    #               1293-1306.
    #
    #  Revision History:
    #  18 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  19 Jul 2011  Wrapped up initial coding; began testing, debugging.  (RSL)
    #  20 Jul 2011  Added code to open, write, close a temporary text file with lat,
    #               lon, wind speed, several RMW values for testing of RMW algorithms.  (RSL)
    #  21 Jul 2011  Changed $becoming to $extra; modified the IF blocks checking for
    #               "BECOMING EXTRATROPICAL" to check only for "EXTRATROPICAL"; added
    #               code to check the RMW so that if the max. sustained wind value is the
    #               same as one of the wind radii, then that mean radius is used for
    #               the RMW -- e.g., if MSW = 50 kt, and RADIUS OF 050 KT WINDS is present,
    #               the RMW is set to the mean of the R50 values.  (RSL)
    #  22 Jul 2011  Modified the regexp when initially setting $long_range to search for
    #               LONG RANGE or EXTENDED.  (RSL)
    #  25 Jul 2011  Modified the temp file name ($temp) to simply prepend 'rmw_' to the input
    #               warning file name; added an IF block to check for a tropical cyclone
    #               formation alert, and to print a message if one is being processed, and
    #               then exit the program.  (RSL)
    #  31 Oct 2011  Added code to obtain TC speed & direction information for WES.EXE
    #               output, including importing Math::Trig, and the associated functions,
    #               and added several arrays & scalars to calculate, store the various
    #               intermediate results.  (RSL)
    #  01 Nov 2011  Fixed a bug in the determination of date for the forecasts: added
    #               use of Date::Calc::Add_Delta_DHMS, using the analysis date/time and
    #               the forecast tau values; added an IF block to the FOR loop in which
    #               @r100_65_50_35 is loaded, to check whether any of the R100 values
    #               are less than RMW, and if so, set to the missing value.  (RSL)
    #  16 Dec 2011  Added support to extract the storm name in $tc_name.  (RSL)
    #  26 Jan 2012  Commented-out debug code for the temp file.  (RSL)
    #  06 Feb 2012  Changed the call to die() to a call to warn() when an TC formation
    #               alert is being processed, as such an alert may still have TC product
    #               information.  (RSL)
    #  04 May 2012  Moved the 'use Math::Trig' and 'use Date::Calc' pragmas to the top
    #               of the main program; moved the storm translation speed & direction
    #               calculations to newly added sub. storm_dir_speed, and added a call
    #               to that sub; added calculation of the date number.  (RSL)
    #  08 May 2012  Added the $i3d, $k, $ind_sio, $nquad parameters, and added code
    #               to store wind radii parameters in the 3-D array @r100_65_50_35_3D
    #               for each of the 4 quadrants where available.  (RSL)
    #  14 May 2012  Fixed a bug in the index finding FOR loop for @r100_65_50_35_3D
    #               for the case of > 1 wind radii.  (RSL)
    #  18 May 2012  Added a call to fileparse() to strip off any leading subdir. from
    #               $infile; modified the call to sub. get_year_month to use the bare
    #               input file name & to obtain the storm ID string; added creation of
    #               a storm warning ID string based on the input file name.  (RSL)
    #  25 May 2012  Added def. of $RMIN (minimum RMW value), and added 2 IF blocks to
    #               check the calculated value of RMW and reset any that are less than
    #               the min. to the min.; added def. of $CP_MIN (minimum central pressure/
    #               MSLP), and added an IF block to check the value of MSLP and reset any
    #               that are < the min. to the min.  (RSL)
    #  13 Sep 2012  Commented-out the call to sub. interp_rads and the last FOR loop
    #               inside the last IF block for Delft3D output, so that interpolation
    #               of the R34/... wind radii is not performed.  (RSL)
    #  19 Sep 2012  Changed the IF block within the main FOR loop, so that the output @rmw is
    #               _always_ calculated based on the ratio of max. wind speed to max. of
    #               R34/R50/R64; changed debug 'print' statements to 'warn' statements &
    #               enclosed them in IF blocks to print only if $DEBUG = 1; changed the
    #                IF block in which the value of $rmw[$nrec] is compared to the value
    #               of $avg_dist[0] so that if it is greater than or equal to that value,
    #               it is set to $avg_dist[0] - 5 (nmi), per Sampson et al. (2010);
    #               modified the code for the RMW value (@rmw) to use the wind speed
    #               ratio squared (see comments below).  (RSL)
    #  27 Sep 2012  Added an IF block about the increment of $nrec to perform only
    #               for certain specific input records; reverted the IF block where the
    #               $long_range flag is set to use the original case (only 'LONG RANGE OUTLOOK').
    #               (RSL)
    #  28 Sep 2012  Modified the IF block in which $nrec is incremented, to check the
    #               contents of $line[0] instead of /EXTRATROPICAL$/, to work correctly
    #               on any system (was behaving unexpectedly on Windows); added an IF block
    #               to decrement $nrec if it was incorrectly incremented; added an IF block
    #               to skip the remainder of the FOR loop if a record with the text
    #               "BECOMING EXTRATROPICAL" is read.  (RSL)
    #  15 May 2013  Added def. of $iend to correctly handle parsing of the "WARNING NR"
    #               record when it includes "RELOCATED" or "CORRECTED" in that record.
    #               (RSL)
    #
    #*******************************************************************************

    #  Parse out the input parameter(s).
    my ($infile) = @_;

    #  Declare, initialize lodal parameters.
    my @line;
    my $is_analysis = 0;
    my ($nrec, $irad, $idist, $long_range, $extra, $iend) = (0) x 6;
    my ($yr0,$mon0,$day0,$hr0,$min0,$sec0,$Dd,$Dh,$Dm,$Ds,$yr_now,$mon_now,$day_now,$hr_now,$min_now,$sec_now) = (0) x 16;
    my ($temp, $rmw_orig, $rmw0, $rmw1, $rmw2, $i, $j, $ilat, $ilon, $sum, $num, @avg_dist, @rad_speed, @dist, @sorted_dist, @valu, $yyyymm);
    #  TC direction, speed parameters:
    my ($inm, @start, @end, @tau, $dir, $km, @distnc, @reverse_s, @reverse_d, @r100, $ir100, $i3d, $k, $ind_sio, $previous);
    my $nm2km = 1.852;      # Conversion for nmi to km
    my @dash = ('-') x 80;
    my $maxlen = -1;
    my $nquad = 4;  # Number of quadrants for R35/R50/...
    my $RMIN = 6/$nm2km;  # Minimum RMW in NM (6 km) per Deltares WES documentation (WES_Sci_tech_v8.doc, p. 7).
    my $CP_MIN = 870;     # Lowest known CP of Pacific Basin (T. Tip, 1979), mb
    #  Wind speed values needed for the corresponding WES.EXE radii, in KT:
    my @spd_needed = (35,50,65,100);
    my @spd_in_order = (100,64,50,34);
    my $MISS = 1.0e30;        # Missing value code, per Deltares WES documentation (WES_UM_13.doc, p. 11).
    my $RMW_MAX = 500/$nm2km; # Maximum radius of max. winds (typ. TS-force wind radius is 483 km), km converted to nmi

    #  Temporary output file for development, debugging (RSL, Jul 2011)
    #$temp = 'rmw_' . $infile;
    #print "TEMP file is $temp\n";
    #open(MYTEMP, ">$temp")  or die "ERROR: Unable to open output file $temp ($!); stopped";

    #  Extract only the file name without a subdir. name.
    my ($fname, $p, $e) = fileparse($infile);

    #  Extract a warning ID from the file name, without the extension.
    ($warn_id, $p, $e) = fileparse($infile, ('.txt',".txt_[0-9]{14}"));

    #  Determine the year and month of this warning file from the file name.
    ($yyyymm, $tc_id) = get_year_month($fname);

    if ($DEBUG) { warn "\n  NOTE:  TC ID string is $tc_id; Warning ID is $warn_id\n\n"; }

    open(INPUT, "<$infile") or die "ERROR: Unable to open input file $infile ($!); stopped";

    #  Loop over input file records....
    while(<INPUT>) {

        #  Check for a TC formation alert.
        if ($_ =~ /TROPICAL CYCLONE FORMATION ALERT/) {
            warn "\n NOTE: The present file, $infile, is a Tropical Cyclone Formation Alert.\n Formation of a significant tropical cyclone is possible within the next 12 to 24 hours.\n\n";
        }

        #  If the record contains the following string, drop to the end of the loop (i.e., skip the rest).
        if ($_ =~ /BECOMING EXTRATROPICAL/) { next; }

        #  Break this record into fields.
        @line = split ' ', $_;

        #  Break out of the loop after the data records.
        if ($line[0] =~ /REMARKS:/) {
            #    if ($line[0] =~ /REMARKS:/ || ($long_range == 1 && $line[0] == "EXTRATROPICAL")) {
            if ($DEBUG) { warn "** NREC = $nrec; Date has ", $#date+1," elems; RMW has ", $#rmw+1," elems.\n"; }
            if (scalar @rmw != scalar @date) {
                $rmw[$nrec] = 35.37 - (0.111 * $msw_kt[$nrec]) + (0.570 * (abs($lat[$nrec]) - 25));  # Gross et al. (2004) in nmi

                #  Check whether this RMW is too small, and reset to the minimum if so.
                if ($rmw[$nrec] < $RMIN) {
                    warn " NOTE: RMW[$nrec] ($rmw[$nrec]) < min. ($RMIN); resetting to min.\n";
                    $rmw[$nrec] = $RMIN;
                }
                if ($DEBUG) { warn " ** Now RMW has", $#rmw+1," elems.\n"; }
            }
            last;
        }

        #  Check for the storm name, with a line that begins with "1." and has the indicated text.
        if ($_ =~ /^1\./ && $_ =~ / WARNING NR /) {
            #  Define the ending index according to whether the last field is a number.
            if ($line[$#line] =~ /^\d{3}$/) {
                #  Last field is a number, so use that.
                $iend = $#line - 3;
            } else {
                #  Last field is not a number, so want to subtract 4.
                $iend = $#line - 4;
            }
            #  Loop over a subset of the fields.
            for $inm (1..$iend) {
                #  Check whether this field contains numerals.
                if ($line[$inm] !~ /\d+/) {
                    #  No numerals, so process the text.
                    #  Check for parentheses.
                    if ($line[$inm] =~ /^\(/) {
                        #  This line has parens, so extract only the text.
                        my $tmp = substr($line[$inm],1,length($line[$inm])-2);
                        #  Append the parens with initialcaps text inside.
                        $tc_name .= '(' . ucfirst(lc($tmp)) . ')';
                    } else {
                        #  Append the initialcaps text field.
                        $tc_name .= ucfirst(lc($line[$inm]));
                    }
                } else {
                    #  A storm number; append it as-is.
                    $tc_name .= $line[$inm];
                }
                #  Append a space if not the last field.
                if ($inm < ($#line-3)) { $tc_name .= ' '; }
            }
        }

        #  Initialize the R100 array with missing value codes.
        for $ir100 ( 0..$#spd_needed ) {
            $r100_65_50_35[$nrec][$ir100] = $MISS;
            for ($i3d = 0; $i3d < $nquad; $i3d++) {
                $r100_65_50_35_3D[$nrec][$ir100][$i3d] = $MISS;
            }
        }

        #  If the long range outlook has been reached, set a flag.
        if ($_ =~ /LONG RANGE OUTLOOK:/) { $long_range = 1; }
        ##if ($_ =~ /(LONG RANGE|EXTENDED) OUTLOOK:/) { $long_range = 1; }
        #if ($_ =~ /(LONG RANGE|EXTENDED) OUTLOOK:/ || $_ =~ /DISSIPATED/ ) { $long_range = 1; }

        #  Find analysis values.
        if ($line[0] =~ /WARNING/ && $line[1] =~ /POSITION:/) {
            $is_analysis = 1;
            $an_fcst[$nrec] = $is_analysis;
        } elsif ($line[0] =~ /FORECASTS:/) {
            $is_analysis = 0;
            $an_fcst[$nrec] = $is_analysis;
        }

        #  If this is a forecast, get the forecast time (i.e., tau) (hr)
        if ($_ =~ /VALID AT:/) {
            $tau[$nrec] = $line[0];
        }

        #  Get date, time, position:
        if ($line[0] =~ /^(\d){6}Z/) {
            if ($DEBUG) { warn "DATE: $_"; }

            #  Time & position of analysis.
            if ($is_analysis) {
                #  Use the analysis date/time.
                $date[$nrec] = sprintf("%d%02d", $yyyymm, substr($line[0],0,2));
                $time[$nrec] = sprintf("%02d%02d", substr($line[0], 2, 2), substr($line[0], 4, 2));
                $yr0 = substr($yyyymm, 0, 4);    # Starting year
                $mon0 = substr($yyyymm, 4, 2);   # Starting month
                $day0 = substr($line[0],0,2);    # Starting day
                $hr0 = substr($line[0], 2, 2);   # Starting hour
                $min0 = substr($line[0], 4, 2);  # Starting minute
                $sec0 = 0;                       # Starting second
            } else {
                #  If this is a forecast (i.e., past the analysis data set), then calculate the date
                #  based on the analysis date plus the tau value.
                ($yr_now,$mon_now,$day_now,$hr_now,$min_now,$sec_now) = Add_Delta_DHMS($yr0,$mon0,$day0,$hr0,$min0,$sec0,
                    $Dd,$tau[$nrec],$Dm,$Ds);
                $date[$nrec] = sprintf("%4d%02d%02d", $yr_now, $mon_now, $day_now);
                #$time[$nrec] = sprintf("%02d%02d", $hr_now, $min_now, $sec_now);
                $time[$nrec] = sprintf("%02d%02d", $hr_now, $min_now);
            }

            #  Get the latitude & longitude, respectively.
            #  Indices of latitude & longitude.
            $ilat = $#line - 1;
            $ilon = $#line;

            #  Latitude
            $lat[$nrec] = substr($line[$ilat], 0, -1) + 0;  # Force to be a number
            #  Reverse the latitude sign if in the southern hemisphere
            if (substr($line[$ilat], -1) =~ /S/) { $lat[$nrec] *= -1; }
            #  Longitude
            $lon[$nrec] = substr($line[$ilon], 0, -1) + 0;  # Force to be a number
            #  Reverse the longitude sign if in the western hemisphere
            if (substr($line[$ilon], -1) =~ /W/) { $lon[$nrec] *= -1; }

            #  If this is past the analysis data, then calculate the TC speed & direction.
            if ($nrec > 0) {
                ($tc_dir[$nrec], $tc_speed[$nrec]) = storm_dir_speed($lon[$nrec-1], $lat[$nrec-1],
                    $lon[$nrec], $lat[$nrec], $tau[$nrec-1], $tau[$nrec]);


                if ($DEBUG) { warn " ^^^^^^vvv TC_Dir = $tc_dir[$nrec]; TC_Speed = $tc_speed[$nrec]; Delta_Time = ", $tau[$nrec] - $tau[$nrec-1] ," hr.\n"; }
            }
        }

        #  Get initial movement values (direction, speed).
        if ($_ =~ /MOVEMENT PAST SIX HOURS/) {
            #  Storm direction, degrees:
            $tc_dir[$nrec] = $line[5];
            #  Storm translation speed, knots:
            $tc_speed[$nrec] = $line[8];
            #  Store the time (should be zero here):
            $tau[$nrec] = 0;

            if ($DEBUG) { warn " ^^^^^^ TC_Dir = $tc_dir[$nrec]; TC_Speed = $tc_speed[$nrec].\n"; }
        }

        #  Get max. sustained winds, radius of maximum winds, min. sea level pressure.
        if ($_ =~ /MAX SUSTAINED WINDS -/) {
            #  Maximum sustained winds (in knots):
            $msw_kt[$nrec] = $line[4] + 0;
            #  Calculate the minimum sea level pressure from the wind speed.
            $mslp[$nrec] = wind_pressure_relationship($msw_kt[$nrec]);

            #  Check whether the MSLP/central pressure is below the known minimum,
            #  reset to the min. if so.
            if ($mslp[$nrec] < $CP_MIN) {
                warn " NOTE: MSLP($nrec) ($mslp[$nrec]) < min. ($CP_MIN); resetting to the min.\n";
                $mslp[$nrec] = $CP_MIN;
            }
            if ($DEBUG) { warn "  --  MSW = $msw_kt[$nrec], MSLP = $mslp[$nrec].\n"; }

            #  Now increment the index.
            #$nrec++;
        }

        #  Obtain wind radius information.
        if ($_ =~ /RADIUS OF/) {
            $rad_speed[$irad] = $line[2] + 0;
            #  Unset the long_range flag if this record has wind radii.
            ###$long_range = 0;
        }
        if ($_ =~ /QUADRANT/) {
            if ($_ =~ /NORTHEAST/) {
                $idist = 0;
            } else {
                $idist++;
            }
            if ($#line > 3) {
                $dist[$irad][$idist] = $line[6];
            } else {
                $dist[$irad][$idist] = $line[0];
            }
            if ($_ =~ /NORTHWEST QUADRANT/) {
                $irad++;
            }
        }
        #  Track the number of wind radii.
        if (scalar @rad_speed > $maxlen) {
            $maxlen = scalar @rad_speed;
        }
        #print "  -- -- MaxLen = $maxlen.\n";

        if ($_ =~ /REPEAT POSIT:/ ||  $_ =~ /VECTOR TO/ ||
            #($long_range == 1 && ($_ =~ /BECOMING/ ||  $_ =~ /EXTRATROPICAL/))) {
            ($long_range == 1 && $_ =~ /EXTRATROPICAL/)) {

            #  Set the extratropical flag if found.
            #if ($long_range == 1 && ($_ =~ /BECOMING/ ||  $_ =~ /EXTRATROPICAL/)) { $becoming = 1; }
            if ($long_range == 1 && $_ =~ /EXTRATROPICAL/) { $extra = 1; }
            if ($DEBUG) { warn "\n\nNOTE: INPUT RECORD IS:\n$_\n"; }
            if ($DEBUG) { warn "LONG_RANGE = $long_range, EXTRA = $extra.\n"; }

            #if ($long_range == 1 && $extra == 1 && $_ =~ /VECTOR TO/) {
            #if ($long_range == 1 && $extra == 1 && $_ =~ /VECTOR TO/ && scalar(@rad_speed) == 0) {
            #  We want to skip incrementing everything.
            #  if ($DEBUG) { warn " * * * * * * * * * * Skipping the increment.... * * * * * * * * * *\n"; }
            #if ($DEBUG) { warn ")))NREC = $nrec, RAD_SPEED has ", $#rad_speed+1, " elems; MSW_KT has ", scalar(@msw_kt), "elems; DIST is ",$#dist+1," x ", $#{$dist[0]}+1," in size.\n"; }
            #next;
            #}


            if ($DEBUG) { warn "NREC = $nrec, RAD_SPEED has ", $#rad_speed+1, " elems; MSW_KT has ", scalar(@msw_kt), "elems; DIST is ",$#dist+1," x ", $#{$dist[0]}+1," in size.\n"; }
            #  Determine the radius of max winds.  Calculate from a formula if no radii were available.
            if (($#rad_speed+1) < 1) {
                #  No wind radii are present.
                if (scalar(@msw_kt) == $nrec) {
                    if ($DEBUG) { warn " * * * * * * * * * * NREC appears to be too large ($nrec); decrementing it.... * * * * * * * * * *\n"; }
                    $nrec--;
                }
                #  Calculate the radius of maximum winds (in km) using a modified Rankine vortex
                #  formulation (Knaff and Zehr, 2007).  This was stated to be developed using Atlantic storms.
                #$rmw[$nrec] = 66.785 - (0.09102 * $msw_kt[$nrec]) + (1.0619 * (abs($lat[$nrec]) - 25));
                #  Formulation of Gross et al. (2004) (in nmi).  This was stated to be developed using Atlantic storms.
                #  This is for the case of having no wind radii, as a fail-safe.
                if ($DEBUG) { warn "RMW by Gross et al. (2004) formula....\n"; }
                $rmw[$nrec] = 35.37 - (0.111 * $msw_kt[$nrec]) + (0.570 * (abs($lat[$nrec]) - 25));
                $rmw0 = 0; $rmw1 = 0; $rmw2 = 0; $rmw_orig = 0;
            } else {
                #  At least 1 set of wind radii are present.

                #  Use both RMW algorithms to check whether either is better than the above.
                $rmw0 = (66.785 - (0.09102 * $msw_kt[$nrec]) + (1.0619 * (abs($lat[$nrec]) - 25)))/$nm2km;  # Knaff & Zehr (2007) in km, converted to nmi
                $rmw1 = 35.37 - (0.111 * $msw_kt[$nrec]) + (0.570 * (abs($lat[$nrec]) - 25));  # Gross et al. (2004) in nmi

                #  Initialize the averaging parameters for this record.
                $sum = 0;
                $num = 0;
                @avg_dist = ();

                #  Check the number of wind radii.
                if (($#rad_speed+1) == 1) {
                    #  Loop over spd_in_order to find the matching index for this wind speed.
                    for $k (0..$#spd_in_order) {
                        #  Check whether this wind speed matches that read from the file.
                        if ($spd_in_order[$k] == $rad_speed[0]) {
                            #  Matching wind speed found, so store it & drop out of loop.
                            $ind_sio = $k;
                            last;
                        }
                    }

                    #  Only 1 set of radii, so use a ratio
                    for $i (0..$#{$dist[0]}) {
                        #  Store this radius in the 3-D array.
                        $r100_65_50_35_3D[$nrec][$ind_sio][$i] = $dist[0][$i] + 0;
                        $sum += $dist[0][$i];
                        $num++;
                    }
                    $avg_dist[0] = $sum/$num;
                    if ($DEBUG) { warn "AVG_Dist1 = $avg_dist[0]; RAD_Speed = $rad_speed[0].\n"; }
                    #$rmw[$nrec] = ($msw_kt[$nrec]/$rad_speed[0]) *  $avg_dist[0];  # Too large;
                    # Ratio of wind radius max. speed to MSW approach
                    #$rmw[$nrec] = ($rad_speed[0]/$msw_kt[$nrec]) *  $avg_dist[0];

                    #  Always use the ratio of the max. forecast wind speed to
                    #  the wind speed w/radii to estimate the radius of max. winds.
                    #  Note that the ratio is squared; cf. Holland et al. (2010), p. 4394.
                    $rmw[$nrec] = (($rad_speed[0]/$msw_kt[$nrec])**2) *  $avg_dist[0];
                    $rmw2 = $rmw[$nrec];

                } else {
                    #  More than 1 set of wind radii are present.
                    #  Extrapolate.
                    for $i (0..$#dist) {
                        #  Loop over spd_in_order to find the matching index for this wind speed.
                        for $k (0..$#spd_in_order) {
                            #  Check whether this wind speed matches that read from the file.
                            if ($spd_in_order[$k] == $rad_speed[$i]) {
                                #  Matching wind speed found, so store it & drop out of loop.
                                $ind_sio = $k;
                                last;
                            }
                        }
                        for $j (0..$#{$dist[$i]}) {
                            #  Store this radius in the 3-D array.
                            $r100_65_50_35_3D[$nrec][$ind_sio][$j] = $dist[$i][$j] + 0;
                            $sum += $dist[$i][$j];
                            $num++;
                        }
                        $avg_dist[$i] = $sum/$num;
                        if ($DEBUG) { warn "AVG_Dist $i = $avg_dist[$i], RAD_Speed = $rad_speed[$i].\n"; }
                        $sum = 0;
                        $num = 0;
                    }  #  end for $i...

                    # Original: $rmw[$nrec] = abs(interp($rad_speed[0],$rad_speed[$#rad_speed],$msw_kt[$nrec],$avg_dist[0],$avg_dist[$#avg_dist]));
                    #$rmw[$nrec] = abs(interp($rad_speed[1],$rad_speed[0],$msw_kt[$nrec],$avg_dist[1],$avg_dist[0]));

                    #  Always use the ratio of the max. forecast wind speed to
                    #  the wind speed w/radii to estimate the radius of max. winds.
                    #  Note that the ratio is squared; cf. Holland et al. (2010), p. 4394.
                    $rmw[$nrec] = (($rad_speed[0]/$msw_kt[$nrec])**2) *  $avg_dist[0];  # Ratio of max. wind radius speed to MSW approach
                    $rmw2 = ($rad_speed[0]/$msw_kt[$nrec]) *  $avg_dist[0];  # Ratio of max. wind radius speed to MSW approach

                }  # end if(($#rad_speed...
                #  Store the calculated RMW.
                $rmw_orig = $rmw[$nrec];
                #  Sort the averaged wind radii, to ensure that the 1st element is the minimum.
                @sorted_dist = sort {$a <=> $b} @avg_dist;
                #  Find a reasonable value for RMW based on calculated values

                #print "   Just before calling find_best_value: RMW2,0,1,nrec = $rmw2, $rmw0, $rmw1, $rmw[$nrec].\n";

                #$rmw[$nrec] = find_best_value($sorted_dist[0], $rmw2, $rmw0, $rmw1, $rmw[$nrec]);

                #print "   Just after calling find_best_value: new RMW = $rmw[$nrec].\n";

                #  Sanity check for whether the MSW is the same as the innermost wind radius;
                #  e.g., if MSW=50 KT and R50 is present, then the RMW is set to the mean of the R50 values.
                #if ($msw_kt[$nrec] == $rad_speed[0]) {
                #  print "    Sanity check -- RMW changed from $rmw[$nrec] to $avg_dist[0].\n";
                #  $rmw[$nrec] = $avg_dist[0];
                #}

                #  Interpolate from NHC/JTWC wind radii to WES.EXE radii.
                @reverse_s = reverse @rad_speed;
                @reverse_d = reverse @avg_dist;
                #@r100 = interp_rads(\@reverse_s, \@spd_needed, \@reverse_d );
                if ($DEBUG) { warn " -- After call to interp_rads: RadSpeed = @reverse_s; SpdNeeded = @spd_needed, AVG_Dist = @reverse_d , R100 = @r100.\n"; }

                #  Fill the R100 array with the interpolated  values.
                #for $ir100 ( 0..$#spd_needed ) {
                #  #  Check whether the R100 value exceeds the RMW; set to missing value if it does not.
                #  if ($r100[$ir100] <= $rmw[$nrec]) {
                #    #  This R100/R65/R50/R35 value is less than RMW, so set it to the missing value.
                #    $r100_65_50_35[$nrec][$ir100] = $MISS;
                #  } else {
                #    #  It is valid, so store it.
                #    $r100_65_50_35[$nrec][$ir100] = $r100[$ir100];
                #  }
                #}
                $ind_sio = 0;

                ###        if ($rmw0 > $sorted_dist[0] || $rmw1 > $sorted_dist[0]

            }  #  end if (($#rad_speed+1) < 1)

            #  Gross error check:
            #  Check whether this RMW is too small, and reset to the minimum if so.
            if ($rmw[$nrec] < $RMIN) {
                if ($DEBUG) { warn " NOTE: Calculated RMW[$nrec] ($rmw[$nrec]) < min. ($RMIN); resetting to min.\n"; }
                $rmw[$nrec] = $RMIN;
                #} elsif ($rmw[$nrec] >= $avg_dist[0]) {
            } elsif (scalar(@avg_dist) > 0 && $rmw[$nrec] >= $avg_dist[0]) {
                #  If calculated RMAX is greater than the (average) minimum of the R34/R50/R64 radii, then it is set
                #  to the minimum - 5 nmi (Sampson et al., 2010, p. 1296).
                if ($DEBUG) { warn "Calculated RMW ($rmw[$nrec]) is greater than the (average) minimum wind radius for this data set ($avg_dist[0]); setting to the minimum - 5 nmi.\n"; }
                $rmw[$nrec] = $avg_dist[0] - 5;
            } elsif ($rmw[$nrec] > $RMW_MAX) {
                if ($DEBUG) { warn "  WARNING: RMW > Maximum \($RMW_MAX km\); resetting to Maximum.\n"; }
                $rmw[$nrec] = $RMW_MAX;
            }

            if ($DEBUG) { warn ":::Date, time, lat, lon, MSW, RMW, MSLP, RMW0,1,2,Orig = $date[$nrec], $time[$nrec], $lat[$nrec], $lon[$nrec], $msw_kt[$nrec], $rmw[$nrec], $mslp[$nrec], $rmw0, $rmw1, $rmw2, $rmw_orig.\n"; }
            #print substr($date[$nrec],0,4),substr($date[$nrec],4,2),substr($date[$nrec],6,2),substr($time[$nrec],0,2),substr($time[$nrec],2,2),substr($time[$nrec],4,2),"\n";
            $datenum[$nrec] = Mktime(substr($date[$nrec],0,4)+0,substr($date[$nrec],4,2)+0,substr($date[$nrec],6,2)+0,
                substr($time[$nrec],0,2)+0,substr($time[$nrec],2,2)+0, 0);
            #print "  MaxLen = $maxlen.\n";
            if ($DEBUG) { warn @dash, "\n"; }

            #  Resize the value array, then fill it with the current values of avg_dist.
            #  The length of value should always be maxlen, which should be constant for
            #  all data records.  Therefore, as avg_dist varies in size, value will be zero-
            #  padded at the end (i.e., the right hand side).  This may not work if avg_dist
            #  increases, then decreases.
            @valu = (0) x $maxlen;
            @valu[0..$#avg_dist] = @avg_dist;

            #  Print data to a temp file.
            #if (($#avg_dist+1) < 1) {
            #  print MYTEMP "$lat[$nrec] $lon[$nrec] $msw_kt[$nrec] $rmw[$nrec] $rmw0 $rmw1 $rmw2 0.0 0.0\n";
            #} elsif (($#avg_dist+1) == 1) {
            #  print MYTEMP "$lat[$nrec] $lon[$nrec] $msw_kt[$nrec] $rmw[$nrec] $rmw0 $rmw1 $rmw2 @avg_dist 0.0\n";
            #} else {
            #print MYTEMP "$lat[$nrec] $lon[$nrec] $msw_kt[$nrec] $rmw[$nrec] $rmw0 $rmw1 $rmw2 @avg_dist\n";
            #}
            #print MYTEMP "$lat[$nrec] $lon[$nrec] $msw_kt[$nrec] $rmw[$nrec] $rmw0 $rmw1 $rmw2 @valu\n";

            #  Now increment the index & clear local arrays.
            #my @increment_line = split ' ', $_;
            #if ($_ =~ /REPEAT POSIT:/ ||  $_ =~ /VECTOR TO/ ||  $_ =~ /EXTRATROPICAL$/) { print "----Incrementing NREC; INPUT = $_\n"; $nrec++; print "----Now, NREC = $nrec\n"; }
            if ($_ =~ /REPEAT POSIT:/ ||  $_ =~ /VECTOR TO/ ||  $line[0] eq "EXTRATROPICAL") {
                my @prev = split ' ', $previous;
                if ($DEBUG) { warn "----Incrementing NREC; Prev = $prev[0], INPUT = $_"; }
                #if ($prev[0] ne "EXTRATROPICAL") {
                $nrec++;
                if ($DEBUG) { warn "----Now, NREC = $nrec\n"; }
                #}
            }
            #$nrec++;
            $irad = 0;
            @rad_speed = ();
            @dist = ();
            @sorted_dist = ();
            @avg_dist = ();
            $rmw0 = 0; $rmw1 = 0; $rmw2 = 0; $rmw_orig = 0;
        }  #  end if ($_ =~ /REPEAT_POSIT....
        $previous = $_;

    }  # end while...
    close(INPUT);
    #  close(MYTEMP);

    return;
}

#*******************************************************************************
sub get_year_month {
    #
    #  Revision History:
    #  18 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  18 May 2012  Added determination of storm ID string, and return of ID; moved the
    #               call to gmtime() to outside the nested IF block of the main ELSE
    #               block.  (RSL)
    #
    #*******************************************************************************

    #  Parse the input argument (input file name).
    my ($fname) = @_;

    #  Declare the date & id strings.
    my $date = '';
    my $id = '';

    #  Check the length of the input file name.
    if (length($fname) == 26) {
        #  It is of the form wpNNYYweb_YYYYMMDDHHMM.txt, so extract directly from
        #  the file name.
        $date = substr($fname, 10, 6);
        #  Extract BBNN and YYYY from the file name
        $id = substr($fname, 0, 4) . substr($fname, 10, 4);
    } else {
        #  Get the current system (GMT/UTC) month, year
        my (undef,undef,undef,undef,$mon,$year,undef,undef,undef) = gmtime;
        #  Original wpNNYYweb.txt or unknown file type -- use current system date in UTC.
        if ($fname =~ /wp(\d){4}(\w)+/) {
            #  File name is of the original wpNNYY* format
            $date =  sprintf("%04d%02d", (substr($fname, 4, 2) + 2000), ++$mon);
            #  Extract BBNN from the file name and create the year from the YY
            $id = sprintf("%s%04d", substr($fname, 0, 4), substr($fname, 4, 2) + 2000);
        } else {
            #  An unknown file name format.
            $date =  sprintf("%04d%02d", ($year + 1900), ++$mon);
            #  Unknown name, use "xx" for the basin abbrev, use 98 for number, and add system year
            $id = sprintf("xx98%04d", ($year + 1900));
        }

    }
    return($date,$id);
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
    #  The spreadsheet file "Knaff_Zehr_WPR.xls" (renamed from "Knaff_Zehr WPR.xls")
    #  was received via email on 16 Feb 2012 from:
    #  LT Matthew Watts
    #  Decision Support Operations Officer
    #  Joint Typhoon Warning Center
    #  Pearl Harbor, Hawaii
    #  808-471-4595
    #  matthew.watts@navy.mil
    #  matthew.watts@navy.smil.mil
    #
    #  Syntax: $mslp = wind_pressure_relationship($msw);
    #  where:  $mslp is the returned MSLP value (in hPa), and
    #          $msw is the user-supplied MSW value (in kt)
    #
    #  Calls:  interp (subroutine)
    #
    #  Called by: main program
    #
    #  Revision History:
    #  18 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  22 Jul 2011  Added an ELSIF block to extrapolate to wind speeds below the
    #               minimum wind speed of the table.  (RSL)
    #  02 Mar 2012  Revised the @WPR array to use the values provided by JTWC personnel
    #               (see above comments).  (RSL)
    #
    #*******************************************************************************

    my ($msw) = @_;
    my $mslp = 0;

    #  Define the wind-pressure relationship for the WestPac basin,
    #  based on Dvorak (1984) as summarized by Velden et al. (2006).
    #  The 1st column is maximum steady wind speed (MSW) in knots;
    #  the 2nd column is minimum sea level pressure (MSLP) in hPa (= mb).
    #my @WPR = (
    #  [30, 1000],
    #  [35,  997],
    #  [45,  991],
    #  [55,  984],
    #  [65,  976],
    #  [77,  966],
    #  [90,  954],
    # [102,  941],
    # [115,  927],
    # [127,  914],
    # [140,  898],
    # [155,  879],
    # [170,  858],
    #);

    #  New values from LT Watts of JTWC, in the spreadsheet file "Knaff_Zehr_WPR.xls",
    #  received via email on 16 Feb 2012.  These values are on the 2nd sheet of that
    #  spreadsheet file.
    my @WPR = (
        [15.00, 1010.510425],
        [16.00, 1009.777727],
        [17.00, 1009.044856],
        [18.00, 1008.311811],
        [19.00, 1007.578592],
        [20.00, 1006.845199],
        [21.00, 1006.111632],
        [22.00, 1005.377891],
        [23.00, 1004.643976],
        [24.00, 1003.909887],
        [25.00, 1003.175624],
        [26.00, 1002.441187],
        [27.00, 1001.706576],
        [28.00, 1000.971790],
        [29.00, 1000.236831],
        [30.00, 999.501698],
        [31.00, 998.766391],
        [32.00, 998.030910],
        [33.00, 997.295255],
        [34.00, 996.559426],
        [35.00, 995.823423],
        [36.00, 995.087245],
        [37.00, 994.350894],
        [38.00, 993.614369],
        [39.00, 992.877670],
        [40.00, 992.140797],
        [41.00, 991.403750],
        [42.00, 990.666528],
        [43.00, 989.929133],
        [44.00, 989.191564],
        [45.00, 988.453821],
        [46.00, 987.715904],
        [47.00, 986.977813],
        [48.00, 986.239547],
        [49.00, 985.501108],
        [50.00, 984.762495],
        [51.00, 984.023708],
        [52.00, 983.284746],
        [53.00, 982.545611],
        [54.00, 981.806302],
        [55.00, 981.066819],
        [56.00, 980.327162],
        [57.00, 979.587330],
        [58.00, 978.847325],
        [59.00, 978.107146],
        [60.00, 977.366793],
        [61.00, 976.626265],
        [62.00, 975.885564],
        [63.00, 975.144689],
        [64.00, 974.403640],
        [65.00, 973.662416],
        [66.00, 972.921019],
        [67.00, 972.179448],
        [68.00, 971.437703],
        [69.00, 970.695783],
        [70.00, 969.953690],
        [71.00, 969.211423],
        [72.00, 968.468981],
        [73.00, 967.726366],
        [74.00, 966.983577],
        [75.00, 966.240614],
        [76.00, 965.497476],
        [77.00, 964.754165],
        [78.00, 964.010680],
        [79.00, 963.267020],
        [80.00, 962.523187],
        [81.00, 961.779180],
        [82.00, 961.034998],
        [83.00, 960.290643],
        [84.00, 959.546114],
        [85.00, 958.801410],
        [86.00, 958.056533],
        [87.00, 957.311482],
        [88.00, 956.566256],
        [89.00, 955.820857],
        [90.00, 955.075284],
        [92.00, 953.583615],
        [93.00, 952.837519],
        [94.00, 952.091250],
        [95.00, 951.344807],
        [96.00, 950.598189],
        [97.00, 949.851398],
        [98.00, 949.104432],
        [99.00, 948.357293],
        [100.00, 947.609980],
        [101.00, 946.862492],
        [102.00, 946.114831],
        [103.00, 945.366995],
        [104.00, 944.618986],
        [105.00, 943.870803],
        [106.00, 943.122445],
        [107.00, 942.373914],
        [108.00, 941.625208],
        [109.00, 940.876329],
        [110.00, 940.127275],
        [111.00, 939.378048],
        [112.00, 938.628646],
        [113.00, 937.879071],
        [114.00, 937.129322],
        [115.00, 936.379398],
        [116.00, 935.629301],
        [117.00, 934.879029],
        [118.00, 934.128584],
        [119.00, 933.377964],
        [120.00, 932.627171],
        [121.00, 931.876203],
        [122.00, 931.125062],
        [123.00, 930.373746],
        [124.00, 929.622257],
        [125.00, 928.870593],
        [126.00, 928.118756],
        [127.00, 927.366744],
        [128.00, 926.614559],
        [129.00, 925.862199],
        [130.00, 925.109666],
        [131.00, 924.356958],
        [132.00, 923.604077],
        [133.00, 922.851021],
        [134.00, 922.097791],
        [135.00, 921.344388],
        [136.00, 920.590810],
        [137.00, 919.837059],
        [138.00, 919.083133],
        [139.00, 918.329034],
        [140.00, 917.574760],
        [141.00, 916.820313],
        [142.00, 916.065691],
        [143.00, 915.310895],
        [144.00, 914.555926],
        [145.00, 913.800782],
        [146.00, 913.045465],
        [147.00, 912.289973],
        [148.00, 911.534307],
        [149.00, 910.778468],
        [150.00, 910.022454],
        [151.00, 909.266267],
        [152.00, 908.509905],
        [153.00, 907.753369],
        [154.00, 906.996660],
        [155.00, 906.239776],
        [156.00, 905.482718],
        [157.00, 904.725487],
        [158.00, 903.968081],
        [159.00, 903.210502],
        [160.00, 902.452748],
        [161.00, 901.694820],
        [162.00, 900.936719],
        [163.00, 900.178443],
        [164.00, 899.419993],
        [165.00, 898.661370],
        [166.00, 897.902572],
        [167.00, 897.143600],
        [168.00, 896.384455],
        [169.00, 895.625135],
        [170.00, 894.865641],
        [171.00, 894.105974],
    );

    #  Loop over wind speeds....
    for my $i (0..$#WPR) {
        #  Check whether this MSW value equals the user-supplied value.
        if ($WPR[$i][0] == $msw) {
            #  The values are equal.  Simply return the corresponding MSLP.
            $mslp = $WPR[$i][1];
        } elsif ($i > 0 && ($WPR[$i][0] > $msw && $WPR[$i-1][0] < $msw) ||
            ($i == $#WPR && $WPR[$i][0] < $msw)) {
            #  Interpolate between wind speeds, or extrapolate above the max. wind speed.
            $mslp = interp($WPR[$i][0],$WPR[$i-1][0],$msw,$WPR[$i][1],$WPR[$i-1][1]);
        } elsif ($i == 0 && $WPR[$i][0] > $msw) {
            #  Extrapolate to low wind speed below the minimum.
            $mslp = interp($WPR[$i+1][0],$WPR[$i][0],$msw,$WPR[$i+1][1],$WPR[$i][1]);
        } else {
            #  Some sort of fail-safe may go here.
        }
    }  # end for my $i....

    return $mslp;
}

#*******************************************************************************
sub find_best_value {
    #
    #  Revision History:
    #  21 Jul 2011  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  22 Jul 2011  Added def. of $div (divisor for determining minimum value).  (RSL)
    #
    #*******************************************************************************

    my ($val, $failsafe, $v1, $v2, $v3) = @_;
    my @candidates = ();
    my ($valu, $cand, $sum);
    my $div = 3;  # Started with 4
    my $small = $val/$div;
    my $good = -9999;

    foreach $valu ($v1, $v2, $v3) {
        if ($valu > $val || $valu < $small) {
            $valu = 0;
        } else {
            push @candidates, ($valu);
        }
    }

    $cand = @candidates;
    print "Candidates has $cand values: @candidates.\n";

    if ($cand == 0) {
        $good = $failsafe;
    } elsif ($cand == 1) {
        $good = $candidates[0];
    } else {
        #  Mean of candidate values.
        foreach $valu (@candidates) { $sum += $valu; }
        #$good = $candidates[0];
        $good = $sum/$cand;
    }

    return $good;
}

#*******************************************************************************
sub check_longitude {
    #
    #  This subroutine checks the array of longitudes for any that have opposite signs
    #  (i.e., cross the date line from +180 to -180 or vice versa), and change the
    #  hemisphere of any that are found.  Changes are made according to the signs of
    #  the majority of longitudes; e.g., if most are positive but a few are negative,
    #  then the negative values are changed.
    #
    #  Syntax: check_longitude();
    #
    #  Calls: [No external routines are used.]
    #
    #  Called by: main program
    #
    #  Revision History:
    #  09 Feb 2012  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #
    #*******************************************************************************

    my ($i, $long, @index_pos, @index_neg);
    my ($npos, $nneg) = (0) x 2;

    #  Loop over longitudes....
    for $i ( 0..$#lon ) {
        #  Check the sign of this longitude.
        if ($lon[$i] >= 0) {
            #  Positive; increment the counter & store the index.
            $npos++;
            push @index_pos, ($i);
        } else {
            #  Negative; increment the counter & store the index.
            $nneg++;
            push @index_neg, ($i);
        }
    }

    #  Check whether there are both signs in the longitudes.
    if ($npos > 0 && $nneg > 0) {
        #  There are, so check which sign is more frequent.
        if ($npos > $nneg) {
            #  More positive values, so translate the negatives.
            #  Loop over negative longitudes....
            for $i ( 0..$#index_neg) {
                #  Convert to the other hemisphere.
                $lon[$index_neg[$i]] += 360;
            }
        } else {
            #  More negative values, so translate the positives.
            #  Loop over positive longitudes....
            for $i ( 0..$#index_pos) {
                #  Convert to the other hemisphere.
                $lon[$index_pos[$i]] -= 360;
            }
        }
    }

    return;
}

#*******************************************************************************
sub read_best_track {
    #
    #  read_best_track
    #
    #  This Perl subroutine reads a tropical cyclone (TC) best track file from the
    #  Joint Typhoon Warning Center (JTWC).  A description of the file format can be
    #  found at this URL:
    #
    #  http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/wpindex.html
    #
    #  the contents of which were used as guidance in the development of this routine.
    #
    #  Syntax: (*envpress, *latt, *long, *RMW, *MSLP, *datestr, *max_wind_kt, *storm_name,
    #           *R34, *R34_Rad, *R50, *R50_Rad, *R64, *R64_Rad, *bt_dir, #bt_speed, *bt_datenum) = read_best_track($infile);
    #  where:  $envpress is the returned environmental pressure (mbar),
    #          @latt is the returned array of latitudes (decimal degrees),
    #          @long is the returned array of longitudes (decimal degrees),
    #          @RMW is the returned array of radii of maximum winds (nmi),
    #          @MSLP is the returned array of minimum sea level pressure (mbar),
    #          @datestr is the returned array of date/time values (YYYYMMDDHHMM),
    #          @max_wind_kt is the returned array of maximum sustained wind speed in knots,
    #          @storm_name is the returned array of storm names,
    #          @R34 is the the returned array of 34-knot radius wind speeds (i.e., 34) [0s if not present],
    #          @R34_Rad is the returned array of 34-knot radii values [0s if not present],
    #          @R50 is the the returned array of 50-knot radius wind speeds (i.e., 50) [0s if not present],
    #          @R50_Rad is the returned array of 50-knot radii values [0s if not present],
    #          @R64 is the the returned array of 64-knot radius wind speeds (i.e., 64) [0s if not present],
    #          @R64_Rad is the returned array of 64-knot radii values [0s if not present],
    #          @bt_dir is the returned array of storm directions,
    #          @bt_speed is the returned array of storm translational speeds,
    #          @bt_datenum is the returned array of date numbers, and
    #          $infile is the name of the input best track file.
    #
    #  Calls: sub. storm_dir_speed; Date::Calc::Mktime
    #
    #  Called by: main program
    #
    #  Revision History:
    #  21 Feb 2012  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  22 Feb 2012  Modified the defs. of @lat, @lon to replace the length() calls
    #               with "-1" in the substr() calls; added comments where needed.  (RSL)
    #  24 Feb 2012  Added in IF block to check for the number of fields when obtaining
    #               the storm name; changed the 'DB' field of the %TC_TYPE hash to
    #               'Tropical Disturbance'; copied from tst_best_track.pl.  (RSL)
    #  03 May 2012  Added variables for R34/R50/R64 wind speed and radii parameters;
    #               added a flag and an IF block to check whether the current record is
    #               a duplicate of the previous according to the date/time string;
    #               added IF blocks to read, store the R34/... parameters; added the
    #               R34... parameters to the return() call; updted the help content.
    #               (RSL)
    #  04 May 2012  Added parameters for the storm translational speed & direction and
    #               date number; added a call to sub. storm_dir_speed; added a call to
    #               Mktime for the date number.  (RSL)
    #  07 May 2012  Added a line to skip records with too few parameters.  (RSL)
    #  08 May 2012  Replaced zero with $MISS_D3D (1.0e30), which is the Delft3D missing
    #               value code, where arrays were initialized with place-holders; added
    #               an IF block to check whether the R34 values are zero, and if so,
    #               to replace those with $MISS_D3D.  (RSL)
    #
    #*******************************************************************************

    #  Parse the input parameter list.
    my ($infile) = @_;

    #  Declare, initialize the local parameters.
    my $starz = ('*') x 80;
    my (@lat, @lon, @rmw, @bgp, @mslp, @vmax_kt, @dtg, $date, $time, @name, @line, @sdir, @sspd);
    my ($i, $nrec, @r34, @r34_dist, @r50, @r50_dist, @r64, @r64_dist, $dir, $spd, @dnum);
    my $BGPRESS = 1013.25;  # Default background (AKA environmental) pressure, mbar
    my $MISS = -9999;       # Missing value code
    my $MISS_D3D = 1.0e30;  # Missing value code for Delft3D
    my $NONAME = "NONAME";  # Storm name if none is given, or not a TC
    my $isdup = 0;          # Flag for duplicate record (i.e., has R50 &/or R64 record)

    #  Declare, initialize the environmental (i.e., background) pressure parameters.
    my ($num, $sum, $env_press) = (0) x 3;

    #  Define a hash of the known storm types with the abbreviation as the key:
    my %TC_TYPE = (
        TD => 'Tropical Depression',
        TS => 'Tropical Storm',
        TY => 'Typhoon',
        ST => 'Super Typhoon',
        HU => 'Hurricane',
        SD => 'Subtropical Depression',
        SS => 'Subtropical Storm',
        EX => 'Extropical System',
        IN => 'Inland',
        DS => 'Dissipating',
        LO => 'Low',
        WV => 'Tropical Wave',
        ET => 'Extrapolated',
        XX => 'Unknown',
        DB => 'Tropical Disturbance',   # This type is undocumented, but is common in the files.
    );

    #  Open the input file.
    open (INPUT, "<$infile") or die "ERROR: Unable to open input file $infile: $!; stopped";

    #  Loop over input file records....
    while (<INPUT>) {
        s/\s//g;                # Remove the spaces in this line.
        @line = split ',', $_;  # Split this comma-delimited record into fields
        next if $#line < 17;    # Skip records with too few parameters.

        #if ($. < 11) { print "Line = "; for my $w (@line) {print " $w\n";}; print "$starz\n"; }
        #  Check whether this is a duplicate record (i.e., has an R50/R64 for the previous record
        #  that had an R34.
        if (scalar @dtg > 1 && $line[2] == $dtg[$#dtg] && $line[11] > 34) {
            $isdup = 1;
        }

        #  Check whether this is a duplicate record.
        if ($isdup == 0) {
            #  It is not, so load the arrays w/the data.

            #  Store the date/time group.
            push @dtg, ($line[2]);
            push @dnum, (Mktime(substr($line[2],0,4)+0,substr($line[2],4,2)+0,substr($line[2],6,2)+0,substr($line[2],8,2)+0,0,0));


            #  Check latitude value (tenths of degrees) to get the correct sign according to the hemisphere.
            if (substr($line[6],-1) eq 'S') {
                push @lat, (substr($line[6], 0, -1)/10 * -1);
            } else {
                push @lat, (substr($line[6], 0, -1)/10);
            }

            #  Check longitude value (tenths of degrees) to get the correct sign according to the hemisphere.
            if (substr($line[7],-1) eq 'W') {
                push @lon, (substr($line[7], 0, -1)/10 * -1);
            } else {
                push @lon, (substr($line[7], 0, -1)/10);
            }

            #  Determine the storm translation speed & direction.
            if (scalar @dtg > 1) {
                ($dir, $spd) = storm_dir_speed($lon[$#lon-1],$lat[$#lat-1],$lon[$#lon],$lat[$#lat],0,6);
                push @sdir, ($dir);
                push @sspd, ($spd);
            } else {
                push @sdir, ($MISS_D3D);
                push @sspd, ($MISS_D3D);
            }

            #  Maximum sustained wind speed in knots: 0 through 300.
            push @vmax_kt, ($line[8]);

            #  Minimum sea level pressure, 1 through 1100 MB.
            push @mslp, ($line[9]);

            #  R34/50/64 parameters.
            #  Store the R34 wind speed (34 or zero)
            push @r34, ($line[11]);

            #  Get the current number of valid records
            $nrec = $#r34;

            #  Load the R50, R64 wind speed values with place holders.
            $r50[$nrec] = $MISS_D3D;
            $r64[$nrec] = $MISS_D3D;

            #  Loop over the the R34 quadrants to store the radii & load the R50, R64 radii with place holders.
            for $i (0..3) {
                #  Check whether a zero is present.
                if ($line[11+$i+2] == 0) {
                    #  A zero is present, so store the missing value code.
                    $r34_dist[$nrec][$i] = $MISS_D3D;
                } else {
                    $r34_dist[$nrec][$i] = $line[11+$i+2];
                }
                $r50_dist[$nrec][$i] = $MISS_D3D;
                $r64_dist[$nrec][$i] = $MISS_D3D;
            }

            #  Add the remaining parameters:
            push @bgp, ($line[17]);   # Pressure in millibars of the last closed isobar, 900 - 1050 mb.
            push @rmw, ($line[19]);   # Radius of max winds, 0 - 999 nm.
            #  Check whether the storm name should be present.
            if ($#line >= 27) {
                #  It should be present.
                push @name, (qq($TC_TYPE{$line[10]} $line[27]));  # Storm name
            } else {
                #  Too few fields, so use the default name.
                push @name, ($NONAME);  # Name
            }

        } else {
            #  It is a duplicate record, so load only the R50 or R64 parameters.
            #  Get the current number of valid records
            $nrec = $#r34;
            #  Check which RNN value is being read.
            if ($line[11] == 50) {
                #  R50, so store that wind speed.
                $r50[$nrec] = $line[11];
                #  Loop over quadrants to store the R50 radii.
                for $i (0..3) {
                    $r50_dist[$nrec][$i] = $line[11+$i+2];
                }

            } elsif ($line[11] == 64) {
                #  R64, so store that wind speed.
                $r64[$nrec] = $line[11];
                #  Loop over quadrants to store the R64 radii.
                for $i (0..3) {
                    $r64_dist[$nrec][$i] = $line[11+$i+2];
                }
            }
        }
        #  Reset the duplicate flag.
        $isdup = 0;
    }
    close(INPUT);

    #  Calculate the mean of the "background" or environmental pressure (if any).
    #  Loop over values....
    for my $i (0..$#bgp) {
        #  Check validity.
        if ($bgp[$i] > $MISS) {
            #  A valid pressure value, so increment the counter and accumulator.
            $num++;
            $sum += $bgp[$i];
        }
    }
    #  Check for the number of valid values.
    if ($num > 0) {
        #  At least 1, so calculate the mean pressure of the last closed isobar.
        $env_press = $sum/$num;
    } else {
        #  No valid values, so use the default background pressure.
        $env_press = $BGPRESS;
    }

    #  print "lats = @lat\n";
    #  print "Env Press = $env_press\n";
    #  print "lons = @lon\n";
    #  print "names = @name\n";

    #(*latt, *long, *RMW, *MSLP, *datestr, *Max_Wind, *StormName) = read_best_track($infile,$outfile);
    return(\$env_press, \@lat, \@lon, \@rmw, \@mslp, \@dtg, \@vmax_kt, \@name,
        \@r34, \@r34_dist, \@r50, \@r50_dist, \@r64, \@r64_dist, \@sdir, \@sspd, \@dnum);
}

#*******************************************************************************
sub combine_data {
    #
    #  combine_data
    #
    #  This subroutine combines data sets from a standard JTWC warning file and a
    #  best track file in chronological order.  The output consists of the original
    #  warning array names with the entire data series in each.
    #
    #  Syntax: combine_data();
    #
    #  Calls: [No external routines are used.]
    #
    #  Called by: main program
    #
    #  Revision History:
    #  04 May 2012  Initial coding.  R.S. Linzell, QNA/NRL Code 7322
    #  08 May 2012  Added the array @r100_65_50_35_3D_tmp, and the parameters $MISS,
    #               $nquad, $indx to provide support for the R34/R50/R64 wind radii
    #               parameters. (RSL)
    #  25 Feb 2013  Added $MISS_EBT, @pres_outer, @rad_outer, @eyed to handle new
    #               extended best track parameters; added lines in the IF block within
    #               the FOR loop to load those new parameters.  (RSL)
    #  08 Mar 2013  Added defs. $time_best_track (date number of current best track
    #               record), $start_time_data (date number of 1st TC warning record),
    #               & added an IF block to use only those best track records that
    #               precede the 1st TC warning record (i.e., the analysis); moved
    #               incrementing of $indx to correctly increment it when the output
    #               arrays are appended.  (RSL)
    #
    #*******************************************************************************

    #  Create local arrays that will store the corresponding input warning arrays.
    my ($i, $key, @date_tmp, @time_tmp, @lat_tmp, @lon_tmp, @mslp_tmp, @msw_kt_tmp);
    my (@rmw_tmp, @tc_dir_tmp, @tc_speed_tmp, @datenum_tmp, @r100_65_50_35_3D_tmp);
    my (@pres_outer, @rad_outer, @eyed, $time_best_track, $start_time_data);
    my $MISS = 1.0e30;
    my $MISS_EBT = -99;
    my $nquad = 4;
    my $indx = 0;

    #  Copy the input warning arrays.
    @date_tmp = @date;
    @time_tmp = @time;
    @lat_tmp = @lat;
    @lon_tmp = @lon;
    @mslp_tmp = @mslp;
    @msw_kt_tmp = @msw_kt;
    @rmw_tmp = @rmw;
    @tc_dir_tmp = @tc_dir;
    @tc_speed_tmp = @tc_speed;
    @datenum_tmp = @datenum;
    @r100_65_50_35_3D_tmp = @r100_65_50_35_3D;
    @pres_outer = @pres_outer_iso_w;
    @rad_outer = @rad_outer_iso_w;
    @eyed = @eye_diam_w;

    #  Clear the input warning arrays.
    @date = (); @time = (); @lat = (); @lon = (); @mslp = (); @msw_kt = (); @rmw = ();
    @tc_dir = (); @tc_speed = (); @datenum = (); @r100_65_50_35_3D = ();
    @eye_diam_w = (); @rad_outer_iso_w = (); @pres_outer_iso_w = ();

    #  Create the date number for the start (i.e., analysis) of the TC warning.
    #warn "Best Track datestr = @datestr\n  TC Warning date = @date_tmp\n  TC Warning time = @time_tmp\n";
    $start_time_data = Mktime(substr($date_tmp[0],0,4),substr($date_tmp[0],4,2),substr($date_tmp[0],6,2),
        substr($time_tmp[0],0,2),substr($time_tmp[0],2,2),0);
    #warn "Start_Time_data = $start_time_data.\n";

    # Loop over sorted dates to store data in the original warning arrays in correct order.
    for $key (sort keys %time_type) {
        #warn "Date: $key  Desc: $time_type{$key}[0]  Index:$time_type{$key}[1] \n";
        #  Index of this record:
        $i = $time_type{$key}[1];
        #  Determine the type of this record.
        if ($time_type{$key}[0] =~ /B/) {
            #  Best track record.
            #  Create the date number for this best track record.
            $time_best_track = Mktime(substr($datestr[$i],0,4),substr($datestr[$i],4,2),substr($datestr[$i],6,2),
                substr($datestr[$i],8,2),0,0);

            #warn " Time_best_track = $time_best_track.\n";

            #  We want only those best track records that occur before the 1st
            #  TC warning record (i.e., the analysis).
            if ($time_best_track < $start_time_data) {
                push @date, (substr($datestr[$i],0,8));
                push @time, (substr($datestr[$i],8,2) . '00');
                push @lat, ($latt[$i]);
                push @lon, ($long[$i]);
                push @mslp, ($MSLP[$i]);
                push @msw_kt, ($max_wind_kt[$i]);
                push @rmw, ($RMW[$i]);
                push @data_type, ('B');
                push @tc_dir, ($bt_dir[$i]);
                push @tc_speed, ($bt_speed[$i]);
                push @datenum, ($bt_datenum[$i]);
                push @pres_outer_iso_w, ($pres_outer_iso[$i]);
                push @rad_outer_iso_w, ($rad_outer_iso[$i]);
                push @eye_diam_w, ($eye_diam[$i]);
                #  Load the r100 array with the corresponding best track parameters.  Indices are,
                #  respectively, record number, wind speed (100, 64/65, 50, 34/35, respectively), and
                #  quadrant (starting at northeast and progressing clockwise).
                push @{ $r100_65_50_35_3D[$indx][0] }, ($MISS) x $nquad;  # R100 (not presently used by JTWC)
                push @{ $r100_65_50_35_3D[$indx][1] }, @{ $R64_Rad[$i] } [0..$nquad-1];  # R64
                push @{ $r100_65_50_35_3D[$indx][2] }, @{ $R50_Rad[$i] } [0..$nquad-1];  # R50
                push @{ $r100_65_50_35_3D[$indx][3] }, @{ $R34_Rad[$i] } [0..$nquad-1];  # R34
                $indx++;
            }
        } else {
            #  Warning record.
            push @date, ($date_tmp[$i]);
            push @time, ($time_tmp[$i]);
            push @lat, ($lat_tmp[$i]);
            push @lon, ($lon_tmp[$i]);
            push @mslp, ($mslp_tmp[$i]);
            push @msw_kt, ($msw_kt_tmp[$i]);
            push @rmw, ($rmw_tmp[$i]);
            push @data_type, 'F';
            push @tc_dir, ($tc_dir_tmp[$i]);
            push @tc_speed, ($tc_speed_tmp[$i]);
            push @datenum, ($datenum_tmp[$i]);
            push @pres_outer_iso_w, ($MISS_EBT);
            push @rad_outer_iso_w, ($MISS_EBT);
            push @eye_diam_w, ($MISS_EBT);
            #  Load the r100 array with the corresponding warning parameters.  Indices are as above.
            push @{ $r100_65_50_35_3D[$indx][0] }, @{ $r100_65_50_35_3D_tmp[$i][0] } [0..$nquad-1];  # R100 (not presently used by JTWC)
            push @{ $r100_65_50_35_3D[$indx][1] }, @{ $r100_65_50_35_3D_tmp[$i][1] } [0..$nquad-1];  # R64
            push @{ $r100_65_50_35_3D[$indx][2] }, @{ $r100_65_50_35_3D_tmp[$i][2] } [0..$nquad-1];  # R50
            push @{ $r100_65_50_35_3D[$indx][3] }, @{ $r100_65_50_35_3D_tmp[$i][3] } [0..$nquad-1];  # R34
            $indx++;
        }
        #$indx++;
    }

    return;
}

#*******************************************************************************
sub storm_dir_speed {
    #
    #  storm_dir_speed
    #
    #  This subroutine calculates the storm direction and translational speed for a given
    #  time step relative to the previous time step.
    #
    #  Syntax: ($dir, $speed) = storm_dir_speed($lon1, $lat1, $lon2, $lat2, $tau1, $tau2);
    #  where:  $dir is the output storm translational direction,
    #          $speed is the output storm translational speed,
    #          $lon1 is the input previous longitude,
    #          $lat1 is the input previous latitude,
    #          $lon2 is the input current longitude,
    #          $lat2 is the input current latitude,
    #          $tau1 is the input previous tau (time step) value, and
    #          $tau2 is the input current tau value.
    #
    #  Calls: Math::Trig::deg2rad, great_circle_distance, great_circle_bearing, rad2deg
    #
    #  Called by: subs. read_jtwc_file, read_best_track
    #
    #  Revision History:
    #  04 May 2012  Initial coding; moved code from sub. read_jtwc_file.  R.S. Linzell,
    #               QNA/NRL Code 7322
    #  22 May 2012  Added an IF block to check whether the input lat & lon values are
    #               the same, and if so, sets the speed and direction to zero.  (RSL)
    #
    #*******************************************************************************

    #  Parse out the input parameters.
    my ($lon1, $lat1, $lon2, $lat2, $tau1, $tau2) = @_;

    #  Create local parameters that will store the corresponding input warning arrays.
    my (@start, @end, $km, $dir, $speed);
    my $nm2km = 1.852;      # Conversion for nmi to km

    #  Check whether the latitudes and longitudes are equal.
    if ($lat1 == $lat2 && $lon1 == $lon2) {
        #  The values are equal, so set speed & dir. to zero.
        $speed = 0.0;
        $dir = 0.0;
    } else {
        #  Starting lon, lat (from previous record) converted to radians:
        @start = (deg2rad($lon1), deg2rad(90 - $lat1));
        #  Ending lon, lat (from the current record), converted to radians:
        @end = (deg2rad($lon2), deg2rad(90 - $lat2));
        #  Great circle distance, in km:
        $km = great_circle_distance(@start, @end, 6378.137);
        #  Bearing from start to end point (radians):
        $dir = great_circle_bearing(@start, @end);
        #  Convert direction to degrees:
        $dir = rad2deg($dir);
        #  TC translation speed is distance / time:
        $speed = ($km/$nm2km)/($tau2 - $tau1);
    }

    return($dir, $speed);

}
