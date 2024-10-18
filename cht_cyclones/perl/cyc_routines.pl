use strict;

#*******************************************************************************
sub read_std_best_track {
    #
    #  read_std_best_track
    #
    #  This Perl subroutine reads a tropical cyclone (TC) best track file from the
    #  Joint Typhoon Warning Center (JTWC).  A description of the file format can be
    #  found at this URL:
    #
    #  http://www.usno.navy.mil/NOOC/nmfc-ph/RSS/jtwc/best_tracks/wpindex.html
    #
    #  the contents of which were used as guidance in the development of this routine.
    #  Best track files in ATCF format from NHC are also supported.
    #
    #  Syntax: (*envpress, *latt, *long, *RMW, *MSLP, *datestr, *max_wind_kt, *storm_name,
    #           *R34, *R34_Rad, *R50, *R50_Rad, *R64, *R64_Rad, *bt_datenum, *pres_outer_iso,
    #           *rad_outer_iso, *storm_id, *eye_diam) = read_std_best_track($all_flag,$infile);
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
    #          @bt_datenum is the returned array of date numbers,
    #          @pres_outer_iso is the returned array of pressure of the outermost close isobar,
    #          @rad_outer_iso is the returned array of radii of the outermost close isobar,
    #          @storm_id is the returned array of storm ID strings,
    #          @eye_diam is the returned array of eye diameter values,
    #          $all_flag is a flag set to 1 to read all records, or to zero for only full records, and
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
    #               Mktime for the date nuber.  (RSL)
    #  07 May 2012  Added a line to skip records with too few parameters.  (RSL)
    #  08 May 2012  Replaced zero with $MISS_D3D (1.0e30), which is the Delft3D missing
    #               value code, where arrays were initialized with place-holders; added
    #               an IF block to check whether the R34 values are zero, and if so,
    #               to replace those with $MISS_D3D.  (RSL)
    #  04 Feb 2013  Copied from parse_jtwc_warning.pl, sub. read_best_track; renamed
    #               to read_std_best_track; removed TC speed & direction calculations;
    #               added reading, storing, returning of pressure & radius of outermost
    #               closed isobar.  (RSL)
    #  07 Feb 2013  Added array @sid to store the storm ID string (e.g., 'AL0112'); updated
    #               the help content.  (RSL)
    #  13 Feb 2013  Added array @eyediam to store & return the storm eye diameter
    #               values; updated the help content.  (RSL)
    #  15 May 2013  Added input parameter $all_flag to control whether all records are
    #               read, or only complete records are read; turned off certain warning
    #               messages that occur when missing data are processed.  (RSL)
    #
    #*******************************************************************************

    #  Parse the input parameter list.
    my ($all_flag,$infile) = @_;

    #  Turn off certain warnings.
    no warnings qw(numeric uninitialized);

    #  Declare, initialize the local parameters.
    my $starz = ('*') x 80;
    my (@lat, @lon, @rmw, @bgp, @mslp, @vmax_kt, @dtg, $date, $time, @name, @line, @rad_outer_iso);
    my ($i, $nrec, @r34, @r34_dist, @r50, @r50_dist, @r64, @r64_dist, $dir, $spd, @dnum, @sid, @eyediam);
    my $BGPRESS = 1013.25;  # Default background (AKA environmental) pressure, mbar
    my $MISS = -999;       # Missing value code
    my $MISS_D3D = 0; # 1.0e30;  # Missing value code for Delft3D
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

        # Skip records with too few parameters if not all of the inputs are wanted.
        if ($#line < 17 && $all_flag == 0) { next; }
        #next if $#line < 17;    

        #if ($. < 11) { print "Line = "; for my $w (@line) {print " $w\n";}; print "$starz\n"; }
        #  Check whether this is a duplicate record (i.e., has an R50/R64 for the previous record
        #  that had an R34.
        if (scalar @dtg > 1 && $line[2] == $dtg[$#dtg] && $line[11] > 34) {
            $isdup = 1;
            #warn " . . . .  Duplicate . . . .\n";
        }

        #  Check whether this is a duplicate record.
        if ($isdup == 0) {
            #  It is not, so load the arrays w/the data.

            #  Store the storm ID (basin, sequential storm number, last 2 digits of year).
            push @sid, ($line[0] . $line[1] . substr($line[2], 2, 2));
            #  Store the date/time group.
            push @dtg, ($line[2]);
            #  Store the calculated date number (sec. since beginning of time)
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
            push @bgp, ($line[17]);            # Pressure in millibars of the last closed isobar, 900 - 1050 mb.
            push @rad_outer_iso, ($line[18]);  # Radius of the last closed isobar in nm, 0 - 9999 nm
            push @rmw, ($line[19]);            # Radius of max winds, 0 - 999 nm.
            push @eyediam, ($line[21]);        # Eye diameter, 0 - 999 nm

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


    return($env_press, \@lat, \@lon, \@rmw, \@mslp, \@dtg, \@vmax_kt, \@name, \@r34, \@r34_dist,
        \@r50, \@r50_dist, \@r64, \@r64_dist, \@dnum, \@bgp, \@rad_outer_iso, \@sid, \@eyediam);
}

#********************************************************************************
sub write_cyc {
    #
    # write_cyc
    # 
    # This perl subroutine writes the cyc file from the best track file and the 
    # forecast file
    #
    # Syntax: (
    #
    #
    # Revision History:
    # 19 May 2022  Initial coding. Jay Veeramony, NRL 7322
    #
    # 
    #*******************************************************************************

    #  Parse the input parameters.  Note that arrays are input as references
    #  (hence the "$$" notation below).
    my ($trkfile,$name,$eye_diam,$envpress,$rad_outer_iso,$year,$month,$day,$time,
        $lat,$lon,$msw,$rmax,$r100,$mslp,$SpiderwebRadius, $NrRadialBins, $NrDirectionalBins ) = @_;

    #  Turn off certain warnings.
    no warnings qw(numeric uninitialized);

    # Set some defaults (These should be set outside of here actually)
    my $wind_profile = "holland2010";
    my $wind_pressure_relation = "holland2008";
    my $rmax_relation = "gross2004";
    my $wind_conversion_factor = 1.0;
    my $background_pressure = 1012;
    my $phi_spiral = 20;

    #  Initialize the indices of the R34/50/64 radii within the 3-D matrix.
    my ($i34,$i50,$i64,$sum34,$sum50,$sum64) = (0) x 6;
    my $i100 = 0;

    #  Define the missing value code (per above documentation).
    my $MISS = -999;

    #  Define the output format.
    my $FMT ="%-8s %-6s %8.3f %8.3f %12.1f %10.1f %10.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f \n";

    #  Get only the storm proper name (e.g., ALICE from Hurricane Alice).
    my @sname = split ' ', $name;
    my $stname = uc($sname[$#sname]);

    #  Open the output file.
    open(TRK, ">$trkfile") or die "ERROR: Unable to open output file $trkfile ($!); stopped";

    #  Print the header data
    printf TRK "Name                   \"" . $stname . "\"\n";
    printf TRK "WindProfile            " . $wind_profile . "\n";
    printf TRK "WindPressureRelation   " . $wind_pressure_relation . "\n";
    printf TRK "RMaxRelation           " . $rmax_relation . "\n";
    printf TRK "Backgroundpressure     " . $background_pressure . "\n";
    printf TRK "PhiSpiral              " . $phi_spiral . "\n";
    printf TRK "WindConversionFactor   " . $wind_conversion_factor . "\n";
    printf TRK "SpiderwebRadius        " . $SpiderwebRadius . "\n";
    printf TRK "NrRadialBins           " . $NrRadialBins . "\n";
    printf TRK "NrDirectionalBins      " . $NrDirectionalBins . "\n";
    printf TRK "#\n";
    printf TRK "#   Date   Time       Lat      Lon   Vmax (kts)   Pc (hPa)  Rmax (NM)  R35(NE)  R35(SE)  R35(SW)  R35(NW)  R50(NE)  R50(SE)  R50(SW)  R50(NW)  R65(NE)  R65(SE)  R65(SW)  R65(NW) R100(NE) R100(SE) R100(SW) R100(NE) \n";
    printf TRK "#\n";

    #  Loop over the number of input data records....
    for my $i (0..$#$lat) {

        #  Determine which of the R34/50/64 wind radii are nonzero by summing each.
        $sum34 = $$r100[$i][3][0] + $$r100[$i][3][1] + $$r100[$i][3][2] + $$r100[$i][3][3];
        $sum50 = $$r100[$i][2][0] + $$r100[$i][2][1] + $$r100[$i][2][2] + $$r100[$i][2][3];
        $sum64 = $$r100[$i][1][0] + $$r100[$i][1][1] + $$r100[$i][1][2] + $$r100[$i][1][3];

        #  Define the indices of the R34/50/64 wind radii, depending on which values are nonzero.
        if ($sum64 == 0 && $sum50 == 0 && $sum34 >= 0) {
            #  R34 values may be present, R50 & R64 are not.
            $i34 = 3;
            $i50 = 2;
            $i64 = 1;
            $i100= 0;
        } elsif ($sum64 == 0 && $sum50 > 0 && $sum34 > 0) {
            #  R34 & R50 values are present, R64 are not.
            $i34 = 2;
            $i50 = 3;
            $i64 = 1;
            $i100= 0;
        } elsif ($sum64 > 0 && $sum50 > 0 && $sum34 > 0) {
            #  All 3 radii are present.
            $i34 = 1;
            $i50 = 2;
            $i64 = 3;
            $i100= 0;
        }

        #  Replace zeros with missing value code.
        if ($$msw[$i] == 0)           { $$msw[$i] = $MISS; }
        if ($$mslp[$i] == 0)          { $$mslp[$i] = $MISS; }
        if ($$rmax[$i] == 0)          { $$rmax[$i] = $MISS; }
        if ($$eye_diam[$i] == 0)      { $$eye_diam[$i] = $MISS; }
        if ($$envpress[$i] == 0)      { $$envpress[$i] = $MISS; }
        if ($$rad_outer_iso[$i] == 0) { $$rad_outer_iso[$i] = $MISS; }

        if (length($$month[$i]) == 1) { $$month[$i] = "0" . $$month[$i]; }
        if (length($$day[$i])   == 1) { $$day[$i]   = "0" . $$day[$i];   }
        if (length($$time[$i])  == 1) { 
            $$time[$i] = "0" . $$time[$i] . "0000";
        }  elsif ( length($$time[$i]) == 2 ) {
            $$time[$i] = $$time[$i] . "0000";
        }  elsif ( length($$time[$i]) == 4 ) {
            $$time[$i] = $$time[$i] . "00";
        }

        my $yyyymmdd = $$year[$i] . $$month[$i] . $$day[$i];
        my $vmax = $$msw[$i];

        printf TRK $FMT, $yyyymmdd, $$time[$i],
        $$lat[$i], $$lon[$i], $vmax, $$mslp[$i], $$rmax[$i],
        $$r100[$i][$i34][0], $$r100[$i][$i34][1] ,$$r100[$i][$i34][2], $$r100[$i][$i34][3],
        $$r100[$i][$i50][0], $$r100[$i][$i50][1], $$r100[$i][$i50][2], $$r100[$i][$i50][3],
        $$r100[$i][$i64][0], $$r100[$i][$i64][1], $$r100[$i][$i64][2], $$r100[$i][$i64][3],
        $$r100[$i][$i100][0], $$r100[$i][$i100][1], $$r100[$i][$i100][2], $$r100[$i][$i100][3];
 
    }
        
    #  Close the output file.
    close(TRK);

    return;
}
1;
