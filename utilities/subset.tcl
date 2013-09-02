#!/users/modellers/pica/Software/bin/convsh

# Convsh script subset.tcl
#
# Takes one file in (e.g. GRIB or NetCDF or PP) and writes out
# a regular gridded subset in x, y, and z.
#
# subset.tcl -i infile -o outfile -of -zm -mm outformat -xs xstart -xe xend -xi xinc
#               -ys ystart -ye yend -yi yinc -levs levs
#
# -i : input file name(s) (and path) - can be more than one.
# -o : output file name (and path)
# -of: output file format (netcdf, grads, drs or utf) - Defaults to NetCDF.
# -zm: create zonal means.
# -mm: create meridional means.
# -p:  parameters selected.
# -xs : x-dimension start - i.e. first grid point (usually 0)
# -xe : x-dimension end - i.e. last grid point (usually 359)
# -xi : x-dimension increment - i.e. spacing between grid points (usually 1)
# -ys, -ye and -yi are the same as above but for the y-dimension.
# -levs : selected output levels (e.g. from model levels,
#         pressure levels or potential temperature levels).
#
# Version 1.0: Written by Ag Stephens, BADC, 22 April 2002.
#       Designed for use with the BADC archive of ECMWF data (e.g. ERA-40).
# Version 1.1: Ag Stephens, BADC, 01 August 2002.
#       Updated to allow level subsetting of BADC archives of ECMWF ERA-15 and Operational data.
# Version 1.2: Ag Stephens, BADC, 30 April 2003.
#       Updated to allow the creation of zonal and meridional means.
#
# Note: to create a subset across the 0 degree Meridian, you should
# use a negative value (i.e. degrees West) as the xstart argument and
# a positive value as the xend argument.
#
# Define the exit procedure if errors encountered
proc exitandreport {problem} {   # Exit and tell user how to use command
  puts "\nExiting because $problem not defined or defined incorrectly...\n"
  puts "Correct usage of 'subset.tcl'"
  puts "-----------------------------"
  puts "subset.tcl -i infile(s) -o outfile \[-zm\] \[-mm\] \[-of outformat\] \[-xs xstart "
  puts "             -xe xend -xi xinc\]  \[-ys ystart -ye yend -yi yinc\] \[-levs levs\]\n"
  puts "Required arguments"
  puts "------------------"
  puts "-i infile(s)  - input file name(s) (and directory paths if desired)."
  puts "              - can be more than one file, globs allowed.\n"
  puts "-o outfile    - output file name (and path).\n"
  puts "-of outformat - output file format (netcdf \[default\], grads, drs or utf).\n"
  puts "Optional arguments"
  puts "------------------"
  puts "-zm        - create zonal means."
  puts "-mm        - create meridional means."
  puts "-xs xstart - x-dimension start - i.e. first grid point."
  puts "-xe xend   - x-dimension end - i.e. last grid point."
  puts "-xi xinc   - x-dimension increment - i.e. spacing between grid points (must be >= 0.5)."
  puts "(If one x-dimension argument is used then all three are required).\n"
  puts "-ys, -ye and -yi are the same as above but for the y-dimension."
  puts "(If one y-dimension argument is used then all three are required).\n"
  puts "-levs levels - selected output levels (e.g. from model levels,"
  puts "               pressure levels or potential temperature levels)."
  puts "             - individual levels can be separated by commas."
  puts "             - for model levels ranges can be used (e.g. 1-7,20-24)."
  puts "-p - parameter codes (can be a comma separated list, e.g. 't,q'). Note that parameter codes that begin with a number are not allowed so you should reference them with the letter character(s) precede the number(s) (e.g. '10u' should be called 'u10')."
  exit
}

#  Write out Netcdf file
set outformat netcdf

#  Automatically work out input file type
set filetype 0

#  Get command line arguments:

set i false
set o false
set of false
set zm false
set mm false
set s false
set xs false; #set xstart "unset"
set xe false
set xi false
set ys false; #set ystart "unset" ; set yend "unset"
set ye false
set yi false
set levs false; set inlevels "unset"
set interp_xy false; set subset_z false;
set params false

foreach arg $argv {
   switch -glob -- $arg {
      -i    {set i true; set o false; set of false; set zm false; set mm false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -o    {set i false; set o true; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -of    {set i false; set o false; set of true; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -zm    {set i false; set o false; set of false; set zm true; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -mm    {set i false; set o false; set of false; set zm false; set mm true; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -xs   {set i false; set o false; set of false; set zm false; set mm false; set xs true; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -xe   {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe true
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params false}
      -xi   {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi true; set ys false; set ye false; set yi false
             set levs false; set params false}
      -ys   {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys true; set ye false; set yi false
             set levs false; set params false}
      -ye   {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye true; set yi false
             set levs false; set params false}
      -yi   {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi true
             set levs false; set params false}
      -levs {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs true; set params false}
      -p    {set i false; set o false; set of false; set zm false; set mm false; set xs false; set xe false
             set xi false; set ys false; set ye false; set yi false
             set levs false; set params true}
      -[a-zA-Z]*  {puts "unknown option $arg"; set i false; set o false
             set of false; set zm false; set mm false; set xs false; set xe false; set xi false; set ys false
             set ye false; set yi false; set levs false; set params false}
      default {
         if {$i} {
            set infile [lappend infile $arg]
         } elseif {$o} {
            set outfile $arg
         } elseif {$of} {
            set outformat $arg
         } elseif {$xs} {
            set xstart $arg
         } elseif {$xe} {
            set xend $arg
         } elseif {$xi} {
            set xinc $arg
         } elseif {$ys} {
            set ystart $arg
         } elseif {$ye} {
            set yend $arg
         } elseif {$yi} {
            set yinc $arg
         } elseif {$levs} {
            set inlevels $arg
         } elseif {$params} {
            set paramarg $arg
         } else {
            puts "unknown option $arg"
         }
      }
   }
}

# Check if necessary arguments have been given
foreach testarg {infile outfile} {
  if {![info exists $testarg]} {
    exitandreport $testarg
  }
}
foreach d {x y} {
  if {[info exists ${d}start]} {
    if {![info exists ${d}end] || ![info exists ${d}inc]} {
      exitandreport "${d}-dimension arguments (${d}s, ${d}e, and ${d}i)"
    }
  }
}
foreach d {x y} {
  if {[info exists ${d}end]} {
    if {![info exists ${d}start] || ![info exists ${d}inc]} {
      exitandreport "${d}-dimension arguments (${d}s, ${d}e, and ${d}i)"
    }
  }
}
foreach d {x y} {
  if {[info exists ${d}inc]} {
    if {![info exists ${d}end] || ![info exists ${d}start]} {
      exitandreport "${d}-dimension arguments (${d}s, ${d}e, and ${d}i)"
    }
    if {$d == "x"} {
      if {$xinc < 0.1} {
        exitandreport "xi argument (which must be greater than 0.1)"
      }
    } elseif {$d == "y"} {
      if {$yinc < 0.1} {
        exitandreport "yi argument (which must be greater than 0.1)"
      }
    }
  }
}

if {[info exists outfile]} {
   if {[llength $outfile] > 1} {
      set outfile [lindex $outfile 0]
      puts "Only one output file can be specified, using $outfile"
   }
} else {
   puts "output file name must be given"
   exit
}

#  Read in each of the input files

foreach file $infile {
   readfile $filetype $file
}

# Set the parameters selected

if {[info exists paramarg]} {
  set parameters [split $paramarg ,]
  set no_fields [llength [listhead]]
  for {set ifield 0} {$ifield < $no_fields} {incr ifield} {
    foreach parameter $parameters {
      set upparam [string toupper $parameter]
      if {[get_shortname $ifield] == $upparam} {
        lappend fieldlist $ifield
      }
    }
  }
} else {
  #  Set all fields if none specified
  set fieldlist -1
}

#  Set all fields if none specified
if {[llength $fieldlist] == 0} {
  set fieldlist -1
}

# Print the selected params
if {$fieldlist == -1} {
  set reportparams all
} else {
  set reportparams $fieldlist
}
puts "Parameter numbers to output are: $reportparams\n"

# Set the output levels in z
if {$inlevels != "unset"} {
  regexp -- {level code = ([0-9]+)} [printhead] dummy levtype
    puts "The level type is: $levtype.\n"
    set levs [split $inlevels ,]
    foreach level $levs {
      if {[regexp {^[0-9]+\-[0-9]+$} $level]} {
          if {$levtype != 109} { # If NOT model levels
            puts "\nPlease only use level ranges in argument (e.g. a-b) for data on model levels.\n";
            puts "Exiting programme.\n\n"
          exit
          }
          set levrange [split $level -]
          set firstlev [lindex $levrange 0]
          set lastlev [lindex $levrange 1]
          incr lastlev
          for {set lcount $firstlev} {$lcount < $lastlev} {incr lcount} {
            set levels [lappend levels $lcount]
          }
      } else {
        set levels [lappend levels $level]
      }
    }
    foreach level $levels {
      if {$levtype == 109} {
        #Should work for any number of model levels e.g. 1-31, 1-60
        set level [expr $level-1]
        set outlevels [lappend outlevels $level]
      } elseif {$levtype == 100} {
        # Get number of pressure levels
        regexp -- {number of selected points in z direction = ([0-9]+)} [printhead] dummy1 numlevs dummy2
        puts "NUMBER OF LEVS is $numlevs\n"
        # Get start and end levels
        regexp -- {z\[0\] = ([0-9]+\.?[0-9]+) z\[([0-9]+)\] = ([0-9]+\.?[0-9]+)} [printhead] dummy1 firstlev dummy2 lastlev
        # Define levels depending on the amount and order.
        if {$numlevs == 15} {
          if {$firstlev == 1000} {
            set levmap "1000 925 850 700 500 400 300 250 200 150 100 70 50 30 10"
          } else {
            set levmap "10 30 50 70 100 150 200 250 300 400 500 700 850 925 1000"
          }
        } elseif {$numlevs == 17} {
          if {$firstlev == 1000} {
            set levmap "1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 10"
          } else {
            set levmap "10 30 50 70 100 150 200 250 300 400 500 600 700 775 850 925 1000"
          }
        } elseif {$numlevs == 21} {
          if {$firstlev == 1000} {
            set levmap "1000 925 850 700 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1"
          } else {
            set levmap "1 2 3 5 7 10 20 30 50 70 100 150 200 250 300 400 500 700 850 925 1000"
          }
        } elseif {$numlevs == 23} {
          if {$firstlev == 1000} {
            set levmap "1000 925 850 775 700 600 500 400 300 250 200 150 100 70 50 30 20 10 7 5 3 2 1"
          } else {
            set levmap "1 2 3 5 7 10 20 30 50 70 100 150 200 250 300 400 500 600 700 775 850 925 1000"
          }
        }
        set level [lsearch -exact $levmap $level]
        set level [expr $level]
        set outlevels [lappend outlevels $level]
      } elseif {$levtype == 113} {
        # Potential temperature levels are 265,275,285,300,315,330,350,370,395,430,475,530,600,700,850
        lappend levmap 265 275 285 300 315 330 350 370 395 430 475 530 600 700 850
        set level [lsearch -exact $levmap $level]
        set level [expr $level-1]
        set outlevels [lappend outlevels $level]
      }
    }
  puts "Levels to output are: $levels.\n"
  set subset_z true
} else {
  puts "\nNo subsetting in z (vertical co-ordinate)...\n"
}

# Set output grid type
setgridtype regular

# Convert to Gaussian grid if input data is Spectral
regexp -- {grid code = ([0-9]+)} [printhead] dummy gridname
if {$gridname == 50} {
  puts "Spectral data so transforming to regular 1 x 1 degree (default) grid...\n"
  setgrid 1 360 0.0 1.0
  setgrid 2 181 90.0 -1.0
  spec_trans $fieldlist
}

# Test values for $xstart then subset
if {[info exists xstart]} {
  if {$xstart >= -360 && $xstart <= 360} {
    # Set the number of grid points in x
    set xnum [expr (($xend-$xstart)/$xinc)+1]
    # Convert $xnum into an integer
    set xnum [lindex [split $xnum .] 0]
    # Set output grid in x
    setgrid 1 $xnum $xstart $xinc
    set interp_xy true
    puts "\nSubsetting in x from $xstart deg to $xend deg at intervals of $xinc deg."
  }
} else {
  puts "Not subsetting in x...\n";
}

# Test for $ystart and if exists then subset
if {[info exists ystart]} {
  if {$ystart < $yend} {
    set ytemp $ystart
    set ystart $yend
    set yend $ytemp
  }
  if {$ystart >= -90 && $ystart <= 90} {
    # Set the number of grid points in y
    set ynum [expr (($ystart-$yend)/$yinc)+1]
    # Convert $ynum into an integer (ignore figures after decimal point)
    set ynum [lindex [split $ynum .] 0]
    # Set output grid in y
    set yinc [expr $yinc*(-1)]
    setgrid 2 $ynum $ystart $yinc
    set interp_xy true
    puts "\nSubsetting in y from $ystart deg to $yend deg at intervals of $yinc deg.\n"
  }
} else {
  puts "Not subsetting in y...\n";
}

# Interpolate to new grid if required
if {$interp_xy} {
  puts "Interpolating to new grid.\n"
  interp_grid $fieldlist
}

# Create zonal means if required
if {$zm} {
  puts "Calculating zonal means.\n"
  zonal_mean $fieldlist
}

# Create meridional means if required
if {$mm} {
  puts "Calculating meridional means.\n"
  meridonal_mean $fieldlist
}

# Now set the levels in the output [if not defined, set to all]
if {$subset_z} {
  setdim 3 $fieldlist "$outlevels"
}
#  Write out all input fields to a single Netcdf file

writefile $outformat $outfile $fieldlist


