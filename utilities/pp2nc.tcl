#!/users/modellers/pica/Software/bin/convsh

#  Convsh script pp2nc.tcl
#
#  Convert input files into single Netcdf file.
#
#  All input files must contain the same fields and have identical dimensions
#  except for the time dimension.  For example to convert UM output files into
#  a Netcdf file use the following command:
#
#      pp2nc.tcl -i xaavaa.pc* -o xaava.nc

#  Write out Netcdf file
set outformat netcdf

#  Automatically work out input file type
set filetype 0

#  Convert all fields in input files to Netcdf
set fieldlist -1

#  Get command line arguments:
#      -i input files (can be more than one file)
#      -o output file (single file only)

set i false
set o false
foreach arg $argv {
   switch -glob -- $arg {
      -i      {set i true ; set o false}
      -o      {set i false ; set o true}
      -*      {puts "unknown option $arg" ; set i false; set o false}
      default {
         if {$i} {
            set infile [lappend infile $arg]
         } elseif {$o} {
            set outfile [lappend outfile $arg]
         } else {
            puts "unknown option $arg"
         }
      }
   }
}

if {! [info exists infile]} {
   puts "input file name must be given"
   exit
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

#  Write out all input fields to a single Netcdf file

writefile $outformat $outfile $fieldlist
