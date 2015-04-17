#!/usr/bin/perl -w

###################################################
#  
# $Id: nrlib_to_flens.pl 534 2009-12-11 10:16:23Z perroe $
#
#
# Copyright (C) Norsk Regnesentral, 2009 
#               Per Roe
#               
# Conversion from NRLib2 format to FLENS CVS format.
# * Changes filenames .hpp => .h
# * Changes filenames .cpp => .cc
# * Extracts template implementations from .h-files and put them
#   into seperate .tcc-files.
#
# Usage: ./flens_to_nrlib.pl [nrlib-dir] [flens-dir]
#
###################################################

use strict;
use File::Spec;
use File::Basename;

if (@ARGV != 2) {
  # Too few arguments
  &usage;
}

my $nrlib_dir = $ARGV[0];
my $flens_dir = $ARGV[1];

# 1. Copy .cpp files and rename to .cc
#    Replace .hpp with .h in includes and use <flens/foo.h> syntax.
my $pattern = File::Spec->catfile($nrlib_dir, "*.cpp");
my @cpp_files = glob($pattern);
foreach my $file (@cpp_files) {
  my $cc_file = $file;
  $cc_file =~ s/\.cpp$/.cc/;
  $cc_file = File::Spec->catfile($flens_dir, basename($cc_file));
  printf "%-25s => %-25s\n", basename($file), basename($cc_file);
  open IN,  "<$file";
  open OUT, ">$cc_file";
  while (<IN>) {
    # Handle includes of FLENS .h-files
    s/\"(\w+)\.hpp\"/<flens\/$1.h>/;
    print OUT $_;
  }
}

# 2. Copy .hpp files and rename to .h
#    Extract .tcc-files from corresponding .h-files.
$pattern = File::Spec->catfile($nrlib_dir, "*.hpp");
my @hpp_files = glob($pattern);
foreach my $file (@hpp_files) {
  my $h_file  = $file;
  $h_file =~ s/\.hpp$/.h/;
  $h_file = File::Spec->catfile($flens_dir, basename($h_file));
  printf "%-25s => %-25s\n", basename($file), basename($h_file);
  open IN,  "<$file";
  open OUT, ">$h_file";
  while (<IN>) {
    # Handle includes of FLENS .h-files
    s/\"(\w+)\.hpp\"/<flens\/$1.h>/;
    if ( /\/\/ -------- TEMPLATE IMPLEMENTATIONS ---------/ ) {
      if (defined($_ = <IN>)) {
        if ( /\/\/ Begin: Including file (\w+)\.tcc/ )  {
          my $tcc_file = File::Spec->catfile($flens_dir, "$1.tcc");
          printf "%-25s => %-25s\n", "", basename($tcc_file);
          open TCC_FILE, ">$tcc_file";
          my $finished = 0;
          while ( !$finished && (defined($_ = <IN>)) ) {
            s/\"(\w+)\.hpp\"/<flens\/$1.h>/;
            if ( /\/\/ End: Including file (\w+)\.tcc/ ) {
              print OUT "#include <flens/$1.tcc>\n";
              $finished = 1;
            }
            else {
              print TCC_FILE $_;
            }
          }
        }
        else {
          print "line = $_\n";
          die "Error parsing file.";
        }                
      }
      else {
        print "line = $_\n";
        die "Error parsing file.";
      }
    }
    else {
      print OUT $_;
    }
  }
}




###################################################################
#
#                           SUBROUTINES
#
###################################################################

# Print out usage

sub usage {
  print "\nusage: flens_to_nrlib.pl [flens-dir] [nrlib-dir]\n\n";

  exit;
}
