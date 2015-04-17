#!/usr/bin/perl -w

###################################################
#  
# $Id: flens_to_nrlib.pl 534 2009-12-11 10:16:23Z perroe $
#
#
# Copyright (C) Norsk Regnesentral, 2009 
#               Per Roe
#               
# Conversion from FLENS CVS format to NRLib2 format.
# * Changes filenames .h => .hpp
# * Adds .tcc files to .h files
#
# Usage: ./flens_to_nrlib.pl [flens-dir] [nrlib-dir] 
#
###################################################

use strict;
use File::Spec;
use File::Basename;

if (@ARGV != 2) {
    # Too few arguments
    &usage;
}

my $flens_dir = $ARGV[0];
my $nrlib_dir = $ARGV[1];

# 1. Copy .cc files and rename to .cpp
#    Replace .h with .hpp in includes and use local path.
my $pattern = File::Spec->catfile($flens_dir, "*.cc");
my @cc_files = glob($pattern);
foreach my $file (@cc_files) {
    my $cpp_file = $file;
    $cpp_file =~ s/\.cc$/.cpp/;
    $cpp_file = File::Spec->catfile($nrlib_dir, basename($cpp_file));
    printf "%-25s => %-25s\n", basename($file), basename($cpp_file);
    open IN,  "<$file";
    open OUT, ">$cpp_file";
    while (<IN>) {
        # Handle includes of FLENS .h-files
        s/\<flens\/(\w+)\.h\>/"$1.hpp"/;
        print OUT $_;
    }
}

# 2. Copy .h files and rename to .hpp
#    Merge .tcc-files with corresponding .h-files.
$pattern = File::Spec->catfile($flens_dir, "*.h");
my @h_files = glob($pattern);
foreach my $file (@h_files) {
    my $hpp_file  = $file;
    $hpp_file =~ s/\.h$/.hpp/;
    $hpp_file = File::Spec->catfile($nrlib_dir, basename($hpp_file));
    printf "%-25s => %-25s\n", basename($file), basename($hpp_file);
    open IN,  "<$file";
    open OUT, ">$hpp_file";
    while (<IN>) {
        # Handle includes of FLENS .h-files
        s/\<flens\/(\w+)\.h\>/"$1.hpp"/;
        if (/\#include.*\<flens\/(\w+)\.tcc\>/) {
            # Copy contents of .tcc - file into .h-file
            print OUT "// -------- TEMPLATE IMPLEMENTATIONS ---------\n";
            print OUT "// Begin: Including file $1.tcc\n";
            my $tcc_file = File::Spec->catfile($flens_dir, "$1.tcc");
            printf "%-25s =>\n", basename($tcc_file);
            open TCC_FILE, "<$tcc_file";
            while (<TCC_FILE>) {
                s/\<flens\/(\w+)\.h\>/"$1.hpp"/;
                print OUT $_;
            }
            print OUT "// End: Including file $1.tcc\n";
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
